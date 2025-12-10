import os
import platform
import sys
import subprocess
import shutil
from setuptools import setup, Extension, find_packages
from setuptools.command.sdist import sdist as _sdist
from setuptools.command.build_ext import build_ext


# --- PATH CONFIGURATION ---
SETUP_DIR = os.path.dirname(os.path.abspath(__file__))
# The "Real" source (Parent directory)
EXTERNAL_SOURCES = [
    os.path.abspath(os.path.join(SETUP_DIR, "..", "..", "..", "External")),
    os.path.abspath(os.path.join(SETUP_DIR, "..", "..", "..", "Internal")),
]
# The "Internal" destination (Inside the package)
INTERNAL_SOURCE = os.path.join(SETUP_DIR, "src", "med", "cpp_code")
BUNDLED_CMAKE_DIR = os.path.join(
    INTERNAL_SOURCE, "Internal", "MedPyExport", "generate_binding"
)

SKIP_FOLDER_NAMES = set(
    [
        ".git",
        "CMakeBuild",
        "Release",
        "build",
        "dist",
        "wheelhouse",
        "__pycache__",
        ".vscode",
    ]
)


class CustomSdist(_sdist):
    """Custom sdist command to copy external C++ source into the package tree."""

    def run(self):
        # 1. Clean up any stale copies
        if os.path.exists(INTERNAL_SOURCE):
            shutil.rmtree(INTERNAL_SOURCE)

        def ignore_recursion(src: str, names: list[str]) -> list[str]:
            ignored = []

            # If we are currently looking at the Root folder (MedPyExport)
            ignored = list(set(names).intersection(SKIP_FOLDER_NAMES))
            if os.path.basename(src) == os.path.basename(SETUP_DIR):
                if "src" in names:
                    ignored.append("src")
            for name in names:
                if name.endswith(".egg-info"):
                    ignored.append(name)
            return ignored

        # 2. Copy the external folder into src/med/cpp_src
        for EXTERNAL_SOURCE in EXTERNAL_SOURCES:
            folder_name = os.path.basename(
                EXTERNAL_SOURCE
            )  # e.g., "Internal" or "External"
            destination = os.path.join(INTERNAL_SOURCE, folder_name)
            print(f"Copying C++ source from {EXTERNAL_SOURCE} -> {destination}")
            # dirs_exist_ok=True allows overwriting if needed (Python 3.8+)
            shutil.copytree(
                EXTERNAL_SOURCE,
                destination,
                dirs_exist_ok=True,
                ignore=ignore_recursion,
            )

        # 3. Run the standard sdist process (which looks at MANIFEST.in)
        super().run()

        # 4. (Optional) Cleanup: Delete the copy so your source tree stays clean
        # If you comment this out, you can verify the folder structure manually
        shutil.rmtree(INTERNAL_SOURCE)


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            subprocess.check_output(["cmake", "--version"])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: "
                + ", ".join(e.name for e in self.extensions)
            )

        for ext in self.extensions:
            self.build_extension(ext)

        self.cleanup_cpp_source()

    def cleanup_cpp_source(self):
        # build_lib is where setuptools assembles the package before zipping
        build_dir = os.path.abspath(self.build_lib)
        cpp_code_dest = os.path.join(build_dir, "med", "cpp_code")

        if os.path.exists(cpp_code_dest):
            print(f"Removing C++ source from wheel: {cpp_code_dest}")
            shutil.rmtree(cpp_code_dest)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        # Ensure the destination exists
        if not os.path.exists(extdir):
            os.makedirs(extdir)

        # Config Arguments
        cmake_args = [
            f"-DCMAKE_POLICY_VERSION_MINIMUM=3.10",
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}",
            f"-DPython3_EXECUTABLE={sys.executable}",
            f"-DCMAKE_BUILD_TYPE=Release",
            "-DDMLC_FORCE_SHARED_CRT=ON",
            # Add your specific flags here if needed
            # '-DSWIG_EXECUTABLE=...' (Usually auto-detected correctly in pyproject environment)
        ]

        build_args = ["--config", "Release"]

        # Parallel build
        if "CMAKE_BUILD_PARALLEL_LEVEL" not in os.environ:
            build_args += [f"-j{os.cpu_count()}"]

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        # 1. RUN CMAKE
        subprocess.check_call(
            ["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp
        )

        # 2. RUN BUILD
        subprocess.check_call(
            ["cmake", "--build", "."] + build_args, cwd=self.build_temp
        )
        print("########## Done CMAKE BUILD ############3")

        # 3. MANUAL COPY (If CMake doesn't put the .py file exactly where we want)
        # CMake will put _medpython.so in destination_dir (because of CMAKE_LIBRARY_OUTPUT_DIRECTORY)
        # But SWIG generates 'medpython.py' which might stay in the temp folder depending on your CMakeLists.
        # Let's find medpython.py and move it to destination_dir/med/

        # Look for the generated python file in the build temp tree
        for root, dirs, files in os.walk(self.build_temp):
            if "medpython.py" in files:
                shutil.copy(os.path.join(root, "medpython.py"), extdir)
                break
        # Generate med.py:
        with open(os.path.join(extdir, "med.py"), "w") as f:
            f.write(
                "from .medpython import *\nfrom . import medpython as _med\n__doc__=_med.__doc__\n__all__=_med.__all__ ;\n"
            )
        with open(os.path.join(extdir, "__init__.py"), "w") as f:
            f.write("from .med import *\n")
        # Remove lightgbm
        if os.path.exists(os.path.join(extdir, "lib_lightgbm.so")):
            os.remove(os.path.join(extdir, "lib_lightgbm.so"))
        # strip _medpython.so
        if platform.system() == "Linux" and shutil.which("strip"):
            if os.path.exists(os.path.join(extdir, "_medpython.so")):
                subprocess.check_call(
                    ["strip", os.path.join(extdir, "_medpython.so")],
                    cwd=self.build_temp,
                )
        winlib_file = os.path.join(extdir, "Release", "_medpython.pyd")
        if os.path.exists(winlib_file):
            shutil.copy(winlib_file, extdir)


if os.path.exists(BUNDLED_CMAKE_DIR):
    SOURCE_DIR = BUNDLED_CMAKE_DIR
    print(f"Building from SDIST source: {SOURCE_DIR}")
else:
    SOURCE_DIR = "."
    print(f"Building from DEV source: {SOURCE_DIR}")

setup(
    # Find packages in 'src' (lib1, lib2, and med)
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    # We define one extension, pointing to the root where CMakeLists.txt is
    ext_modules=[CMakeExtension("med.medpython", sourcedir=SOURCE_DIR)],
    cmdclass={
        "build_ext": CMakeBuild,
        "sdist": CustomSdist,
    },
    zip_safe=False,
    include_package_data=True,
)
