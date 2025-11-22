import os
import sys
import subprocess
import shutil
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext


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

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))

        # This is where we want the .so and .py files to end up
        # We target the 'med' folder inside the build directory
        destination_dir = os.path.join(extdir, "med")

        # Ensure the destination exists
        if not os.path.exists(destination_dir):
            os.makedirs(destination_dir)

        # Config Arguments
        cmake_args = [
            f"-DCMAKE_POLICY_VERSION_MINIMUM=3.10",
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={destination_dir}",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            f"-DCMAKE_BUILD_TYPE=Release",
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

        # 3. MANUAL COPY (If CMake doesn't put the .py file exactly where we want)
        # CMake will put _medpython.so in destination_dir (because of CMAKE_LIBRARY_OUTPUT_DIRECTORY)
        # But SWIG generates 'medpython.py' which might stay in the temp folder depending on your CMakeLists.
        # Let's find medpython.py and move it to destination_dir/med/

        # Look for the generated python file in the build temp tree
        for root, dirs, files in os.walk(self.build_temp):
            if "medpython.py" in files:
                shutil.copy(os.path.join(root, "medpython.py"), destination_dir)
                break
        # Generate med.py:
        with open(os.path.join(destination_dir, "med.py"), "w") as f:
            f.write(
                "from medpython import * ; import medpython as _med ; __doc__=_med.__doc__ ; __all__=_med.__all__ ;\n"
            )
        # strip _medpython.so


setup(
    name="medpython",
    version="1.0.0",
    # Find packages in 'src' (lib1, lib2, and med)
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    # We define one extension, pointing to the root where CMakeLists.txt is
    ext_modules=[CMakeExtension("med/medpython", sourcedir=".")],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
)
