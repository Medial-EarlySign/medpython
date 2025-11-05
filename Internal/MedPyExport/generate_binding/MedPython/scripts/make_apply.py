#!/usr/bin/env python

from __future__ import print_function
import re
import sys

"""
some constants
"""

parser = re.compile(
    r"MEDPY_NP_(?P<iotype>INPUT|OUTPUT|INPLACE|VARIANT_OUTPUT)\(\s*(?P<params>.*?)\s*\)"
)
parser2 = re.compile(r"MEDPY_(?P<proptype>GET|SET)_(?P<propname>[A-Z0-9a-z_]+)")
parser3 = re.compile(r"^\s*class\s+(?P<classname>[A-Za-z0-9]+)")
parser4 = re.compile(r"^\s*MEDPY_DOC_Dyn\s*\(\"(?P<classname>[^\"]+\")")
first_word = re.compile(r"^[a-zA-Z_]*")
io_conv_tbl = {
    "INPUT": "%apply ({at}* IN_ARRAY1, unsigned long long DIM1) {p}",
    "OUTPUT": "%apply ({at}** ARGOUTVIEWM_ARRAY1, unsigned long long* DIM1) {p}",
    "VARIANT_OUTPUT": "%apply (void** ARGOUTMVAR_ARRAY1, unsigned long long* DIM1, int* NPYDTC1) {p}",
    "INPLACE": "%apply ({at}* INPLACE_ARRAY1, unsigned long long DIM1) {p}",
}


"""
classes already renamed in this set
"""
renamed_class = set()

"""
Parse arguments and open files
"""


def main(args):
    from os.path import (
        isfile,
        isdir,
        splitext,
        join as joinpath,
        abspath,
        exists as pathexists,
    )
    from os import walk

    if len(args) <= 2:
        print(
            args[0]
            + " - A tool to generate function definitions for numpy.i using a macro in the code"
        )
        print(
            "usage: " + args[0] + " [input_file_name|input_dir_name] [output_file_name]"
        )
        return
    in_fname, out_fname = args[1:3]
    if isfile(in_fname):
        print("parsing " + in_fname)
        file_parse(in_file=open(in_fname), out_file=open(out_fname, "w"))
    elif isdir(in_fname):
        _, _, filenames = next(walk(in_fname), (None, None, []))
        outf = open(out_fname, "w")
        for fname in filenames:
            fname = abspath(joinpath(in_fname, fname))
            _, ext = splitext(fname)
            if ext != ".h":
                continue
            print("parsing " + fname)
            file_parse(in_file=open(fname), out_file=outf)
    elif not pathexists(in_fname):
        print("cannot find " + in_fname)
    else:
        print("undefined action for " + in_fname)


"""
Parse a single input file
"""


def file_parse(in_file, out_file):
    apply_list = {}
    prop_list = {}
    cur_class = ""
    for line in in_file:
        if "class" in line:
            m = parser3.search(line)
            if m:
                d = m.groupdict()
                cur_class = d["classname"]
                if cur_class not in renamed_class and cur_class.startswith("MP"):
                    print(
                        "%rename({}) {};".format(cur_class[2:], cur_class),
                        file=out_file,
                    )
                    renamed_class.add(cur_class)
        if "MEDPY_NP_" in line and not line.startswith("#define"):
            mlst = parser.finditer(line)
            if mlst:
                for m in mlst:
                    d = m.groupdict()
                    array_type = first_word.search(d["params"]).group(0)
                    key = (d["iotype"], array_type)
                    apply_list.setdefault(key, []).append(d["params"])
        if "MEDPY_DOC_Dyn" in line and not line.startswith("#define"):
            # Ignore dynamic docstrings for now
            doc_text = parser4.findall(line)
            if len(doc_text) == 0:
                continue
            doc_text = doc_text[0].strip().strip('"')
            prop_list.setdefault(cur_class, []).append(("DOC", doc_text))
        if ("MEDPY_GET_" in line or "MEDPY_SET_" in line) and not line.startswith(
            "#define"
        ):
            m = parser2.search(line)
            if m and cur_class:
                d = m.groupdict()
                prop_type = first_word.search(d["proptype"]).group(0)
                prop_name = d["propname"]
                prop_list.setdefault(cur_class, []).append((prop_type, prop_name))
    for (io_type, array_type), params in apply_list.items():
        for i in range(len(params)):
            params[i] = "({})".format(params[i])
        params_str = "{" + ", ".join(params) + "}"
        print(io_conv_tbl[io_type].format(at=array_type, p=params_str), file=out_file)
    for cur_class in prop_list:
        print("%extend {}{{\n".format(cur_class), file=out_file)
        for ptype, doc_text in prop_list[cur_class]:
            if ptype != "DOC":
                continue
            print("{} {{}}".format(doc_text), file=out_file)
        print("%pythoncode %{{".format(cur_class), file=out_file)
        prop_def = {}
        for ptype, pname in prop_list[cur_class]:
            if ptype == "DOC":
                continue
            prop_def[pname] = ["None", "None"]
        for ptype, pname in prop_list[cur_class]:
            if ptype == "DOC":
                continue
            # if ptype.lower()=='get':
            #  continue
            prop_func_name = "MEDPY_" + ptype + "_" + pname
            # print('    __swig_{}methods__["{}"] = {}'.format(ptype.lower(), pname, prop_func_name), file=out_file)
            prop_def[pname][0 if ptype == "GET" else 1] = prop_func_name
        # print('    if _newclass:', file=out_file)
        for pname in prop_def:
            print(
                "      {} = property("
                "{}"
                ","
                "{}"
                ")".format(pname, prop_def[pname][0], prop_def[pname][1]),
                file=out_file,
            )
        print("  %}\n};", file=out_file)
    


if __name__ == "__main__":
    main(sys.argv)
