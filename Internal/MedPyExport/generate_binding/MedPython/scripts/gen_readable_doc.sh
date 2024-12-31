#!/bin/bash
echo '%define MPDOCSTRING'
echo '"'
python -c "import medpython; help(medpython)" | egrep -v '(__from_df_adaptor|_swig_|__getattr__|__setattr__|__del__ lambda self|__dict__|__weakref__|_o_get|dictionary for instance variables|list of weak references to|^     \\| $|Proxy of C++|self, \*args\)$|\|  __init__\(self\)$|\|      delete_|     \|  [a-zA-Z_0-9]+\(self)' | sed 's/MEDPY_GET_\(\w\+\)[^-]\+\(.*\)$/\1 \2   (property getter)/;s/MEDPY_SET_\(\w*\).* self, \(\w\+\).*$/\1 <- \2   (property setter)/;s/\(\s\+class \w\+\).*$/\n\n\n\1/;s/ defined here:$/:/;s/\s\+|  -\+$//;s/     |  /        /;/staticmethod/,+11d;/^DATA$/,+100d;0,/^CLASSES$/d;0,/Time/d;s/(\w\+ self,* */(/;s/MEDPY__to_df/to_df/;s/MEDPY__from_df/from_df/;/__to_df_imp/d;/__from_df_imp/d;' |  tr '\n' '\r' | sed 's/\r\s*Data descriptors:\r\r/\r/g;s/\r\s*Data and other attributes:\r\r/\r/g;s/^\r*//;s/\r\s*\(\w\+\)\r\(\s*\)\1 -> /\r\2\1 -> /g;s/\r\s*\(\w\+\)\r\(\s*\)\1    (property/\r\2\1    (property/g' |  tr '\r' '\n'
echo '"'
echo '%enddef'
