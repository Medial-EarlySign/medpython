#!/usr/bin/env python

from __future__ import print_function
import re
import sys

"""
some constants
"""

def fix_docstr(dstr):
  newdoc = re.sub(r'(\w+ )?self ?,? ?','',dstr)
  newdoc = re.sub(r'(std::)?vector< (std::)?string > (& )?','list_str_',newdoc)
  newdoc = re.sub(r'(, |\()(std::)?vector< (float|int|double) > (& )?','\\1list_\\3_',newdoc)
  newdoc = re.sub(r'(String|Int|Float|Double)Vector','list_\\1',newdoc)
  newdoc = re.sub(r'(, |\()(float|int|double) \* ','\\1numpy_arr_\\2_',newdoc)
  newdoc = re.sub(r'(, |\()(float|int|double) ','\\1\\2_',newdoc)
  newdoc = re.sub(r'(, |\()Features features','\\1features',newdoc)
  newdoc = re.sub(r'(, |\()PidRepository rep','\\1repository',newdoc)
  newdoc = re.sub(r'(, |\()Samples samples','\\1samples',newdoc)
  newdoc = re.sub(r'-> (std::)string','-> string',newdoc)
  newdoc = re.sub(r'(const )?(std::)?string (const )?(& )?','str_',newdoc)
  newdoc = re.sub(r'MEDPY_GET_(\w+)\(\)','\\1 ; property(read)',newdoc)
  newdoc = re.sub(r'MEDPY_SET_(\w+)\((.*?)\)','\\1 <- \\2 ; property(write)',newdoc)
  return newdoc


"""
Parse arguments and open files
"""
def main(args):
    lines = file(args[1]).readlines()
    for i in range(len(lines)):
        line = lines[i]
        if line.lstrip().startswith('"""') and line.rstrip().endswith('"""') and line.strip() != '"""':
            #print('\nold:',line.strip())
            #print('new:',fix_docstr(line).strip())
            lines[i] = fix_docstr(line)
            
    file(args[1], 'w').writelines(lines)

if __name__ == "__main__":
  main(sys.argv)
