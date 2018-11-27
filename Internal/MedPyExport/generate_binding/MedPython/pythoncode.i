
%pythoncode %{

"""
Enable stderr capturing under ipynb
"""
class _stderr_fix:
  STDERR_FD = 2
  prevfd = None
  tmp = None
  reader = None
  is_in_ipynb = None

  def in_ipynb(self):
    try:
        if str(type(get_ipython())) == "<class 'ipykernel.zmqshell.ZMQInteractiveShell'>":
            return True
        else:
            return False
    except NameError:
        return False

  def __init__(self):
    from os import dup,dup2
    from tempfile import NamedTemporaryFile
    self.is_in_ipynb = self.in_ipynb()
    if self.is_in_ipynb:
      self.tmp = NamedTemporaryFile()
      self.prevfd = dup(self.STDERR_FD)
      dup2(self.tmp.fileno(), self.STDERR_FD)
      self.reader = open(self.tmp.name)

  def __del__(self):
    if self.is_in_ipynb:
      from os import dup2
      dup2(self.prevfd, self.STDERR_FD)

  def get_cerr(self):
    if self.is_in_ipynb:
      return self.reader.read()
    return ''


_stderr_fix_instance = _stderr_fix()


def cerr():
  return _stderr_fix_instance.get_cerr()


"""
Enable iterators on vector and map adaptor objects
"""


class MapAdaptorKeyIter:
    def __init__(self,o):
        self.obj = o.keys()
        self.i = 0
        self.prev_i = None
        self.max_i = len(self.obj)
    def next(self):
        if self.i >= self.max_i:
            raise StopIteration
        else:
            self.prev_i, self.i = self.i, self.i+1
            return self.obj[self.prev_i]


class IntIndexIter:
    def __init__(self,o):
        self.obj = o
        self.i = 0
        self.prev_i = None
        self.max_i = len(o)
    def next(self):
        if self.i >= self.max_i:
            raise StopIteration
        else:
            self.prev_i, self.i = self.i, self.i+1
            return self.obj[self.prev_i]

def ___fix_vecmap_iter():
    from inspect import isclass
    glob = globals()
    for i in glob:
        o = glob[i]
        try:
          if (isclass(o) and '__len__' in dir(o) and '__getitem__' in dir(o) and not '__iter__' in dir(o) 
              and i.endswith('VectorAdaptor')) :
              setattr(o, '__iter__', lambda x: IntIndexIter(x))
          elif (isclass(o) and '__getitem__' in dir(o) and 'keys' in dir(o) and not '__iter__' in dir(o) 
              and i.endswith('MapAdaptor')) :
              setattr(o, '__iter__', lambda x: MapAdaptorKeyIter(x))
        except: pass

___fix_vecmap_iter()



"""
External Methods in addition to api
"""

def __export_to_pandas(self, sig_name_str, translate=True, pids=None):
    """get_sig(signame [, translate=True][, pids=None]) -> Pandas DataFrame
         translate : If True, will decode categorical fields into a readable representation in Pandas
         pid : If list is provided, will load only pids from the given list
               If 'All' is provided, will use all available pids
    """
    import pandas as pd
    use_all_pids = 0
    if isinstance(pids, str) and pids.upper()=='ALL':
      use_all_pids = 1
      pids=list()
    if pids is None: pids=list()
    sigexporter = self.export_to_numpy(sig_name_str, pids, use_all_pids)
    df = pd.DataFrame.from_dict(dict(sigexporter))
    if not translate:
      return df
    for field in sigexporter.get_categorical_fields():
        df[field] = df[field].astype('category').cat.rename_categories(dict(sigexporter.get_categorical_field_dict(field)))
    return df


def __bind_external_methods():
    setattr(globals()['PidRepository'],'get_sig', __export_to_pandas)

__bind_external_methods()

%}
