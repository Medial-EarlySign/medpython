
%pythoncode %{

__all__ = list(_medpython.get_public_objects())

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
    def __next__(self):
        if self.i >= self.max_i:
            raise StopIteration
        else:
            self.prev_i, self.i = self.i, self.i+1
            return self.obj[self.prev_i]
    def next(self):
        return self.__next__()


class IntIndexIter:
    def __init__(self,o):
        self.obj = o
        self.i = 0
        self.prev_i = None
        self.max_i = len(o)
    def __next__(self):
        if self.i >= self.max_i:
            raise StopIteration
        else:
            self.prev_i, self.i = self.i, self.i+1
            return self.obj[self.prev_i]
    def next(self):
        return self.__next__()

def __to_df_imp(self):
    import pandas as pd
    return pd.DataFrame.from_dict(dict(self.MEDPY__to_df()))

def __from_df_imp(self, df):
    import re
    adaptor = self.MEDPY__from_df_adaptor()
    type_requirements = dict(adaptor.type_requirements)
    for col_name in df.columns: 
        for col_req_name in type_requirements: 
            if re.match('^'+col_req_name+'$',col_name):
                if str(df[col_name].dtype) != type_requirements[col_req_name]:
                    df[col_name] = df[col_name].astype(type_requirements[col_req_name],copy=False)
        adaptor.import_column(col_name ,df[col_name].values)
    self.MEDPY__from_df(adaptor)


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
          if (isclass(o) and 'MEDPY__from_df' in dir(o) and 'MEDPY__from_df_adaptor' in dir(o)):
              setattr(o, 'from_df', __from_df_imp)
          if (isclass(o) and 'MEDPY__to_df' in dir(o)):
              setattr(o, 'to_df', __to_df_imp)
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

def __features__to_df_imp(self):
    import pandas as pd
    featMatFull = Mat()
    self.get_as_matrix(featMatFull)
    featMatFullNew = featMatFull.get_numpy_view_unsafe()
    col_names = self.get_feature_names()
    dfFeatures2 = pd.DataFrame(data = featMatFullNew, columns = col_names )
    
    samps = Samples()
    self.get_samples(samps)
    samps_df = samps.to_df()
    out = pd.concat([samps_df,dfFeatures2], axis=1, copy=False)
    return out

def __features__from_df_imp(self, features_df):
    # Dataframe to MedFeatures:
    
    mat = Mat()
    samples = Samples() 
    ind_sampes = features_df.columns.str.contains('pred_\\d+') | features_df.columns.isin(['id', 'split', 'time', 'outcome', 'outcomeTime']) 
    featuresNames = features_df.columns[~(ind_sampes)]
    # Build data matrix
    mat.set_signals(StringVector(list(featuresNames)))
    mat.load_numpy(features_df.loc[:,features_df.columns[~(ind_sampes)]].values)
    self.set_as_matrix(mat)
    # append Samples
    samples.from_df(features_df.loc[:,features_df.columns[ind_sampes]])
    self.append_samples(samples)



def __bind_external_methods():
    setattr(globals()['PidRepository'],'get_sig', __export_to_pandas)
    setattr(globals()['Features'],'to_df', __features__to_df_imp)
    setattr(globals()['Features'],'from_df', __features__from_df_imp)

__bind_external_methods()

"""
Remove SWIG's global variable access which makes issues for reflection actions
#if 'cvar' in __dict__: del cvar
"""
try:
    del cvar
except: pass

%}
