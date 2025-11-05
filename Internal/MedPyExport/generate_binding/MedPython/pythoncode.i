
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
def __export_to_pandas(self, sig_name_str:str, translate:bool=True, pids:list[int]|str|None=None, float32to64:bool=True, free_signal:bool=True, regex_str:str|None = None, regex_filter:str = '') -> 'pd.DataFrame':
    """get_sig(signame [, translate=True][, pids=None, float32to64=True][,regex_str = None]) -> Pandas DataFrame
         translate : If True, will decode categorical fields into a readable representation in Pandas
         pid : If list is provided, will load only pids from the given list
               If 'All' is provided, will use all available pids
         float32to64 : If True, will convert all float32 columns to float64
         free_signal : If True, will free signal memory as soon as it is loaded into export arrays
         regex_str :  string (if string, defualt column is 'val') or dictionary between column name (e.g, 'val') to regex string to filter
         regex_filter :  string to filter categories paretns this match regex
    """
    import pandas as pd
    import numpy as np
    use_all_pids = 0
    if isinstance(pids, str) and pids.upper()=='ALL':
      use_all_pids = 1
      pids=list()
    if pids is None: pids=list()
    df = self.export_to_numpy(sig_name_str, pids, use_all_pids, int(translate), int(free_signal), regex_filter)
    

    cat_dict = dict()
    cat_dict_int = dict()
    if translate:
      for fld in df.get_categorical_fields():
        cat_dict[fld] = df.get_categorical_field_dict(fld)      
        cat_dict_int[fld] = df.get_categorical_field_dict_int(fld)
    df = dict(df)
    df = pd.DataFrame.from_dict(df)

    if not (regex_str is None):
      dict_id = self.dict_section_id(sig_name_str)
      if (dict_id == 0):
        print("Invalid sig name - skip filter")
      else:
        if not(type(regex_str) == dict):
          regex_str = {'val' : regex_str}
        
        for fld in regex_str:
          lut = self.get_lut_from_regex(dict_id, regex_str[fld])
          lut_np = np.array(lut)
          if translate:
            lut_np = lut_np[np.array(cat_dict_int[fld])]
          tt = df[fld].astype(int).values
          df = df[lut_np[tt - 1]].reset_index(drop=True) # values start from 1 
    if float32to64:
      for column_name in df:
        if df[column_name].dtype == np.float32:
           df[column_name] = df[column_name].astype(np.float64)
    if not translate:
      return df
    for fld in cat_dict:
        df[fld] = pd.Categorical.from_codes(codes=df[fld],categories=cat_dict[fld])
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

def __bootstrapResult_to_df(self):
    import pandas as pd
    dict_obj={'Cohort' : [], 'Measurement': [], 'Value': []}
    for cohort in self.keys():
        cohort_res=self[cohort]
        for measure in cohort_res.keys():
            dict_obj['Cohort'].append(cohort)
            dict_obj['Measurement'].append(measure)
            dict_obj['Value'].append(cohort_res[measure])
    df=pd.DataFrame(dict_obj)
    return df

def convert_to_bootstrap_input(x, arg_name=None):
    import pandas as pd
    import numpy as np
    if type(x) is list:
        x=np.array(x)
    if type(x) is np.ndarray:
        x=x.astype(float)
    if type(x) is pd.Series:
        x=x.astype(float).to_numpy()
    if np.isnan(x).sum()>0:
        if arg_name is None:
            raise NameError('Error - input array has nan inside')
        else:
            raise NameError('Error - input array %s has nan inside'%(arg_name))
    return x

def __bootstrap_wrapper(self, preds, labels):
    import pandas as pd
    import numpy as np
    preds=convert_to_bootstrap_input(preds, 'preds')
    labels=convert_to_bootstrap_input(labels, 'labels')
    res = self._bootstrap(preds, labels)
    return res

def __bootstrap_pid_wrapper(self, pids, preds, labels):
    import pandas as pd
    import numpy as np
    pids=convert_to_bootstrap_input(pids, 'pids')
    preds=convert_to_bootstrap_input(preds, 'preds')
    labels=convert_to_bootstrap_input(labels, 'labels')
    res = self._bootstrap_pid(pids, preds, labels)
    return res

def __btsimple_to_df(self):
    import pandas as pd
    dict_obj={'Measurement': [], 'Value': []}
    for measure in self.keys():
        dict_obj['Measurement'].append(measure)
        dict_obj['Value'].append(self[measure])
    df=pd.DataFrame(dict_obj)
    return df

def __bind_external_methods():
    setattr(globals()['PidRepository'],'get_sig', __export_to_pandas)
    setattr(globals()['Features'],'to_df', __features__to_df_imp)
    setattr(globals()['Features'],'from_df', __features__from_df_imp)
    setattr(globals()['StringBtResultMap'],'to_df', __bootstrapResult_to_df)
    setattr(globals()['StringFloatMapAdaptor'],'to_df', __btsimple_to_df)
    setattr(globals()['Bootstrap'],'bootstrap', __bootstrap_wrapper)
    setattr(globals()['Bootstrap'],'bootstrap_pid', __bootstrap_pid_wrapper)

__bind_external_methods()

"""
Remove SWIG's global variable access which makes issues for reflection actions
#if 'cvar' in __dict__: del cvar
"""
try:
    del cvar
except: pass

Global = GlobalClass()

%}
