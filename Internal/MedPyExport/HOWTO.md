# User Guide


## Note: stderr in Jupyter notebooks
Our API sometimes prints messages to stderr which are not displayed in jupyter notebooks. If you execute the python script in command line you will see the massages however in Jupyter use the cerr() utility function to flush the buffer and print it in the notebook's output window:
```
med.cerr()
```

## getting started

### import the extension:
```python
import medpython as med
```
### Load the repository:
```python
rep = med.PidRepository()
rep.read_all('/home/Repositories/THIN/thin_jun2017/thin.repository',[],['GENDER','BYEAR','DEATH','BDATE'])
print med.cerr()
```

## Use utility function get_sig() to get a signal as Pandas DataFrame

```
albumin = rep.get_sig('Albumin') 
albumin.head(5)
```

### You might want to clean around and change stuff:
The following functions might be useful:

```python
import pandas as pd

def fix_ts_m1900(df, col):
    """ change minutes timestamp counter since 1900 to a pandas datetime object"""
    import datetime as dt
    df[col] = dt.datetime(1900,1,1)+pd.TimedeltaIndex(df[col], unit='m')

def fix_type(df, col, newtype):
    """change type of a given column"""
    df[[col]] = df[[col]].astype(newtype, copy=False)

def fix_date_ymd(df, col): 
    """parse a YYYYMMDD integer/string format to a pandas datetime object"""
    df[col] = pd.to_datetime(df[col], format='%Y%m%d')

def fix_name(df, old_col, new_col): 
    """rename a column"""
    df.rename(columns={old_col: new_col}, inplace=True)

```

These could be used like this:

```python
bdate = rep.get_sig('BDATE')
fix_type(bdate,'val', int)
bdate['val'] = bdate['val'] + 1
fix_date_ymd(bdate, 'val')
fix_name(bdate,'val','BDate')
bdate.head(5)
```

### Inspect the available function and get help (Documentation is work in progress)

```python
help(med.PidRepository)
```

# Addional examples
* [Example of Test+Train model][1] - A partial translation of Flow's TestTrain class to Python
* [Iterating over a Signal][2] - This is not the proper way to work with Python, since Python loops are very slow, however, it's nice as a poc. 

[1]: http://node-02:8888/notebooks/Shlomi/NewAPI/learn.ipynb
[2]: http://node-02:8888/notebooks/Shlomi/NewAPI/USig_Iterate.ipynb

