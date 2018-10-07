
# coding: utf-8

# In[1]:

import sys
import numpy as np
import pandas as pd
#sys.path.insert(0,'/nas1/UsersData/shlomi/MR/Libs/Internal/MedPyExport/generate_binding/CMakeBuild/Linux/Release/MedPython')

import medpython as med

rep = med.PidRepository()
rep.read_all('/home/Repositories/THIN/thin_jun2017/thin.repository',[],['WBC','Cancer_Location',
                                                                          'GENDER','BYEAR',
                                                                          'Cancer_Location',
                                                                          'DEATH','ENDDATE',
                                                                          'Albumin','BDATE'])

print med.cerr()




# In[13]:

def fix_ts_m1900(df, col):
    import datetime as dt
    df[col] = dt.datetime(1900,1,1)+pd.TimedeltaIndex(df[col], unit='m')

def fix_type(df, col, newtype): df[[col]] = df[[col]].astype(newtype, copy=False)

def fix_date_ymd(df, col): df[col] = pd.to_datetime(df[col], format='%Y%m%d')

def fix_name(df, old_col, new_col): df.rename(columns={old_col: new_col}, inplace=True)
    
albumin = rep.get_sig('Albumin')
albumin = albumin[albumin['date'] % 100 != 0]
fix_name(albumin,'val','Albumin')
fix_date_ymd(albumin, 'date')

albumin.head(10)


# In[3]:

gender = rep.get_sig('GENDER')
fix_type(gender,'val', int)
fix_name(gender,'val','Gender')
gender.head(5)


# In[4]:

bdate = rep.get_sig('BDATE')
fix_type(bdate,'val', int)
bdate['val'] = bdate['val'] + 1
fix_date_ymd(bdate, 'val')
fix_name(bdate,'val','BDate')
bdate.head(5)


# In[5]:

mortality = rep.get_sig('DEATH')
fix_type(mortality,'val', int)
mortality = mortality[(mortality['val'] % 100 <= 31) & (mortality['val'] % 100 > 0)]
fix_date_ymd(mortality, 'val')
fix_name(mortality,'val','MortDate')
mortality.head(10)


# In[6]:

data = reduce(lambda left,right: pd.merge(left, right, on='pid', how="left", sort=False),[albumin,gender,bdate,mortality])
data.head()


# In[7]:

td2year = pd.Timedelta(days=365*2)
data['age'] = data['date'].dt.year - data['BDate'].dt.year
data['mort_2y'] = data['MortDate']-data['date'] <= td2year
data.head(20)


# In[8]:

data.dtypes


# In[9]:

old_len = len(data)
data = data[data['Albumin']<8]
data = data[data['age']>0]
print "{} record filtered out".format(old_len - len(data))


# In[10]:

get_ipython().magic(u'matplotlib inline')
import matplotlib.pyplot as plt
plt.plot(data.groupby('age')['Albumin'].mean())


# In[11]:

data2 = data[data['mort_2y'] == True].groupby('Albumin').count()
data2 = data2.reset_index()
data2.head()


# In[12]:

plt.plot(data2['Albumin'],data2['mort_2y'])


# In[ ]:




