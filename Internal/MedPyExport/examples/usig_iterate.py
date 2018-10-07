
# coding: utf-8

# ## Doan't do this, using python loops is slow!

# In[1]:

import sys
#sys.path.insert(0,'/nas1/UsersData/shlomi/MR/Libs/Internal/MedPyExport/generate_binding/CMakeBuild/Linux/Release/MedPython')
import medpython as med

signame = 'Albumin'

rep = med.PidRepository()
rep.read_all('/home/Repositories/THIN/thin_jun2017/thin.repository',[],[signame])

print med.cerr()

for pid in rep.pids[20:30]:
  usv = rep.uget(pid, rep.sig_id(signame))
  if len(usv)>0:
    print('\n\n')
    for rec in usv:
      print("Patient {} had {}={:.2} at {}".format(pid, signame, rec.val(), rec.date()))


# In[ ]:

