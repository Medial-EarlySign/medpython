#!/bin/env python3

# Apply a model    

import med
import sys
import argparse
    
default_repository = '/home/Repositories/THIN/thin_jun2017/thin.repository'

def eprint(*args,**kwargs):
  kwargs['file']=sys.stderr
  print(*args,**kwargs)

ap = argparse.ArgumentParser()
ap.add_argument("--config",dest="configFile",type=str, default=default_repository,
  help=f"repository config file\n(default:{default_repository})")
ap.add_argument("--samples",dest='samplesFile',type=str, required=True,
  help="test samples file")
ap.add_argument("--inModel",type=str, required=True,
  help="model file")
ap.add_argument("--outSamples",type=str, required=True,
  help="Sample Ootput file")
ap.add_argument("--outPerformance",type=str, default='',
  help="optional performance file")
ap.add_argument("--preProcessors",type=str, default='',
  help="optional pre-processors json file")
ap.add_argument("--inpatient", default=False, action='store_true',
  help="indicate that relevant data is in-patient data")
p = ap.parse_args()

med.Global.RNG.srand(12131415)

if p.inpatient:
  eprint("Working on in-patient data");
  med.Global.default_time_unit = med.Time.Minutes
  med.Global.default_windows_time_unit = med.Time.Minutes


# Read Model 
eprint("Reading model")
model = med.Model()
model.read_from_file(p.inModel);

# Add pre-processors
if p.preProcessors != "":
  model.add_pre_processors_json_string_to_model("", p.preProcessors)
    
eprint("Getting required signals")
signalNamesSet = model.get_required_signal_names()

eprint("Reading samples")
samples = med.Samples()
try: samples.read_from_file(p.samplesFile)
except: 
  eprint("Cannot read samples from file %s", samplesFile.c_str());
  raise

eprint("Getting ids");
ids = samples.get_ids()
    
eprint("Reading Repository");
rep = med.PidRepository()
try: rep.read_all(p.configFile, ids, signalNamesSet)
except: 
  eprint(f"Read repository from {p.configFile} failed")
  raise

print(med.cerr())

eprint('Predicting')
model.apply(rep, samples)

# Write to file
eprint("Saving")
try: samples.write_to_file(outSamples)
except: 
  eprint(f'Cannot write predictinos to {p.outSamples}')
  raise

    
# Perfomance
if p.outPerformance != '':
  med.CommonLib.print_auc_performance(samples, range(samples.nSplits()), p.outPerformance)
    
