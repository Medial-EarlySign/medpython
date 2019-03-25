#!/bin/env python3
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
  help="training samples file")
ap.add_argument("--inModel",type=str, required=True,
  help="model json file")
ap.add_argument("--jsonAlt",type=str, dest='json_alt', nargs='+', default=[],
  help="(multiple) alterations to model json file")
ap.add_argument("--matching",type=str, default='',
  help="optional matching parameters")
ap.add_argument("--outPerformance",type=str, default='',
  help="performance output file")
ap.add_argument("--outSamples",type=str, default='',
  help="samples output file")
ap.add_argument("--inpatient", default=False, action='store_true',
  help="indicate that relevant data is in-patient data")
ap.add_argument("--nfolds", type=int, default=-1,
  help="if given, replace given splits with id-index%%nfolds")
ap.add_argument("--folds",type=str, default='',
  help="comma-separated list of folds to actually take for test. If not given - take all")
p = ap.parse_args()

med.Global.RNG.srand(12131415)
med.Global.logger.init_level(med.Global.logger.LOG_MEDALGO_, 6)


if p.inpatient:
  med.Global.default_time_unit = med.Time.Minutes
  med.Global.default_windows_time_unit = med.Time.Minutes
  eprint("Working on in-patient data");


eprint("Reading model");
model = med.Model()
model.init_from_json_file_with_alterations(p.inModel, p.json_alt)

eprint("Getting required signals")
signalNamesSet = model.get_required_signal_names()

eprint("Reading samples")
samples = med.Samples()
try: samples.read_from_file(p.samplesFile)
except:
  eprint("Cannot read samples from file %s", samplesFile.c_str());
  raise

# Override splits
# (this takes too much time in python)
"""
if p.nfolds!=-1:
  nfolds = int(p.nfolds)
  idx=0
  for idSample in samples.idSamples:
    split = idx%nfolds
    idSample.split = split
    for sample in idSample.samples:
      sample.split = split
    idx+=1
"""
if p.nfolds!=-1:
  samples.override_splits(int(p.nfolds))

# Folds to take
nSplits = samples.nSplits();
folds=[]

# get_folds()
if p.folds!="":
  for fold in [int(x) for x in p.folds.split(',')]:
    if fold >= nSplits or fold < 0:
      raise Exception(f'Illegal fold \'{fold}\' in folds')
    if fold in folds:
      raise Exception(f'Duplicate fold \'{fold}\' in folds \'{p.folds}\'')
    folds.append(fold)
else:
  for i in range(nSplits):
    folds.append(i)
    
eprint("Getting ids");
ids = samples.get_ids()

eprint("Reading Repository");
rep = med.PidRepository()
try: rep.read_all(p.configFile,ids, signalNamesSet)
except:
  eprint(f'Read repository from {p.configFile} failed')
  raise

#print(med.cerr())

#Match
matcher = None
if p.matching != '':
  matcher = med.SampleFilter.from_name_params("match", p.matching);

# Cross Validation
allTestSamples = med.Samples()
for iFold in folds:
  eprint(f'Working on split {iFold+1} / {nSplits}')

  #split
  trainSamples = med.Samples()
  testSamples = med.Samples()

  for idSamples in samples.idSamples:
    if idSamples.split == iFold:
      testSamples.idSamples.append(idSamples)
    else:
      trainSamples.idSamples.append(idSamples)

  # Match
  if matcher is not None:
    matcher.filter(rep, trainSamples)
    
  # Learn + apply
  model.clear()
  model.init_from_json_file_with_alterations(p.inModel, p.json_alt)
  model.learn(rep, trainSamples)
  model.apply(rep, testSamples)

  # Append to allTestSamples
  for idSamples in testSamples.idSamples:
    allTestSamples.idSamples.append(idSamples)

# Perfomance
if p.outPerformance!= '':
  med.CommonLib.print_auc_performance(allTestSamples, folds, p.outPerformance)

# Write Samples
if p.outSamples!='':
  try: allTestSamples.write_to_file(p.outSamples)
  except:
    perror(f'Cannot write predictinos to {p.outSamples}')
    raise


