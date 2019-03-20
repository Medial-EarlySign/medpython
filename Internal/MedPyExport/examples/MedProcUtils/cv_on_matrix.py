#!/bin/env python3
import med
import sys
import argparse


def eprint(*args,**kwargs):
  kwargs['file']=sys.stderr
  print(*args,**kwargs)

# Read A matrix from csv/bin
def readMatrix(features, vm):
  if vm.inCsv!='' and vm.inBin != '':
    raise Exception('Exactly one of inCsv and inBin should be given')

  if vm.inCsv != '':
    try: features.read_from_csv_mat(vm.inCsv)
    except:
      eprint(f'Cannot read features from csv file {vm.inCsv}')
  else:
    try: features.asSerializable().read_from_file(vm.inBin)
    except: 
        eprint(f'Cannot read features from bin file {vm.inBin}')
        raise

ap = argparse.ArgumentParser()
ap.add_argument('--predictor',type=str, default='',
  help='predictor type')
ap.add_argument('--params',type=str, default='',
  help='predictor parameters')
ap.add_argument('--inModel',type=str, default='',
  help='model json file')
ap.add_argument('--json_alt',type=str, dest='json_alt', nargs='+', default=[],
  help='(multiple) alterations to model json file')
ap.add_argument('--inCsv',type=str, default='',
  help='input matrix as csv file')
ap.add_argument('--inBin',type=str, default='',
  help='input matrix as bin file')
ap.add_argument('--outPerformance',type=str, default='',
  help='performance output file')
ap.add_argument("--outSamples",type=str, default='',
  help="samples output file")
ap.add_argument("--nfolds", type=int, default=-1,
  help="if given, replace given splits with id-index%%nfolds")
ap.add_argument("--folds",type=str, default='',
  help="comma-separated list of folds to actually take for test. If not given - take all")
p = ap.parse_args()


med.Global.RNG.srand(12131415)
med.Global.logger.init_level(med.Global.logger.LOG_MEDALGO_, 6)

# Sanity
if not(p.model or (p.predictor and p.params)):
    raise Exception('Either model or predictor+params must be given')

# Read Features
eprint('Reading Features')
features = med.Features()
readMatrix(features, p)

# Input options
model = med.Model()
if p.inModel:
    # Initialize predictor
    eprint('Initializing')

    # Read Model json
    model.init_from_json_file_with_alterations(p.inModel, p.jason_alt)
else:
    if not(p.predictor and p.params):
        raise Exception('Both predictor and params must be given')

    # Initialize predictor
    eprint('Initializing Predictor')
    model.make_predictor(p.predictor, p.params)
    

# Override splits
if p.nfolds:
  features.samples.override_splits(int(p.nfolds))
    
# Folds to take
nSplits = features.samples.nSplits()

folds=[]
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
        
# Cross Validation

dummyRep = med.PidRepository()
dummySamples = med.Samples()
allTestSamples = med.SampleVectorAdaptor()
for iFold in folds:
    eprint(f'Working on split {iFold+1}/{nSplits}')

    # Learn
    features.split_by_fold(model.features, iFold, True)
    eprint(f'{len(model.features.data)} X {len(model.features.samples)}')
    model.learn(dummyRep, dummySamples, med.ModelStage.LEARN_FTR_PROCESSORS, med.ModelStage.END)

    # Apply
    features.split_by_fold(model.features, iFold, False)
    eprint(f'{len(model.features.data)} X {len(model.features.samples)}')
    model.apply(dummyRep, dummySamples, med.ModelStage.APPLY_FTR_PROCESSORS, med.ModelStage.APPLY_PREDICTOR)

    # Append to allTestSamples
    allTestSamples.append_vec(model.features.samples)

# Perfomance
samples = med.Samples()
samples.import_from_sample_vec(allTestSamples)
if p.outPerformance:
    print_auc_performance(samples, folds, p.outPerformance)

# Write Samples
if p.outSamples:
    try: samples.write_to_file(p.outSamples)
    except: 
        eprint(f'Cannot write predictinos to {p.outSamples}')
        raise


