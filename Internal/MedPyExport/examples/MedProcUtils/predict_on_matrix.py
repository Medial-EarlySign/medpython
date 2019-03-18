#!/bin/env python3
import med
import sys
import argparse


def eprint(*args,**kwargs):
  kwargs['file']=sys.stderr
  print(*args,**kwargs)

    
ap = argparse.ArgumentParser()
ap.add_argument('--predictor',type=str, default='',
  help='predictor file')
ap.add_argument('--model',type=str, default='',
  help='model file')
ap.add_argument('--inCsv',type=str, default='',
  help='input matrix as csv file')
ap.add_argument('--inBin',type=str, default='',
  help='input matrix as bin file')
ap.add_argument('--outPred',type=str, required=True,
  help='predictions output file')
p = ap.parse_args()

    
    
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

med.Global.RNG.srand(12131415)
med.Global.logger.init_level(med.Global.logger.LOG_MEDALGO_, 6)

class Namespace(object): pass
p=Namespace()

p.params = '/home/Repositories/THIN/thin_jun2017/thin.repository'
p.predictor = '/nas1/UsersData/shlomi/MR/for_shlomi/mi_train.samples'
p.inModel = '/nas1/UsersData/shlomi/MR/for_shlomi/mi_model.json'
p.model = '/tmp/outfile.model'
p.json_alt = []

# Sanity
if bool(p.model) != bool(p.predictor):   # logical xor
    raise Exception('Either model or predictor must be given')

# Input options
samples = med.Samples()
if p.model:
    # Read Model
    eprint('Reading Model')
    model = med.Model()
    try: model.read_from_file(p.model)
    except:
        eprint(f'Cannot read model from binary file {p.model}')
        raise
    
    # Read Features
    eprint('Reading Features')
    readMatrix(model.features, p)
    
    # Apply
    eprint('Applying')
    dummyRep = med.PidRepository()
    model.features.get_samples(samples)
    model.apply(dummyRep, samples, med.ModelStage.APPLY_FTR_PROCESSORS, med.ModelStage.END)


    names = model.features.get_feature_names()
    for name in names: eprint(name)

else:

    # Read predictor
    eprint('Reading Predictor')
    pred = med.Predictor()
    med.CommonLib.read_predictor_from_file(pred, p.predictor)

    # Read Features
    eprint('Reading Features')
    features = med.Features()
    readMatrix(features, p)

    # Apply
    eprint('Predicting')
    pred.predict(features)
    features.get_samples(samples)

    # Write to samples file
    eprint('Writing')
    try: samples.write_to_file(p.outPred)
    except: 
        eprint(f'Cannot write predictions to file {p.outPred}')
        raise


