#!/bin/env python3
import med
import sys
import argparse


def eprint(*args,**kwargs):
  kwargs['file']=sys.stderr
  print(*args,**kwargs)

ap = argparse.ArgumentParser()
ap.add_argument("--predictor",type=str, default='',
  help='predictor type')
ap.add_argument("--params",type=str, default='',
  help='predictor parameters')
ap.add_argument("--inModel",type=str, required=True,
  help="model json file")
ap.add_argument("--json_alt",type=str, nargs='+', default=[],
  help="(multiple) alterations to model json file")
ap.add_argument("--inCsv",type=str, default='',
  help='input matrix as csv file')
ap.add_argument("--inBin",type=str, default='',
  help='input matrix as bin file')
ap.add_argument("--out", required=True, type=str,
  help='predictor/model output file')
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

if not(p.model or (p.predictor and p.params)):
    raise Exception('Either model or predictor+params must be given')

if p.inModel:
    # Initialize predictor
    eprint('Initializing')
    
    # Read Model json
    model = med.Model();
    model.init_from_json_file_with_alterations(p.inModel, p.json_alt)
    
    # Read Features
    eprint('Reading')
    readMatrix(model.features, p)
    
    # Learn
    eprint('Learning')
    dummyRep = med.PidRepository()
    samples = med.Samples()
    model.features.get_samples(samples)
    model.learn(dummyRep, samples, med.ModelStage.LEARN_FTR_PROCESSORS, med.ModelStage.END)

    # Write to file
    eprint("Saving")
    
    try: model.write_to_file(p.out)
    except: 
        eprint(f'Cannot write model to binary file {p.out}')
        raise

else:
    if not(p.predictor and p.params):
        raise Exception('Both predictor and params must be given')


    # Initialize predictor
    eprint('Initializing Predictor')

    pred = med.Predictor(p.predictor, p.params)

    # Read Features
    eprint('Reading Features')
    features = med.Features()
    readMatrix(features, p)

    # Learn
    eprint('Learning')
    pred.learn(features)

    # Write to file
    eprint('Saving')
    write_predictor_to_file(pred, p.out);
    
