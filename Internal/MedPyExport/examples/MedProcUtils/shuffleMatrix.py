#!/bin/env python3
import med
import sys
import argparse

default_repository = '/home/Repositories/THIN/thin_jun2017/thin.repository'

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
ap.add_argument("--inCsv",type=str, default='',
  help='input csv matrix file')
ap.add_argument("--inBin",type=str, default='',
  help='input bin matrix file')
ap.add_argument("--outCsv",type=str, default='',
  help='output csv matrix file')
ap.add_argument("--outBin",type=str, default='',
  help='output bin matrix file')
p = ap.parse_args()

med.Global.RNG.srand(12131415)
med.Global.logger.init_level(med.Global.logger.LOG_MEDALGO_, 6)

# Initialize
if p.outCsv=='' and p.outBin == '':
  raise Exception('At least one of OutCsv and OutBin should be given')

# Read
matrix = med.Features()
readMatrix(matrix, p)

# Shuffle
eprint('Shuffling')
med.CommonLib.shuffleMatrix(matrix)

# Write to file
eprint('Saving')
if p.outCsv!='':
  try: matrix.write_as_csv_mat(p.outCsv)
  except: 
    eprint(f'Cannot write matrix to csv file {p.outCsv}')
    raise

if p.outBin!='':
  try: matrix.asSerializable().write_to_file(p.outBin)
  except:
    eprint(f'Cannot write matrix to bin file {p.outBin}')
    raise
