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
ap.add_argument("--matching",type=str, default='',
  help="optional matching parameters")
ap.add_argument("--samples",dest='samplesFile',type=str, required=True,
  help="training samples file")
ap.add_argument("--inModel",type=str, required=True,
  help="model json file")
ap.add_argument("--outModel",type=str, required=True,
  help="model output file")
ap.add_argument("--inpatient", default=False, action='store_true',
  help="indicate that relevant data is in-patient data")
ap.add_argument("--jsonAlt",type=str, dest='json_alt', nargs='+', default=[],
  help="(multiple) alterations to model json file")
p = ap.parse_args()

med.Global.RNG.srand(12131415)
med.Global.logger.init_level(med.Global.logger.LOG_MEDALGO_, 6)

"""
class Namespace(object): pass
p=Namespace()

p.inpatient = True
p.configFile = "/home/Repositories/THIN/thin_final/thin.repository"
p.samplesFile = "/nas1/UsersData/shlomi/MR/for_shlomi/mi_train.samples"
p.inModel = "/nas1/UsersData/shlomi/MR/for_shlomi/mi_model.json"
p.outModel = "/tmp/outfile.model"
p.json_alt = []
p.matching = ""
"""

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

eprint("Getting ids");
ids = samples.get_ids()

eprint("Reading Repository");
rep = med.PidRepository()
rep.read_all(p.configFile,ids, signalNamesSet)

#print(med.cerr())

#Match
if p.matching != '':
  matcher = med.SampleFilter.from_name_params("match", matchingParams);
  matcher.filter(rep, samples);


eprint("Learning");
model.learn(rep, samples)

#print(med.cerr())

eprint("Saving");

try: model.write_to_file(p.outModel)
except: eprint("Cannot write model to bin file %s\n", p.outModel);



    
