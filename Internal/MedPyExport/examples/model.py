#!/usr/bin/python
import sys
import argparse
import numpy as np
sys.path.insert(0,'/nas1/UsersData/avi/MR/Libs/Internal/MedPyExport/generate_binding/Release/rh-python27')
import med

#
# example program to do some simple modeling operations (json, rep, samples, learn, predict, get_matrix, etc)
#

#----------------------------------------------------------------------------------------------------------
def parse_args():
	parser = argparse.ArgumentParser(description='model.py arguments parser')
	parser.add_argument('--rep', type=str, help='repository to work with' )
	parser.add_argument('--samples', type=str, help='samples for chosen action (learn, predict, get_mat, get_json_mat)' )
	parser.add_argument('--csv', type=str, help='output csv for get_mat, get_json_mat' )
	parser.add_argument('--json', type=str, help='json for learn, get_json_mat' )
	parser.add_argument('--model', type=str, help='model for predict, get_mat' )
	parser.add_argument('--preds', type=str, help='output predictions for predict' )
	
	parser.add_argument('--learn', action='store_true', default=False, help='learn json on samples/rep, output model' )
	parser.add_argument('--predict', action='store_true', default=False, help='predict model on samples/rep, output preds' )
	parser.add_argument('--get_mat', action='store_true', default=False, help='get matrix for model on samples/rep, output csv' )
	parser.add_argument('--get_json_mat', action='store_true', default=False, help='get matrix for json on samples/rep output csv' )
	parser.add_argument('--get_mat_as_np', action='store_true', default=False, help='get matrix for model on samples/rep , gets it as np, and prints some stat' )
	parser.add_argument('--resave_model', action='store_true', default=False, help='rewriting a model in order to get rid of new versions ' )

	
	if len(sys.argv)==1:
		parser.print_help(sys.stderr)
		sys.exit(1)
	
	args = parser.parse_args()
	
	return args
	

#####################################################################################################################	
def main(argv):
	print "#####################"
	args = parse_args()
	print("arguments---->", args)
	
	model = med.Model()
	samples = med.Samples()
	
	if (args.learn or args.get_json_mat):
		print "Reading Json"
		model.init_from_json_file_with_alterations(args.json, '')
	else:
		print "Reading Model '%s'" % args.model
		model.read_from_file(args.model)
		
	if (args.resave_model):
		print "Resaving Model"
		model.write_to_file(args.model)
		sys.exit(1)
		
	signalNamesSet = model.get_required_signal_names()
	
	print "Reading Samples"
	samples.read_from_file(args.samples)
	ids = samples.get_ids()
	
	print "Reading Repository"
	rep = med.PidRepository()
	rep.read_all(args.rep, ids, signalNamesSet)
	

	if (args.learn):
		print "Learning model: "
		model.learn(rep, samples)
		model.write_to_file(args.model)
	
	if (args.predict):
		print "Predicting"
		model.apply(rep, samples)
		samples.write_to_file(args.preds)
	
	if (args.get_mat):
		print "Get Mat (+ preds)"
		model.apply(rep, samples)
		model.features.write_as_csv_mat(args.csv)
	
	if (args.get_json_mat):
		print "Get Json Mat"
		model.learn(rep, samples, med.ModelStage.LEARN_REP_PROCESSORS, med.ModelStage.APPLY_FTR_PROCESSORS)
		model.features.write_as_csv_mat(args.csv)
		
	if (args.get_mat_as_np):
		print "Get Mat As NP:"
		model.apply(rep, samples)
		names = model.features.get_feature_names()
		for name in names:
			vec = model.features.data[name]
			print name , vec.mean(), np.histogram(vec,10)
		fsamples = med.Samples()
		fsamples.import_from_sample_vec(model.features.samples)
		df_samples = fsamples.to_df()
		print df_samples

if __name__ == "__main__":
	main(sys.argv[1:])
