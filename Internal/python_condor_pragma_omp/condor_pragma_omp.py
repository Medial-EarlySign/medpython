#!/usr/bin/python

from __future__ import print_function
from medial_tools.medial_tools import eprint, read_file, write_for_sure, run_cmd_and_print_outputs
import os 
import stat

def before_condor_loop(do_condor_flag, condor_folder, run_id, condor_runner_template):
	if not os.path.exists(condor_folder + "/" + run_id + "/"):
		eprint('creating [{0}]'.format(condor_folder + "/" + run_id + "/"))
		os.makedirs(condor_folder + "/" + run_id + "/")

	condor_filename = condor_folder + "/" + run_id + "/" + run_id + ".condor_runner"
	eprint ("generating [{0}]".format(condor_filename))
	condor_file = open(condor_filename, "w")
	if do_condor_flag == 1:
		for l in open(condor_runner_template):
			print(l, file=condor_file)
		my_folder = os.path.dirname(os.path.realpath(__file__))
		print("executable = " + my_folder + "/condor_it.sh ", file=condor_file)
	else:
		print("#!/bin/bash", file=condor_file)
		print("set -e", file=condor_file)
		print("set -x", file=condor_file)
		print("set -o pipefail", file=condor_file) # makes sure a bad exit code is returned when appropriate
	return condor_filename, condor_file

def start_condor_loop(do_condor_flag, condor_folder, nice_name, run_id):
	cmd_filename = condor_folder + "/" + run_id + "/" + nice_name + "_" + run_id + ".sh"
	cmd_file = open(cmd_filename, 'w')
	print("#!/bin/bash", file=cmd_file)
	print("set -e", file=cmd_file)
	print("set -x", file=cmd_file)
	return cmd_filename, cmd_file

def end_condor_loop(do_condor_flag, condor_folder, nice_name, run_id, cmd_filename, cmd_file, condor_filename, condor_file):
	cmd_file.close()
	st = os.stat(cmd_filename)
	os.chmod(cmd_filename, st.st_mode | stat.S_IEXEC)
	my_folder = os.path.dirname(os.path.realpath(__file__))
	if do_condor_flag == 1:
		print("error = " + condor_folder + "/" + run_id + "/" + nice_name + "_" + run_id + "_$(Cluster).$(Process).err", file=condor_file)
		print("output = " + condor_folder + "/" + run_id + "/" + nice_name + "_" + run_id + "_$(Cluster).$(Process).out", file=condor_file)
		print("arguments = " + cmd_filename, file=condor_file)
		print("queue", file=condor_file)
		print("#####", file=condor_file)
	else:
		print (my_folder + "/condor_it.sh {0} | tee {1}".format(cmd_filename, condor_folder + "/" + run_id + "/" + nice_name + "_" + run_id + ".out"), file=condor_file)

def after_condor_loop(do_condor_flag, condor_folder, run_id, condor_filename, condor_file):
	condor_file.close()
	my_folder = os.path.dirname(os.path.realpath(__file__))
	if do_condor_flag == 1:
		run_cmd_and_print_outputs("python {my_folder}/my_condor_submit_and_wait.py --submit_file {0} --log_file {1}".format(condor_filename, condor_folder + "/" + run_id + "/" + run_id + ".condor_log",
			my_folder=my_folder))
	else:
		run_cmd_and_print_outputs(condor_filename)

