#!/usr/bin/python

from __future__ import print_function
import argparse
import subprocess, os
from collections import defaultdict 
from  medial_tools.medial_tools import eprint, read_file, write_for_sure
import re, time
from subprocess import Popen, PIPE, STDOUT
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--submit_file', default='condor_runner')
parser.add_argument('--seconds_between_attempts', type = int, default=600)
parser.add_argument('--log_file', required=True)
args = parser.parse_args()

log_file = args.log_file
cmd = "condor_submit " + args.submit_file + ' -append "log = {0}" '.format(log_file,)
print(cmd)
p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
output = p.stdout.read()
print(output)
sys.stdout.flush()
p.communicate()
if p.returncode != 0:
	print("condor_submit returncode == [{0}], exiting!".format(p.returncode,))
	exit(1)
if not os.path.exists(log_file):
	print("condor_submit did not generate the expected log_file [{0}], exiting!".format(log_file,))
	exit(1)
cmd = "condor_wait " + log_file + " -wait " + str(args.seconds_between_attempts)
done_waiting = False
while not done_waiting: 
	try:
		print(cmd)
		subprocess.check_output(cmd, shell=True)
		done_waiting = True
	except subprocess.CalledProcessError:
		
		print("still not done...")
		p = Popen('condor_q -wide| tail', stdout = PIPE, stderr = STDOUT, shell = True)
		for line in p.stdout:
			if "0 jobs; 0 completed, 0 removed, 0 idle, 0 running, 0 held, 0 suspended" in line:
				print("detected that all condor jobs are completed even though not all jobs are listed as terminated in [{0}]! done waiting".format(log_file,))
				done_waiting=True
			print(line.replace('\n', ''))
		sys.stdout.flush()

regex = r"\((\d+)\.(\d+)\..+\).+\n.+\(return value (\d+)\)"
with open (log_file, "r") as myfile:
	data = myfile.read()
	matches = re.finditer(regex, data)
	l = []
	for matchNum, match in enumerate(matches):
		matchNum = matchNum + 1
		e = {"cluster": int(match.group(1)), "job_id": int(match.group(2)), "rc": int(match.group(3))}
		l.append(e)		
	latest_cluster = l[-1]["cluster"]
	rcs = defaultdict(list)
	for e in l[::-1]:
		if e["cluster"] == latest_cluster:
			rcs[e["rc"]].append(e["job_id"])
		else:
			break
	print("latest cluster was:", latest_cluster)
	print("jobs return codes:")
	s = 0
	for k, v in rcs.items():
		print(k, ":", len(v), v)
		s += len(v)
	print("% of jobs with return code [0]:", 1.0 * len(rcs[0]) / s)
	if 1.0 * len(rcs[0]) / s < 0.8:
		raise Exception("too many failed jobs!")
	
	
