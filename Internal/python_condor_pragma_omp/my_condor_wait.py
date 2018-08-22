#!/usr/bin/python

from __future__ import print_function
import argparse
import subprocess
from collections import defaultdict 
from  medial_tools.medial_tools import eprint, read_file, write_for_sure
import re
from subprocess import Popen, PIPE, STDOUT
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--log_file', default='condor.log')
parser.add_argument('--seconds_between_attempts', type = int, default=600)
args = parser.parse_args()

cmd = "condor_wait " + args.log_file + " -wait " + str(args.seconds_between_attempts)
done_waiting = False
while not done_waiting: 
	try:
		print(cmd)
		subprocess.check_output(cmd, shell=True)
		done_waiting = True
	except subprocess.CalledProcessError:
		print("still not done...")
		p = Popen('condor_q | tail', stdout = PIPE, stderr = STDOUT, shell = True)
		for line in p.stdout:
			if "0 jobs; 0 completed, 0 removed, 0 idle, 0 running, 0 held, 0 suspended" in line:
				print("detected that all condor jobs are completed even though not all jobs are listed as terminated in [{0}]! done waiting".format(args.log_file,))
				done_waiting=True
			print(line.replace('\n', ''))
		sys.stdout.flush()

regex = r"\((\d+)\.(\d+)\..+\).+\n.+\(return value (\d+)\)"
with open (args.log_file, "r") as myfile:
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
	
	
