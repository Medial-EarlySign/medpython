#!/bin/bash
set -e
set -x
set -o pipefail # makes sure a bad exit code is returned when appropriate

export OMP_NUM_THREADS=7
export PATH
echo "hostname:" `hostname`
echo "whoami:" `whoami`
echo "SHELL: $SHELL"
echo "PATH: $PATH"
echo "which python: " `which python`
echo "python --version: " 
python --version
exec $1 "${@:2}" |&  awk '{ print strftime("[%Y-%m-%d %H:%M:%S]"), $0 }'
