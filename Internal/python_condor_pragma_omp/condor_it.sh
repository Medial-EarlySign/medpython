#!/bin/bash
set -e
set -x
set -o pipefail # makes sure a bad exit code is returned when appropriate

export OMP_NUM_THREADS=7
export PATH
#export LD_LIBRARY_PATH=/server/Work/Libs/Boost/latest/stage/lib
echo "hostname:" `hostname`
echo "whoami:" `whoami`
echo "SHELL: $SHELL"
echo "PATH: $PATH"
echo "MR_ROOT: [$MR_ROOT]"
echo "LD_LIBRARY_PATH: [$LD_LIBRARY_PATH]"
echo "which python: " `which python`
echo "python --version: " 
python --version
exec $1 "${@:2}" |&  awk '{ print strftime("[%Y-%m-%d %H:%M:%S]"), $0 }'
