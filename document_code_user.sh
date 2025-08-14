#!/bin/bash
#sudo doxygen Doxyfile 2>&1 | egrep -v "External|QGDict|Please use a more specific name by including|Eigen::|VowpalWabbit|dmlc::|xgboost::|rabit::" > documentation_creation.log

cd ${0%/*}
OUTPUT_DIR=/home/$USER/html/Libs
mkdir -p ${OUTPUT_DIR}
( cat Doxyfile ; echo "OUTPUT_DIRECTORY=${OUTPUT_DIR}" ) | doxygen - 2>&1 | egrep -v "External|QGDict|Please use a more specific name by including|Eigen::|VowpalWabbit|dmlc::|xgboost::|rabit::" > documentation_creation_user.log
