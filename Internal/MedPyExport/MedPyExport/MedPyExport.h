#ifndef __MED_PY_EXPORT_H
#define __MED_PY_EXPORT_H


/* 
 *
 * Every file that should be exported by SWIG should appear in this include list:
 *
 *
 */

#include "MPPidRepository.h"
#include "MPDictionary.h"
#include "MPSigExporter.h"
#include "MPModel.h"
#include "MPSplit.h"
#include "MPTime.h"
#include "MPFeatures.h"
#include "MPSamples.h"
#include "MPMat.h"
#include "MPGlobal.h"
#include "MPLogger.h"
#include "MPSampleFilter.h"
#include "MPSerializableObject.h"
#include "MPPredictor.h"
#include "MPBootstrap.h"

#ifndef SWIG
#define PUBLIC_OBJECTS "Model", \
 "Sample", \
 "PidRepository", \
 "Dictionary", \
 "FeatureAttr", \
 "Features", \
 "IdSamples", \
 "Mat", \
 "ModelStage", \
 "Samples", \
 "Sig", \
 "Split", \
 "Time", \
 "cerr", \
 "logger_use_stdout", \
 "Global", \
 "Logger", \
 "SampleFilter", \
 "CommonLib", \
 "PredictorTypes", \
 "Predictor", \
 "SampleVectorAdaptor", \
 "MedConvert", \
"Bootstrap"

#endif //SWIG

std::vector<std::string> get_public_objects();

#endif // !__MED_PY_EXPORT_H

