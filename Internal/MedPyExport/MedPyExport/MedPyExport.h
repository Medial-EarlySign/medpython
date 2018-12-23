#ifndef __MED_PY_EXPORT_H
#define __MED_PY_EXPORT_H


/* 
 *
 * Every file that should be exported by SWIG should appear in this include list:
 *
 *
 */

#include "MedPyExportExample.h"
#include "MPPidRepository.h"
#include "MPDictionary.h"
#include "MPSigExporter.h"
#include "MPSampleVecExporter.h"
#include "MPModel.h"
#include "MPSplit.h"
#include "MPTime.h"
#include "MPFeatures.h"
#include "MPSamples.h"
#include "MPMat.h"

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
 "Time"
#endif //SWIG

std::vector<std::string> get_public_objects();

#endif // !__MED_PY_EXPORT_H

