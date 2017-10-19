#pragma once

//===============================================================================
// AlgoMarker.h
//-------------------------------------------------------------------------------
//
// Wrapper for a DLL that will contain :
// (1) All the needed API's to perform as an AlgoMarker
// (2) The option to run full MedProcessTools models with MedRepository inputs.
// (3) All in one single DLL managing it all.
//
//===============================================================================


// AM_DLL_EXPORT is defined only in the matching .cpp file to handle the dll building
// apps just include this h file and hence will work in import mode.

#if defined AM_DLL_IMPORT
#define DLL_WORK_MODE __declspec(dllimport)
#else
#define DLL_WORK_MODE __declspec(dllexport)
#endif


//
// includes of Medial Internal Libraries
//
#include "AlgoMarkerInternal.h"

//===============================================================================
// MedAlgoData - a framework to 
//===============================================================================
class DLL_WORK_MODE MedAlgoMarker {

private:
	MedAlgoMarkerInternal ma;

public:


};