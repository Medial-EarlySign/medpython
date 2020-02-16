#pragma once

#include <MedProcessTools/MedProcessTools/FeatureGenerator.h>
#include <SerializableObject/SerializableObject/SerializableObject.h>

/**
* calculate drug coverage of prescription time in defined the time window. a value between 0 to 1.
*/
class DiabetesFinderGenerator : public FeatureGenerator {
public:

	// Constructor/Destructor
	DiabetesFinderGenerator() : FeatureGenerator() { 
		generator_type = FTR_GEN_DIABETES_FINDER; 
		//names.push_back("df"); 
		req_signals.push_back("Glucose");
	};
	~DiabetesFinderGenerator() {};

	/// The parsed fields from init command.
	int init(map<string, string>& mapper);

	// Naming
	void set_names() { if (names.empty()) names.push_back("FTR_" + int_to_string_digits(serial_id, 6) + ".DiabetesFinder"); tags.push_back("Diabetes"); }

	// Copy
	virtual void copy(FeatureGenerator *generator) { *this = *(dynamic_cast<DiabetesFinderGenerator *>(generator)); }

	// generate a new feature
	int _generate(PidDynamicRec& rec, MedFeatures& features, int index, int num, vector<float *> &_p_data);
	//float get_value(PidDynamicRec &rec, int idx, int time, int sig_outcomeTime);

	// Serialization
	ADD_CLASS_NAME(DiabetesFinderGenerator)
		ADD_SERIALIZATION_FUNCS(generator_type, names, tags, iGenerateWeights, req_signals)
};

MEDSERIALIZE_SUPPORT(DiabetesFinderGenerator);
