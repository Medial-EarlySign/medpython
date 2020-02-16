#include "DiabetesFinderGenerator.h"

#define LOCAL_SECTION LOG_FTRGNRTR
#define LOCAL_LEVEL	LOG_DEF_LEVEL

//=======================================================================================
// DiabetesFinderGenerator : Implementation for the Diabetes Finder AM 
//=======================================================================================


// Generate
//.......................................................................................
int DiabetesFinderGenerator::_generate(PidDynamicRec& rec, MedFeatures& features, int index, int num, vector<float *> &_p_data) {

	float *p_feat = _p_data[0] + index;

	UniversalSigVec usv;
	rec.uget("Glucose", 0, usv);
	for (int i = 0; i < num; i++)
		p_feat[i] = 0;

	for (int i = 0; i < num; i++) {
		int s_time = features.samples[index + i].time;
		int last_glucose_idx = -1;
		for (int j = usv.len-1; j >= 0 ; --j) {
			if (usv.Time(j) <= s_time)
			{
				last_glucose_idx = j;
				break;
			}
		}
		if (last_glucose_idx == -1)
			continue;
		float glucose = usv.Val(last_glucose_idx);

		if (glucose < 100.0f)
		{
			p_feat[i] = 1.0f;
			features.samples[index + i].str_attributes["DiabetesStatus"] = "Non-diabetic";
			features.samples[index + i].str_attributes["Urgency"] = "None";
		}
		else if (glucose < 110) {
			p_feat[i] = 2.0f;
			features.samples[index + i].str_attributes["DiabetesStatus"] = "Pre-diabetic-level1";
			features.samples[index + i].str_attributes["Urgency"] = "Low";
		}
		else if (glucose < 120) {
			p_feat[i] = 3.0f;
			features.samples[index + i].str_attributes["DiabetesStatus"] = "Pre-diabetic-level2";
			features.samples[index + i].str_attributes["Urgency"] = "Intermediate";
		}
		else { 
			p_feat[i] = 4.0f; 
			features.samples[index + i].str_attributes["DiabetesStatus"] = "Diabetic";
			features.samples[index + i].str_attributes["Urgency"] = "High";
		}

	}

	return 0;
}

// Init
//.......................................................................................
int DiabetesFinderGenerator::init(map<string, string>& mapper) {

	for (auto entry : mapper) {
		string field = entry.first;
		if (field == "tags")
			boost::split(tags, entry.second, boost::is_any_of(","));
		else if (field != "fg_type")
			MTHROW_AND_ERR("Unknown parameter [%s] for DiabetesFinderGenerator\n", field.c_str());
	}
	set_names();

	req_signals.clear();
	req_signals.push_back("Glucose");

	return 0;
}


