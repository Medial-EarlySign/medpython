#include "AlgoMarker.h"

#include <Logger/Logger/Logger.h>
#include <MedTime/MedTime/MedTime.h>
#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL	LOG_DEF_LEVEL



//===========================================================================================================
// SimpleExampleEGFRAlgoMarker Implementations ::
// Follows is an implementation of the simplest AlgoMarker.
// This AlgoMarker calculates a single (not supporting batch) eGFR score
//===========================================================================================================
//------------------------------------------------------------------------------------------
int SimpleExampleEGFRAlgoMarker::AddData(int patient_id, const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, float* Values)
{
	if (Values_len != 1 || signalName == NULL)
		return AM_FAIL_RC;

	if (pid >= 0 && patient_id != pid)
		return AM_FAIL_RC;

	pid = patient_id;
	if (string(signalName) == "Creatinine")
		creatinine = Values[0];
	if (string(signalName) == "Age" || string(signalName) == "AGE")
		age = Values[0];
	if (string(signalName) == "Gender" || string(signalName) == "GENDER")
		gender = Values[0];
	if (string(signalName) == "Ethnicity" || string(signalName) == "ethnicity")
		ethnicity = (int)Values[0];

	return AM_OK_RC;

}

//------------------------------------------------------------------------------------------
int SimpleExampleEGFRAlgoMarker::Calculate(AMRequest *request, AMResponses *responses)
{
	// sanities
	if (responses == NULL)
		return AM_FAIL_RC;

	AMMessages *shared_msgs = responses->get_shared_messages();

	if (request == NULL) {
		string msg = "Error :: (" + to_string(AM_MSG_NULL_REQUEST) + " ) NULL request in Calculate()";
		shared_msgs->insert_message(AM_GENERAL_FATAL, msg.c_str());
		return AM_FAIL_RC;
	}


	if (request->get_n_score_types() !=1 || string(request->get_score_type(0)) != "Raw") {
		string msg = "Error (ExampleEGFR) :: Supporting only single score per request, and only Raw score type";
		shared_msgs->insert_message(AM_GENERAL_FATAL, msg.c_str());
		return AM_FAIL_RC;
	}

	// eligibility
	if (pid < 0 || age < 0 || (gender!=1 && gender !=2) || creatinine < 0.2 || creatinine > 100) {
		string msg = "Error (ExampleEGFR) :: Non eligible set: pid " + to_string(pid) + " age " + to_string(age) + " gender " + to_string(gender) + " creatinine " + to_string(creatinine);
		shared_msgs->insert_message(AM_GENERAL_FATAL, msg.c_str());
		return AM_FAIL_RC;
	}

	// scoring
	double eGFR_CKD_EPI = pow(0.993, (double)age);

	if (ethnicity == 1)
		eGFR_CKD_EPI *= 1.159;

	if (gender == 1) {
		// Male
		eGFR_CKD_EPI *= 141.0;
		if (creatinine <= 0.9)
			eGFR_CKD_EPI *= pow(creatinine/0.9, -0.441);
		else
			eGFR_CKD_EPI *= pow(creatinine/0.9, -1.209);
	}
	else {
		// Female
		eGFR_CKD_EPI *= 144.0;
		if (creatinine <= 0.7)
			eGFR_CKD_EPI *= pow(creatinine/0.7, -0.329);
		else
			eGFR_CKD_EPI *= pow(creatinine/0.7, -1.209);
	}


	// insert score to response
	int _pid = request->get_pid(0);
	long long _ts = request->get_timestamp(0);

	// create a response
	responses->set_request_id(request->get_request_id());
	char *stype = request->get_score_type(0);
	MLOG("stype = %s\n", stype);
	responses->insert_score_types(&stype, 1);

	AMResponse *res = responses->create_point_response(_pid, _ts);
	res->init_scores(1);
	int _n_score_types;
	char **stypes;
	responses->get_score_types(&_n_score_types, &stypes);
	MLOG("n %d stypes[0] = %s\n", _n_score_types, stypes[0]);

	res->set_score(0, (float)eGFR_CKD_EPI, stypes[0],"");

	return AM_OK_RC;
}

