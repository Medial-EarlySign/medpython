#define AM_DLL_IMPORT

#include <AlgoMarker/AlgoMarker/AlgoMarker.h>
#include <AlgoMarker/DynAMWrapper/DynAMWrapper.h>

#define DYN(s) ((DynAM::initialized()) ? DynAM::s : s)

#include "CmdArgs.h"


void load_amlib(const string &amlib)
{
	// upload dynamic library if needed
	if (amlib != "") {
		load_am(amlib.c_str());
	}

	if (DynAM::initialized())
		MLOG("Dynamic %s library loaded\n", amlib.c_str());
	else
		MLOG("Dynamic library not loaded\n");

}

void initialize_algomarker(const string &amconfig, AlgoMarker *&test_am)
{
	// Load
	MLOG("Loading AM\n");
	int rc = DYN(AM_API_Load(test_am, amconfig.c_str()));
	if (rc != AM_OK_RC) MTHROW_AND_ERR("ERROR: Failed loading algomarker %s with config file %s ERR_CODE: %d\n", test_am->get_name(), amconfig.c_str(), rc);
	MLOG("Name is %s\n", test_am->get_name());
}

void init_and_load_data(const string &input_json_path, AlgoMarker *am) {
	DYN(AM_API_ClearData(am));

	string in_jsons;
	char * out_messages;
	if (read_file_into_string(input_json_path, in_jsons) < 0)
		MTHROW_AND_ERR("Error on loading file %s\n", in_jsons.c_str());
	MLOG("read %d characters from input jsons file %s\n", in_jsons.length(), input_json_path.c_str());
	int load_status = DYN(AM_API_AddDataByType(am, in_jsons.c_str(), &out_messages));
	if (out_messages != NULL) {
		string msgs = string(out_messages); //New line for each message:
		MLOG("AddDataByType has messages:\n");
		MLOG("%s\n", msgs.c_str());
	}
	DYN(AM_API_Dispose(out_messages));
	MLOG("Added data from %s\n", input_json_path.c_str());
	if (load_status != AM_OK_RC)
		MERR("Error code returned from calling AddDataByType: %d\n", load_status);
}

int get_response_score_into_sample(AMResponse *response, int resp_rc)
{
	int n_scores;
	int pid;
	long long ts;
	char *_scr_type = NULL;

	DYN(AM_API_GetResponseScoresNum(response, &n_scores));
	//int resp_rc = AM_API_GetResponse(resp, i, &pid, &ts, &n_scr, &_scr, &_scr_type);
	//MLOG("resp_rc = %d\n", resp_rc);
	//MLOG("i %d , pid %d ts %d scr %f %s\n", i, pid, ts, _scr, _scr_type);

	DYN(AM_API_GetResponsePoint(response, &pid, &ts));

	if (resp_rc == AM_OK_RC && n_scores > 0) {
		float _scr;
		resp_rc = DYN(AM_API_GetResponseScoreByIndex(response, 0, &_scr, &_scr_type));
		//MLOG("i %d , pid %d ts %d scr %f %s\n", i, pid, ts, _scr, _scr_type);
	}

	return n_scores;
}


int get_preds_from_algomarker_single(AlgoMarker *am,
	const string &sjreq, bool calc_by_type,
	int pred_id, int pred_time, bool print_resp = false)
{

	char *stypes[] = { "Raw" };
	int pid_id = pred_id;
	long long _timestamp = pred_time;
	char *jreq = (char *)(sjreq.c_str());
	char *jresp;
	//MLOG("Before Calculate jreq len %d\n", sjreq.length());
	MedTimer timer;	timer.start();

	if (calc_by_type) {
		DYN(AM_API_CalculateByType(am, JSON_REQ_JSON_RESP, jreq, &jresp));
		if (print_resp) {
			string s = string(jresp);
			MLOG("%s\n", s.c_str());
		}
		DYN(AM_API_Dispose(jresp));
	}
	else {
		AMResponses *resp;
		DYN(AM_API_CreateResponses(&resp));

		AMRequest *req;
		int req_create_rc = DYN(AM_API_CreateRequest("test_request", stypes, 1, &pid_id, &_timestamp, 1, &req));
		if (req_create_rc > 0)
			MTHROW_AND_ERR("Faield to create req\n");

		DYN(AM_API_Calculate(am, req, resp)); // calculate

		DYN(AM_API_GetResponsesNum(resp));

		AMResponse *response;
		int resp_rc = DYN(AM_API_GetResponseAtIndex(resp, 0, &response));
		get_response_score_into_sample(response, resp_rc);

		if (print_resp) {
			//string s = string(response);
			//MLOG("%s\n", s.c_str());
		}

		DYN(AM_API_DisposeRequest(req));
		DYN(AM_API_DisposeResponses(resp));
	}


	//((MedialInfraAlgoMarker *)am)->CalculateByType(JSON_REQ_JSON_RESP, jreq, &jresp);	
	//MLOG("After Calculate jresp len %d\n", strlen(jresp));


	timer.take_curr_time();

	return 0;
}

int main(int argc, char *argv[]) {
	ProgramArgs args;
	if (args.parse_parameters(argc, argv) < 0)
		return -1;

	load_amlib(args.amlib);
	AlgoMarker *test_am;

	MLOG("Creating AM\n");

	if (DYN(AM_API_Create((int)AM_TYPE_MEDIAL_INFRA, &test_am)) != AM_OK_RC)
		MTHROW_AND_ERR("ERROR: Failed creating test algomarker\n");

	if (args.test_stage > 0) //1 and up
		initialize_algomarker(args.amconfig, test_am);

	string sjreq = "";
	if (!args.jreq.empty() && read_file_into_string(args.jreq, sjreq) < 0)
		MTHROW_AND_ERR("Error can't read jreq %s\n", args.jreq.c_str());

	if (sjreq.empty()) {
		sjreq = "{ \"type\" : \"request\", \"request_id\" : \"my test\", \"export\" : {\"prediction\" : \"pred_0\"}, \"requests\" : [{ \"patient_id\": \"" + to_string(args.pid_id) +
			"\", \"time\" : \"" + to_string(args.prediction_time) + "\" }] }";
	}

	if (args.test_stage > 1) //2 and up
		init_and_load_data(args.in_jsons, test_am);

	//Apply 3+:
	if (args.test_stage > 2) //3 and up
		get_preds_from_algomarker_single(test_am, sjreq, args.calc_by_type, args.pid_id, args.prediction_time, true);

	if (args.test_stage > 1)
		DYN(AM_API_ClearData(test_am));

	DYN(AM_API_DisposeAlgoMarker(test_am));

	return 0;
}