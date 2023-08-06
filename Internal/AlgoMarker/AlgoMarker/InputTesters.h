#pragma once
#include <string>
#include <MedProcessTools/MedProcessTools/SampleFilter.h>
#include <MedProcessTools/MedProcessTools/MedModel.h>
#include <SerializableObject/SerializableObject/SerializableObject.h>

using namespace std;

typedef enum {
	INPUT_TESTER_TYPE_UNDEFINED = 0,
	INPUT_TESTER_TYPE_SIMPLE = 1,
	INPUT_TESTER_TYPE_ATTR = 2,
	INPUT_TESTER_TYPE_FEATURE_JSON = 3
} InputTesterType;

typedef enum {
	TESTER_STAGE_UNDEFINED = 0,
	TESTER_STAGE_BEFORE_MODEL = 1,
	TESTER_STAGE_AFTER_MODEL = 2
} TesterStage;

//==============================================================================================================
// InputTester : holds a single tester - this is the base class
//==============================================================================================================
class InputTester : public SerializableObject {

public:
	// the type of the tester
	int type = (int)INPUT_TESTER_TYPE_UNDEFINED;

	int stage = (int)TESTER_STAGE_UNDEFINED;

	// the tester can be defined as a warning only
	int is_warning = 0;

	// return code and messages to return in case of not passing the test
	int externl_rc = 0;	 // rcs -1 and 0 are reserved 
	int internal_rc = 0; // rcs -1 and 0 are reserved 
	string err_msg = "";
	string cant_evel_msg = ""; ///< message when can't evalute test. result < 0

	int max_outliers_flag = 0; // use or not use the tester to accumulate outliers counts

	bool stop_processing_more_errors = false; ///< if true will stop process more errors

	string tester_params; // params for the internal tester

	// initialize from string 
	virtual void input_from_string(const string &in_str) { return; };

	// testing the tester on a given rep for a certain pid,timestamp
	// returns: 1: passes the test , 0: did not pass , -1: could not test
	// also returns: nvals (if relevant): number of tests in the window time defined in the test
	//               noutliers (if relevant) : number of outliers found
	virtual int test_if_ok(MedPidRepository &rep, int pid, long long timestamp, int &nvals, int &noutliers) { nvals = 0; noutliers = 0; return -1; }

	virtual int test_if_ok(MedSample &sample) { return -1; };				// 1: good to go 0: did not pass -1: could not test

	virtual int test_if_ok(int pid, long long timestamp,
		const unordered_map<string, unordered_set<string>> &dict_unknown) { return -1; }

	// 1: good to go 0: did not pass -1: could not test
	int test_if_ok(MedPidRepository &rep, int pid, long long timestamp) {
		int nvals, noutliers;
		return test_if_ok(rep, pid, timestamp, nvals, noutliers);
	}

	void print();

	// get a new InputTester
	static InputTester *make_input_tester(int it_type);
	static int name_to_input_tester_type(const string &name);

	virtual ~InputTester() {};
};
//==============================================================================================================

//==============================================================================================================
// InputTesterSimple : an implementation that is able to test one of the following tests:
// (1) test that the signal actually exist in name (in the signals list in the repository)
// (2) within a given window: minimal number of tests
// (3) within a given window: maximal number of outliers
// (4) count outliers within a given window
// (5) within a given window: maximal number of tests
//
// Does this using the object SanitySimpleFilter defined in MeProcessTools/SampleFilter.h
//==============================================================================================================
class InputTesterSimple : public InputTester {

public:
	SanitySimpleFilter sf;
	string err_message_template;

	InputTesterSimple() {
		type = (int)INPUT_TESTER_TYPE_SIMPLE;
		stage = (int)TESTER_STAGE_BEFORE_MODEL;
	}

	void input_from_string(const string &in_str);
	int test_if_ok(MedPidRepository &rep, int pid, long long timestamp, int &nvals, int &noutliers); // 1: good to go 0: did not pass -1: could not test

	int test_if_ok(int pid, long long timestamp, const unordered_map<string, unordered_set<string>> &dict_unknown);
};
//==============================================================================================================
// InputTesterAttr : an implementation that is able to test the attributes created in the samples file of a model
// (1) test that the attribute exists (it should be there or an error will be reported)
// (2) test its value is below some bound (<=)
//
// Does this by directly testing the given sample
//==============================================================================================================
class InputTesterAttr : public InputTester {

public:
	string attr_name;
	float attr_max_val;

	InputTesterAttr() {
		type = (int)INPUT_TESTER_TYPE_ATTR;
		stage = (int)TESTER_STAGE_AFTER_MODEL;
	}

	void input_from_string(const string &in_str);
	int init(map<string, string>& mapper);
	int test_if_ok(MedSample &sample);						// 1: good to go 0: did not pass -1: could not test

	int test_if_ok(int pid, long long timestamp, const unordered_map<string, unordered_set<string>> &dict_unknown) { return 1; }
};

///==============================================================================================================
/// InputTesterJsonFeature : an implementation that is able to test feature created by MedModel json
/// (1) Calculate a feature
/// (2) test its value is below some bound (<=)
///
/// Does this by directly testing the given sample
///==============================================================================================================
class InputTesterJsonFeature : public InputTester {
private:
	MedModel feature_generator;
	bool _learned = false;
	string resolved_feat_name = "";
	vector<string> req_signals;
public:
	bool is_binary_model = false; ///< if true it is trained model
	string json_model_path = ""; ///< realative path to am config, in same folder
	string feature_name = ""; ///< feature name to look for
	float feat_min_val = MED_MAT_MISSING_VALUE; ///< when missing value, no limit
	float feat_max_val = MED_MAT_MISSING_VALUE; ///< when missing value, no limit
	bool verbose_learn = true; ///< can control output to screen in first time
	bool verbose_apply = false; ///< can control output to screen on apply
	unordered_set<string> allow_missing_signals; ///< list of allowed signal to miss

	InputTesterJsonFeature() { 
		type = (int)INPUT_TESTER_TYPE_FEATURE_JSON; 
		stage = (int)TESTER_STAGE_BEFORE_MODEL;
	}

	void input_from_string(const string &in_str);
	int init(map<string, string>& mapper);
	/// 1: good to go 0: did not pass -1: could not test
	int test_if_ok(MedPidRepository &rep, int pid, long long timestamp, int &nvals, int &noutliers);

	int test_if_ok(int pid, long long timestamp, const unordered_map<string, unordered_set<string>> &dict_unknown) { return 1; }
};
//==============================================================================================================

struct InputSanityTesterResult {
	int external_rc = 0;
	int internal_rc = 0;
	string err_msg = "";
};
//==============================================================================================================
// InputSanityTester : able to read a config file containing several tests and test them.
// Format of config file:
// # comment lines start with #
// NAME <name of tester : for debug prints, etc>
// # each filter defined using:
// FILTER	<filter type>|<filter params>|warning_or_error|use_for_max_outliers_flag|external_rc|internal_rc|err_msg
// warining_or_error: values are WARNING or ERROR
// use_for_max_outliers_flag: ACC=yes or ACC=no
// # max_overall_outliers config
// MAX_OVERALL_OUTLIERS	<number>
//==============================================================================================================
class InputSanityTester {

public:
	vector<InputTester *> testers;
	int max_overall_outliers = (int)1e9;
	string name = "";


	~InputSanityTester() { clear(); }

	int read_config(const string &f_conf);

	int test_if_ok(int pid, long long timestamp, const unordered_map<string, unordered_set<string>> &dict_unknown, vector<InputSanityTesterResult> &res);

	// tests all simple testers (Before running model)
	int test_if_ok(MedPidRepository &rep, int pid, long long timestamp, int &nvals, int &noutliers, vector<InputSanityTesterResult> &res); // tests and stops at first cardinal failed test

	 // tests and stops at first cardinal failed test 
	int test_if_ok(MedPidRepository &rep, int pid, long long timestamp, vector<InputSanityTesterResult> &res) {
		int nvals, noutliers;
		return test_if_ok(rep, pid, timestamp, nvals, noutliers, res);
	}

	// tests all attr testers for a single given sample (After running model)
	int test_if_ok(MedSample &sample, vector<InputSanityTesterResult> &res);

	void clear() {
		for (auto &p_it : testers)
			if (p_it != NULL) delete p_it;
		testers.clear();
		max_overall_outliers = (int)1e9;
		name = "";
	}
};
//==============================================================================================================
