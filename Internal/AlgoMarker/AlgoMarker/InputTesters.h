#pragma once
#include <string>
#include <MedProcessTools/MedProcessTools/SampleFilter.h>
#include <MedProcessTools/MedProcessTools/SerializableObject.h>

using namespace std;

typedef enum {
	INPUT_TESTER_TYPE_UNDEFINED = 0,
	INPUT_TESTER_TYPE_SIMPLE = 1
} InputTesterType;


//==============================================================================================================
// InputTester : holds a single tester - this is the base class
//==============================================================================================================
class InputTester : public SerializableObject {

public:
	// the type of the tester
	int type = (int)INPUT_TESTER_TYPE_UNDEFINED;

	// return code and messages to return in case of not passing the test
	int externl_rc;	 // rcs -1 and 0 are reserved 
	int internal_rc; // rcs -1 and 0 are reserved 
	string err_msg;

	string tester_params; // params for the internal tester

	// initialize from string 
	virtual void input_from_string(const string &in_str) { return; };

	// testing the tester on a given rep for a certain pid,timestamp
	// returns: 1: passes the test , 0: did not pass , -1: could not test
	// also returns: nvals (if relevant): number of tests in the window time defined in the test
	//               noutliers (if relevant) : number of outliers found
	virtual int test_if_ok(MedRepository &rep, int pid, long long timestamp, int &nvals, int &noutliers) { return -1; }

	// 1: good to go 0: did not pass -1: could not test
	int test_if_ok(MedRepository &rep, int pid, long long timestamp) {
		int nvals, noutliers;
		return test_if_ok(rep, pid, timestamp, nvals, noutliers);
	}

	// get a new InputTester
	static InputTester *make_input_tester(InputTesterType it_type);
	static InputTesterType name_to_input_tester_type(const string &name);
};
//==============================================================================================================

//==============================================================================================================
// InputTesterSimple : an implementation that is able to test one of the following tests:
// (1) test that the signal actually exist in name (in the signals list in the repository)
// (2) within a given window: minimal number of tests
// (3) within a given window: maximal number of outliers
// (4) count outliers within a given window
//
// Does this using the object SanitySimpleFilter defined in MeProcessTools/SampleFilter.h
//==============================================================================================================
class InputTesterSimple : public InputTester {

public:
	SanitySimpleFilter sf;

	InputTesterSimple() { type = (int)INPUT_TESTER_TYPE_SIMPLE; }

	void input_from_string(const string &in_str);
	int test_if_ok(MedRepository &rep, int pid, long long timestamp, int &nvals, int &noutliers); // 1: good to go 0: did not pass -1: could not test

};
//==============================================================================================================

//==============================================================================================================
// InputSanityTester : able to read a config file containing several tests and test them.
// Format of config file:
// # comment lines start with #
// NAME <name of tester : for debug prints, etc>
// # each filter defined using:
// FILTER	<filter type>|<filter params>|external_rc|internal_rc|err_msg
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
	int test_if_ok(MedRepository &rep, int pid, long long timestamp, int &nvals, int &noutliers); // tests and stops at first cardinal failed test

	 // tests and stops at first cardinal failed test 
	int test_if_ok(MedRepository &rep, int pid, long long timestamp) {
		int nvals, noutliers;
		return test_if_ok(rep, pid, timestamp, nvals, noutliers);
	}

	void clear() {
		for (auto &p_it : testers)
			if (p_it != NULL) delete p_it;
		testers.clear();
		max_overall_outliers = (int)1e9;
		name = "";
	}
};
//==============================================================================================================
