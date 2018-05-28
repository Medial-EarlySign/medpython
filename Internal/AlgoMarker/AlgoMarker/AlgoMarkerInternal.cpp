#include "AlgoMarkerInternal.h"
#include "AlgoMarkerErr.h"
#include <Logger/Logger/Logger.h>
#include <boost/algorithm/string.hpp>
#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL	LOG_DEF_LEVEL

//-------------------------------------------------------------------------------------------------------------------------
InputTester *InputTester::make_input_tester(int it_type)
{
	if (it_type == (int)INPUT_TESTER_TYPE_SIMPLE)
		return new InputTesterSimple;
	if (it_type == (int)INPUT_TESTER_TYPE_ATTR)
		return new InputTesterAttr;

	return NULL;
}

//-------------------------------------------------------------------------------------------------------------------------
int InputTester::name_to_input_tester_type(const string &name)
{
	if ((name == "simple") || (name == "SIMPLE") || (name == "Simple"))
		return (int)INPUT_TESTER_TYPE_SIMPLE;
	if ((name == "attr") || (name == "ATTR") || (name == "Attr"))
		return (int)INPUT_TESTER_TYPE_ATTR;

	return (int)INPUT_TESTER_TYPE_UNDEFINED;
}

//-------------------------------------------------------------------------------------------------------------------------
void InputTesterSimple::input_from_string(const string &in_str)
{
	sf.init_from_string(in_str);
}

//-------------------------------------------------------------------------------------------------------------------------
void InputTesterAttr::input_from_string(const string &in_str)
{
	this->init_from_string(in_str);
}


//-------------------------------------------------------------------------------------------------------------------------
int InputTesterAttr::init(map<string, string>& mapper)
{
	for (auto entry : mapper) {
		string field = entry.first;
		if (field == "attr_name" || field == "name") { attr_name = entry.second; }
		else if (field == "max") attr_max_val = stof(entry.second);
	}

	return 0;
}


//-------------------------------------------------------------------------------------------------------------------------
// 1: good to go 0: did not pass -1: could not test
int InputTesterSimple::test_if_ok(MedRepository &rep, int pid, long long timestamp, int &nvals, int &noutliers)
{
	MedSample s;

	s.id = pid;
	s.time = (int)timestamp;

	int rc = sf.test_filter(s, rep, nvals, noutliers);

	if (rc == SanitySimpleFilter::Passed) return 1;
	if (rc > 0) return 0;

	// if we are here the test could not be performed for some reason and we fail making the test, returning -1 in this case

	return -1;
}

//-------------------------------------------------------------------------------------------------------------------------
// (1) test that the attribute exists (it should be there or an error will be reported)
// (2) test its value is below some bound (<=)
// Does this by directly testing the given sample
// returns -1: can't test (no such attr) 0: failed test 1: all ok.
//-------------------------------------------------------------------------------------------------------------------------
int InputTesterAttr::test_if_ok(MedSample &sample)
{
	if (sample.attributes.find(attr_name) == sample.attributes.end())
		return -1;

	if (sample.attributes[attr_name] <= attr_max_val)
		return 1;
	
	return 0;
}

//-------------------------------------------------------------------------------------------------------------------------
int InputSanityTester::read_config(const string &f_conf)
{
	ifstream inf(f_conf);

	if (!inf)
		return -1;

	//MLOG("initializing sanity tester from file %s\n", f_conf.c_str());
	string curr_line;
	while (getline(inf, curr_line)) {
		if ((curr_line.size() > 1) && (curr_line[0] != '#')) {

			if (curr_line[curr_line.size()-1] == '\r')
				curr_line.erase(curr_line.size()-1);

			vector<string> fields;
			split(fields, curr_line, boost::is_any_of("|\t"));

			// Format of config file:
			// # comment lines start with #
			// TESTER_NAME <name of tester : for debug prints, etc>
			// # each filter defined using:
			// FILTER	<filter type>|<filter params>|warning_or_error|use_for_max_outliers_flag|external_rc|internal_rc|err_msg
			// warining_or_error: values are WARNING or ERROR 
			// use_for_max_outliers_flag: ACC=yes or ACC=no
			// # max_overall_outliers config
			// MAX_OVERALL_OUTLIERS	<number>

			if (fields[0] == "FILTER") {
				if (fields.size() >= 2) {
					int type = InputTester::name_to_input_tester_type(fields[1]);
					if (type == (int)INPUT_TESTER_TYPE_UNDEFINED) {
						MERR("ERROR: (1) : InputSanityTester::read_config() parsing filter %s\n", curr_line.c_str());
						return -1;
					}
					InputTester *i_test = InputTester::make_input_tester(type);
					if (i_test == NULL) {
						MERR("ERROR: (2) : InputSanityTester::read_config() parsing filter %s\n", curr_line.c_str());
						return -1;
					}
					if (fields.size() >= 3) i_test->tester_params = fields[2];
					if (fields.size() >= 4) {
						if (boost::iequals(fields[3], "WARNING") || boost::iequals(fields[3], "WARN"))
							i_test->is_warning = 1;
						else if (boost::iequals(fields[3], "ERROR") || boost::iequals(fields[3], "ERR"))
							i_test->is_warning = 0;
						else
							i_test->is_warning = 0; // default case if having format problems
					}

					if (fields.size() >= 5) {
						if (boost::iequals(fields[4], "ACCUMULATE=1") || boost::iequals(fields[4], "ACC=1"))
							i_test->max_outliers_flag = 1;
						else if (boost::iequals(fields[4], "ACCUMULATE=0") || boost::iequals(fields[4], "ACC=0"))
							i_test->max_outliers_flag = 0;
						else
							i_test->max_outliers_flag = 0; // default is having format problems;
					}
					if (fields.size() >= 6) i_test->externl_rc = stoi(fields[5]);
					if (fields.size() >= 7) i_test->internal_rc = stoi(fields[6]);
					if (fields.size() >= 8) i_test->err_msg = fields[7];

					i_test->input_from_string(i_test->tester_params);

					testers.push_back(i_test);
				}
				else {
					MERR("ERROR: (3) : InputSanityTester::read_config() parsing filter %s\n", curr_line.c_str());
					return -1;
				}

			}
			else if (fields[0] == "TESTER_NAME") name = fields[1];
			else if (fields[0] == "MAX_OVERLALL_OUTLIERS") max_overall_outliers = stoi(fields[1]);
		}
	}

	return 0;
}


//-------------------------------------------------------------------------------------------------------------------------
// tests all simple testers
//-------------------------------------------------------------------------------------------------------------------------
int InputSanityTester::test_if_ok(MedRepository &rep, int pid, long long timestamp, int &nvals, int &noutliers, vector<InputSanityTesterResult> &Results)
{
	int outliers_count = 0;
	int n_warnings = 0;
	int n_errors = 0;
	for (auto &test : testers) {

		if (test->type == INPUT_TESTER_TYPE_ATTR) continue; // these are tested elsewhere

		InputSanityTesterResult res;
		res.external_rc = 0;
		res.internal_rc = 0;
		res.err_msg = "";

		int t_nvals, t_noutliers;
		int rc = test->test_if_ok(rep, pid, timestamp, t_nvals, t_noutliers);
//		MLOG("###>>> pid %d time %d test %s : nvals %d nout %d rc %d\n", pid, timestamp, test->tester_params.c_str(), t_nvals, t_noutliers, rc);
		if (test->max_outliers_flag) outliers_count += t_noutliers;
		if (rc < 0) {
			res.external_rc = AM_ELIGIBILITY_ERROR;
			res.internal_rc = -2;
			res.err_msg = "Could not run filter on sample. ["+name+"]";

			Results.push_back(res);
			n_errors++;
		}

		if (rc == 0) {

			// we failed the test for a good reason and get out
			res.external_rc = test->externl_rc;
			res.internal_rc = test->internal_rc;
			res.err_msg = test->err_msg + "["+name+"]";

			Results.push_back(res);

			if (test->is_warning)
				n_warnings++;
			else
				n_errors++;
		}

		// rc == 1 : nothing to do - passed the test
	}

	if (outliers_count > max_overall_outliers) {
		InputSanityTesterResult res;
		res.external_rc = AM_ELIGIBILITY_ERROR;
		res.internal_rc = -3;
		res.err_msg = "Too many outliers detected (" + to_string(outliers_count) + ") ["+name+"]";
		Results.push_back(res);
		n_errors++;
	}


	// MLOG("###>>> pid %d n_errors %d n_warnings %d\n", pid, n_errors, n_warnings);
	if (n_errors > 0)
		return 0;

	return 1;
}


//-------------------------------------------------------------------------------------------------------------------------
// tests all attr testers on a given sample
//-------------------------------------------------------------------------------------------------------------------------
int InputSanityTester::test_if_ok(MedSample &sample, vector<InputSanityTesterResult> &Results)
{
	int n_warnings = 0;
	int n_errors = 0;

	for (auto &test : testers) {

		if (test->type != INPUT_TESTER_TYPE_ATTR) continue; // only these tested here

		InputSanityTesterResult res;
		res.external_rc = 0;
		res.internal_rc = 0;
		res.err_msg = "";

		int rc = test->test_if_ok(sample);
		if (rc < 0) {
			res.external_rc = AM_ELIGIBILITY_ERROR;
			res.internal_rc = -2;
			res.err_msg = "Could not find attribute " + ((InputTesterAttr *)test)->attr_name + ". Are you sure you're using a model that generates it?";

			Results.push_back(res);
			n_errors++;
		}

		if (rc == 0) {

			// we failed the test
			res.external_rc = test->externl_rc;
			res.internal_rc = test->internal_rc;
			res.err_msg = test->err_msg + "["+name+"]";

			Results.push_back(res);

			if (test->is_warning)
				n_warnings++;
			else
				n_errors++;
		}

		// rc == 1 : nothing to do - passed the test
	}

	if (n_errors > 0)
		return 0;

	return 1;
}