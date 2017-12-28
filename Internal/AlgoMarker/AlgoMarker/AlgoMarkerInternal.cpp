#include "AlgoMarkerInternal.h"
#include "AlgoMarkerErr.h"
#include <Logger/Logger/Logger.h>
#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL	LOG_DEF_LEVEL

//-------------------------------------------------------------------------------------------------------------------------
InputTester *InputTester::make_input_tester(int it_type)
{
	if (it_type == (int)INPUT_TESTER_TYPE_SIMPLE)
		return new InputTesterSimple;

	return NULL;
}

//-------------------------------------------------------------------------------------------------------------------------
int InputTester::name_to_input_tester_type(const string &name)
{
	if ((name == "simple") || (name == "SIMPLE") || (name == "Simple"))
		return (int)INPUT_TESTER_TYPE_SIMPLE;

	return (int)INPUT_TESTER_TYPE_UNDEFINED;
}

//-------------------------------------------------------------------------------------------------------------------------
void InputTesterSimple::input_from_string(const string &in_str)
{
	sf.init_from_string(in_str);
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
	if (rc == SanitySimpleFilter::Failed) return 0;
	if (rc == SanitySimpleFilter::Failed_Max_Nvals) return 0;
	if (rc == SanitySimpleFilter::Failed_Min_Nvals) return 0;
	if (rc == SanitySimpleFilter::Failed_Outliers) return 0;

	// if we are here the test could not be performed for some reason and we fail making the test, returning -1 in this case

	return -1;
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
			// FILTER	<filter type>|<filter params>|use_for_max_outliers_flag|external_rc|internal_rc|err_msg
			// # max_overall_outliers config
			// MAX_OVERALL_OUTLIERS	<number>

			if (fields[0] == "FILTER") {
				if (fields.size() >= 2) {
					int type = InputTester::name_to_input_tester_type(fields[1]);
					if (type == (int)INPUT_TESTER_TYPE_UNDEFINED)
						return -1;
					InputTester *i_test = InputTester::make_input_tester(type);
					if (i_test == NULL)
						return -1;
					if (fields.size() >= 3) i_test->tester_params = fields[2];
					if (fields.size() >= 4) i_test->max_outliers_flag = stoi(fields[3]);
					if (fields.size() >= 5) i_test->externl_rc = stoi(fields[4]);
					if (fields.size() >= 6) i_test->internal_rc = stoi(fields[5]);
					if (fields.size() >= 7) i_test->err_msg = fields[6];

					i_test->input_from_string(i_test->tester_params);

					testers.push_back(i_test);
				}
				else
					return -1;

			}
			else if (fields[0] == "TESTER_NAME") name = fields[1];
			else if (fields[0] == "MAX_OVERLALL_OUTLIERS") max_overall_outliers = stoi(fields[1]);
		}
	}

	return 0;
}


//-------------------------------------------------------------------------------------------------------------------------
// tests and stops at first cardinal failed test
int InputSanityTester::test_if_ok(MedRepository &rep, int pid, long long timestamp, int &nvals, int &noutliers, InputSanityTesterResult &res)
{
	res.external_rc = 0;
	res.internal_rc = 0;
	res.err_msg = "";

	int outliers_count = 0;
	for (auto &test : testers) {

		int t_nvals, t_noutliers;
		int rc = test->test_if_ok(rep, pid, timestamp, t_nvals, t_noutliers);
		//MLOG("###>>> pid %d time %d test %s : nvals %d nout %d rc %d\n", pid, timestamp, test->tester_params.c_str(), t_nvals, t_noutliers, rc);
		if (test->max_outliers_flag) outliers_count += t_noutliers;
		if (rc < 0) {
			res.external_rc = AM_ELIGIBILITY_ERROR;
			res.internal_rc = -2;
			res.err_msg = "Could not run filter on sample. ["+name+"]";

			return -1; // Failed running the test itself
		}

		if (rc == 0) {

			// we failed the test for a good reason and get out
			res.external_rc = test->externl_rc;
			res.internal_rc = test->internal_rc;
			res.err_msg = test->err_msg + "["+name+"]";

			return 0;
		}

	}

	if (outliers_count > max_overall_outliers) {
		res.external_rc = AM_ELIGIBILITY_ERROR;
		res.internal_rc = -3;
		res.err_msg = "Too many outliers detected (" + to_string(outliers_count) + ") ["+name+"]";
		return 0;
	}

	return 1;
}
