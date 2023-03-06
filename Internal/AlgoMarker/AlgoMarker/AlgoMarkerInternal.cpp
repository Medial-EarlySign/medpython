#include "AlgoMarkerInternal.h"
#include "AlgoMarkerErr.h"
#include <Logger/Logger/Logger.h>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL	LOG_DEF_LEVEL

//-------------------------------------------------------------------------------------------------------------------------
InputTester *InputTester::make_input_tester(int it_type)
{
	if (it_type == (int)INPUT_TESTER_TYPE_SIMPLE)
		return new InputTesterSimple;
	if (it_type == (int)INPUT_TESTER_TYPE_ATTR)
		return new InputTesterAttr;
	if (it_type == (int)INPUT_TESTER_TYPE_FEATURE_JSON)
		return new InputTesterJsonFeature;

	return NULL;
}

//-------------------------------------------------------------------------------------------------------------------------
int InputTester::name_to_input_tester_type(const string &name)
{
	if ((name == "simple") || (name == "SIMPLE") || (name == "Simple"))
		return (int)INPUT_TESTER_TYPE_SIMPLE;
	if ((name == "attr") || (name == "ATTR") || (name == "Attr"))
		return (int)INPUT_TESTER_TYPE_ATTR;
	if ((name == "feature") || (name == "feature_json"))
		return (int)INPUT_TESTER_TYPE_FEATURE_JSON;

	return (int)INPUT_TESTER_TYPE_UNDEFINED;
}

//-------------------------------------------------------------------------------------------------------------------------
void InputTester::print()
{
	MLOG("InputTester: type %d is_warning %d external_rc %d internal_rc %d max_outliers_flag %d err \'%s\' \n",
		type, is_warning, externl_rc, internal_rc, max_outliers_flag, err_msg.c_str());
}

//-------------------------------------------------------------------------------------------------------------------------
void InputTesterSimple::input_from_string(const string &in_str)
{
	sf.init_from_string(in_str);
	err_message_template = err_msg;
}

//-------------------------------------------------------------------------------------------------------------------------
void InputTesterAttr::input_from_string(const string &in_str)
{
	this->init_from_string(in_str);
}
//-------------------------------------------------------------------------------------------------------------------------
void InputTesterJsonFeature::input_from_string(const string &in_str)
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
int InputTesterSimple::test_if_ok(MedPidRepository &rep, int pid, long long timestamp, int &nvals, int &noutliers)
{
	err_msg = err_message_template;
	MedSample s;

	s.id = pid;
	s.time = (int)timestamp;

	int rc = sf.test_filter(s, rep, nvals, noutliers);

	if (rc == SanitySimpleFilter::Passed) return 1;
	if (rc > 0) return 0;

	// if we are here the test could not be performed for some reason and we fail making the test, returning -1 in this case

	return -1;
}

int InputTesterSimple::test_if_ok(int pid, long long timestamp, const unordered_map<string, unordered_set<string>> &dict_unknown) {
	if (sf.values_in_dictionary) {
		//Will be actually tested in here:
		if (dict_unknown.find(sf.sig_name) == dict_unknown.end())
			return 1;
		else {
			//Has bad codes:
			//change err_msg based on err_message_template: replace $SIGNAL with signal and $VALUE with missing value:
			string base_msg = boost::replace_all_copy(err_message_template, "$SIGNAL", sf.sig_name);
			
			bool start = true;
			if (base_msg.find("$VALUE") != string::npos) {
				stringstream ss;
				for (const string &val : dict_unknown.at(sf.sig_name)) {
					if (!start)
						ss << ",";
					ss << val;
					start = false;
				}
				err_msg = boost::replace_all_copy(base_msg, "$VALUE", ss.str());
			}
			else
				err_msg = base_msg;
			
			return 0;
		}
	}
	return 1;
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

int InputTesterJsonFeature::test_if_ok(MedPidRepository &rep, int pid, long long timestamp, int &nvals, int &noutliers)
{
	MedSample s;
	s.id = pid;
	s.time = (int)timestamp;
	int rc = -1;
	nvals = 0;
	noutliers = 0;
	//Assume MedModel was initialized once by init 
	//=> in first test we will learn the mode if not learned, and prepare for apply - to increase speed.
	//So first test is slower and rest are faster. "Cold start"
	//Not thread safe - except using OMP

	//Check all signals exists in repository:
	bool fail_signals = false;
	for (const string &req_s : req_signals) {
		int sid = rep.sigs.sid(req_s);
		if (sid < 0) {
			fail_signals = true;
			MLOG("Missing signal %s model eligibility testing\n", req_s.c_str());
		}
		//only when not allowed - not in allowed list:
		if (allow_missing_signals.find(req_s) == allow_missing_signals.end())

			if ((!rep.in_mem_mode_active() && !rep.index.index_table[sid].full_load) ||
				(rep.in_mem_mode_active() && rep.in_mem_rep.data.find(pair<int, int>(pid, sid)) == rep.in_mem_rep.data.end())) {
				fail_signals = true;
				MLOG("Signal %s not loaded for model eligibility testing\n", req_s.c_str());
			}

		if (fail_signals)
			break;
	}

	if (fail_signals)
		return -1;

	if (!_learned) {
#pragma omp critical(InputTesterJsonFeature)
		{
			MedSamples samples;
			//samples can be empty - can't use components that needs data and training
			if (!verbose_learn)
				medial::print::medmodel_logging(false);
			if (!is_binary_model)
				feature_generator.learn(rep, samples, MedModelStage::MED_MDL_LEARN_REP_PROCESSORS, MED_MDL_APPLY_FTR_PROCESSORS);
			try {
				feature_generator.init_model_for_apply(rep,
					MedModelStage::MED_MDL_LEARN_REP_PROCESSORS, MedModelStage::MED_MDL_APPLY_FTR_PROCESSORS);

				if (is_binary_model)  //need to generate matrix
					feature_generator.no_init_apply(rep, samples, MedModelStage::MED_MDL_LEARN_REP_PROCESSORS, MedModelStage::MED_MDL_APPLY_FTR_PROCESSORS);

				if (!verbose_learn)
					medial::print::medmodel_logging(true);

				vector<string> all_names;
				feature_generator.features.get_feature_names(all_names);
				resolved_feat_name = all_names[find_in_feature_names(all_names, feature_name, true)];

				_learned = true;
			}
			catch (...) {
				MERR("Error failed in Learn\n");
			}
		}
	}

	if (!_learned)
		return -2;

	//Model is prepared for apply and generation of feature and stores exact name in resolved_feat_name
	MedSamples smps;
	smps.import_from_sample_vec({ s });
	if (!verbose_apply)
		medial::print::medmodel_logging(false);
	try {
		feature_generator.no_init_apply(rep, smps, MedModelStage::MED_MDL_LEARN_REP_PROCESSORS, MedModelStage::MED_MDL_APPLY_FTR_PROCESSORS);
	}
	catch (...) {
		return -2;
	}
	if (!verbose_apply)
		medial::print::medmodel_logging(true);
	const vector<float> &vals = feature_generator.features.data.at(resolved_feat_name);
	float f_res = vals[0];

	rc = int(
		((feat_min_val == MED_MAT_MISSING_VALUE) || (f_res >= feat_min_val)) &&
		((feat_max_val == MED_MAT_MISSING_VALUE) || (f_res <= feat_max_val))
		);

	return rc;
}

int InputTesterJsonFeature::init(map<string, string>& mapper) {
	string base_path;
	for (auto &it : mapper)
	{
		if (it.first == "feature_name")
			feature_name = it.second;
		else if (it.first == "base_path")
			base_path = it.second;
		else if (it.first == "feat_min_val")
			feat_min_val = med_stof(it.second);
		else if (it.first == "feat_max_val")
			feat_max_val = med_stof(it.second);
		else if (it.first == "json_model_path")
			json_model_path = it.second;
		else if (it.first == "is_binary_model")
			is_binary_model = med_stoi(it.second) > 0;
		else if (it.first == "verbose_learn")
			verbose_learn = med_stoi(it.second) > 0;
		else if (it.first == "verbose_apply")
			verbose_apply = med_stoi(it.second) > 0;
		else if (it.first == "allow_missing_signals")
			boost::split(allow_missing_signals, it.second, boost::is_any_of(","));
		else
			MTHROW_AND_ERR("Error - unknown arg %s\n", it.first.c_str());
	}

	if (feature_name.empty())
		MTHROW_AND_ERR("Error - InputTesterJsonFeature::init - must specify feature_name\n");
	if (json_model_path.empty())
		MTHROW_AND_ERR("Error - InputTesterJsonFeature::init - must specify json_model_path\n");
	//init MedModel from json:
	if (!file_exists(json_model_path))
		json_model_path = base_path + path_sep() + json_model_path;
	if (!is_binary_model)
		feature_generator.init_from_json_file(json_model_path);
	else
		feature_generator.read_from_file(json_model_path);

	feature_generator.get_required_signal_names(req_signals);
	return 0;
}

//-------------------------------------------------------------------------------------------------------------------------
int InputSanityTester::read_config(const string &f_conf)
{
	ifstream inf(f_conf);

	if (!inf)
		return -1;
	string base_path = boost::filesystem::path(f_conf.c_str()).parent_path().string();

	//MLOG("initializing sanity tester from file %s\n", f_conf.c_str());
	string curr_line;
	while (getline(inf, curr_line)) {
		if ((curr_line.size() > 1) && (curr_line[0] != '#')) {

			if (curr_line[curr_line.size() - 1] == '\r')
				curr_line.erase(curr_line.size() - 1);

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
				//MLOG("Parsing FILTER line: %s\n", curr_line.c_str());
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
					if (fields.size() >= 9) i_test->cant_evel_msg = fields[8];
					if (fields.size() >= 10) i_test->stop_processing_more_errors = stoi(fields[9]) > 0;

					i_test->input_from_string("base_path=" + base_path + ";" + i_test->tester_params);

					//i_test->print();
					//MLOG("testers size is: %d\n", testers.size());

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
int InputSanityTester::test_if_ok(MedPidRepository &rep, int pid, long long timestamp, int &nvals, int &noutliers, vector<InputSanityTesterResult> &Results)
{
	int outliers_count = 0;
	int n_warnings = 0;
	int n_errors = 0;
	for (auto &test : testers) {

		if (test->stage != TESTER_STAGE_BEFORE_MODEL) continue; // these are tested elsewhere

		InputSanityTesterResult res;
		res.external_rc = 0;
		res.internal_rc = 0;
		res.err_msg = "";

		int t_nvals, t_noutliers;
		int rc = test->test_if_ok(rep, pid, timestamp, t_nvals, t_noutliers);
		//MLOG("###>>> pid %d time %d test %s : nvals %d nout %d rc %d\n", pid, timestamp, test->tester_params.c_str(), t_nvals, t_noutliers, rc);
		if (test->max_outliers_flag) outliers_count += t_noutliers;
		if (rc < 0) {
			res.external_rc = AM_ELIGIBILITY_ERROR;
			res.internal_rc = -2;
			//res.err_msg = "Could not run filter on sample. ["+name+"]";
			if (!test->cant_evel_msg.empty())
				res.err_msg = test->cant_evel_msg;
			else
				res.err_msg = "Could not run filter on sample."; // no name in err message (a request...)

			Results.push_back(res);
			n_errors++;
			if (test->stop_processing_more_errors)
				break;
		}

		if (rc == 0) {

			// we failed the test for a good reason and get out
			res.external_rc = test->externl_rc;
			res.internal_rc = test->internal_rc;
			//res.err_msg = test->err_msg + "["+name+"]";
			res.err_msg = test->err_msg; // no name in err message (a request...)

			//MLOG("###>>>>>> found an error: %d %d %s\n", res.external_rc, res.internal_rc, res.err_msg.c_str());

			Results.push_back(res);

			if (test->is_warning)
				n_warnings++;
			else
				n_errors++;

			if (test->stop_processing_more_errors)
				break;
		}

		// rc == 1 : nothing to do - passed the test
	}

	if (outliers_count > max_overall_outliers) {
		InputSanityTesterResult res;
		res.external_rc = AM_ELIGIBILITY_ERROR;
		res.internal_rc = -3;
		//res.err_msg = "Too many outliers detected (" + to_string(outliers_count) + ") ["+name+"]";
		res.err_msg = "Too many outliers detected (" + to_string(outliers_count) + ")"; // no name in err message (a request...)
		Results.push_back(res);
		n_errors++;
	}


	//MLOG("###>>> pid %d n_errors %d n_warnings %d\n", pid, n_errors, n_warnings);
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

		if (test->stage != TESTER_STAGE_AFTER_MODEL) continue; // only these tested here

		InputSanityTesterResult res;
		res.external_rc = 0;
		res.internal_rc = 0;
		res.err_msg = "";

		int rc = test->test_if_ok(sample);
		if (rc < 0) {
			res.external_rc = AM_ELIGIBILITY_ERROR;
			res.internal_rc = -2;
			if (!test->cant_evel_msg.empty())
				res.err_msg = test->cant_evel_msg;
			else
				res.err_msg = "Could not find attribute " + ((InputTesterAttr *)test)->attr_name + ". Are you sure you're using a model that generates it?";

			Results.push_back(res);
			n_errors++;
			if (test->stop_processing_more_errors)
				break;
		}

		if (rc == 0) {

			// we failed the test
			res.external_rc = test->externl_rc;
			res.internal_rc = test->internal_rc;
			//res.err_msg = test->err_msg + "["+name+"]";
			res.err_msg = test->err_msg; // removed name from message... (previous code in comment above).

			Results.push_back(res);

			if (test->is_warning)
				n_warnings++;
			else
				n_errors++;
			if (test->stop_processing_more_errors)
				break;
		}

		// rc == 1 : nothing to do - passed the test
	}

	if (n_errors > 0)
		return 0;

	return 1;
}

int InputSanityTester::test_if_ok(int pid, long long timestamp,
	const unordered_map<string, unordered_set<string>> &dict_unknown,
	vector<InputSanityTesterResult> &Results) {
	int n_warnings = 0;
	int n_errors = 0;
	for (auto &test : testers) {

		if (test->stage != TESTER_STAGE_BEFORE_MODEL) continue; // these are tested elsewhere

		InputSanityTesterResult res;
		res.external_rc = 0;
		res.internal_rc = 0;
		res.err_msg = "";

		int rc = test->test_if_ok(pid, timestamp, dict_unknown);
		//MLOG("###>>> pid %d time %d test %s : nvals %d nout %d rc %d\n", pid, timestamp, test->tester_params.c_str(), t_nvals, t_noutliers, rc);
		if (rc < 0) {
			res.external_rc = AM_ELIGIBILITY_ERROR;
			res.internal_rc = -2;
			//res.err_msg = "Could not run filter on sample. ["+name+"]";
			if (!test->cant_evel_msg.empty())
				res.err_msg = test->cant_evel_msg;
			else
				res.err_msg = "Could not run filter on sample."; // no name in err message (a request...)

			Results.push_back(res);
			n_errors++;
			if (test->stop_processing_more_errors)
				break;
		}

		if (rc == 0) {

			// we failed the test for a good reason and get out
			res.external_rc = test->externl_rc;
			res.internal_rc = test->internal_rc;
			//res.err_msg = test->err_msg + "["+name+"]";
			res.err_msg = test->err_msg; // no name in err message (a request...)

										 //MLOG("###>>>>>> found an error: %d %d %s\n", res.external_rc, res.internal_rc, res.err_msg.c_str());

			Results.push_back(res);

			if (test->is_warning)
				n_warnings++;
			else
				n_errors++;

			if (test->stop_processing_more_errors)
				break;
		}

		// rc == 1 : nothing to do - passed the test
	}

	//MLOG("###>>> pid %d n_errors %d n_warnings %d\n", pid, n_errors, n_warnings);
	if (n_errors > 0)
		return 0;

	return 1;
}