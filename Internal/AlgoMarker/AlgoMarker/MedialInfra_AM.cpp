#include "AlgoMarker.h"

#include <Logger/Logger/Logger.h>
#include <MedTime/MedTime/MedTime.h>
#include <MedUtils/MedUtils/MedUtils.h>
#include <json/json.hpp>
#include "AlgoMarkerErr.h"

#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL LOG_DEF_LEVEL

#define PREDICTION_SOURCE_UNKNOWN 0
#define PREDICTION_SOURCE_ATTRIBUTE 1
#define PREDICTION_SOURCE_ATTRIBUTE_AS_JSON 2
#define PREDICTION_SOURCE_JSON 3
#define PREDICTION_SOURCE_PREDICTIONS 4
// #define AM_TIMING_LOGS

// Load data stauts
#define LOAD_DATA_STATUS_SUCC 0
#define LOAD_DATA_STATUS_NON_FATAL 1
#define LOAD_DATA_STATUS_FATAL 2

class json_req_export
{
public:
	string field;
	int type = PREDICTION_SOURCE_UNKNOWN;
	int pred_channel = -1; // relevant only if type is PREDICTION_SOURCE_PREDICTIONS

	bool operator==(const json_req_export &other) const
	{
		if (this->field == other.field && this->pred_channel == other.pred_channel &&
			this->type == other.type)
			return true;
		else
			return false;
	}

	struct HashFunction
	{
		size_t operator()(const json_req_export &other) const
		{
			size_t xHash = std::hash<int>()(other.type);
			size_t yHash = std::hash<int>()(other.pred_channel) << 1;
			size_t zHash = std::hash<string>()(other.field) << 2;
			return xHash ^ yHash ^ zHash;
		}
	};
};

class json_req_info
{
public:
	int sample_pid = -1;
	long long sample_time = -1;
	int load_data = 0;
	unordered_map<string, json_req_export> exports;
	string flag_threshold = "";
	float flag_threshold_numeric = MED_MAT_MISSING_VALUE;

	int conv_time = -1;		// this one is calculated
	int sanity_test_rc = 0; // calculated, keeping eligibility testing result
	int sanity_caught_err = 0;
	vector<InputSanityTesterResult> sanity_res;
	MedSample *res = NULL;
};

// local helper functions (these are CalculateByType helpers)
void add_to_json_array(json &js, const string &key, const string &s_add);
void add_to_json_array(nlohmann::ordered_json &js, const string &key, const string &s_add);
void json_to_char_ptr(json &js, char **jarr);
void json_to_char_ptr(nlohmann::ordered_json &js, char **jarr);
bool json_verify_key(json &js, const string &key, int verify_val_flag, const string &val);
bool json_verify_key(nlohmann::ordered_json &js, const string &key, int verify_val_flag, const string &val);
int json_parse_request(json &jreq, json_req_info &defaults, json_req_info &req_i, string &error_message);
void add_flag_response(nlohmann::ordered_json &js, float score, const MedAlgoMarkerInternal &ma,
					   const string &flag_threshold, float flag_threshold_numeric);

//===========================================================================================================
//===========================================================================================================
// MedialInfraAlgoMarker Implementations ::
// Follows is an implementation of an AlgoMarker , which basically means filling in the:
// Load , Unload, ClearData, AddData and Calculate APIs. ( + private internal functions)
// This specific implementation uses medial internal infrastructure for holding data, models, and getting
// predictions.
//===========================================================================================================
//===========================================================================================================
//-----------------------------------------------------------------------------------
// Load() - reading a config file and initializing repository and model
//-----------------------------------------------------------------------------------
int MedialInfraAlgoMarker::Load(const char *config_f)
{
	int rc;

	// read config and check some basic sanities
	rc = read_config(string(config_f));

	if (rc != AM_OK_RC)
		return rc;

	if (type_in_config_file != "MEDIAL_INFRA")
		return AM_ERROR_LOAD_NON_MATCHING_TYPE;

	size_t max_name_len = 5000;
	size_t str_size = strnlen(get_name(), max_name_len);
	if (str_size == 0 || str_size == max_name_len)
	{
		MERR("ERROR: Name is missing\n");
		return AM_ERROR_LOAD_BAD_NAME;
	}

	// loading tester file if needed
	if (input_tester_config_file != "")
	{
		if (ist.read_config(input_tester_config_file) < 0)
		{
			MERR("ERROR: Could not read testers config file %s\n", input_tester_config_file.c_str());
			return AM_ERROR_LOAD_BAD_TESTERS_FILE;
		}
	}

	// prepare internal ma for work: set name, rep and model
	ma.set_name(get_name());
	ma.set_model_end_stage(model_end_stage);

	try
	{
		if (ma.init_rep_config(rep_fname.c_str()) < 0)
			return AM_ERROR_LOAD_READ_REP_ERR;
	}
	catch (...)
	{
		return AM_ERROR_LOAD_READ_REP_ERR;
	}

	ma.set_time_unit_env(get_time_unit());

	try
	{
		if (ma.init_model_from_file(model_fname.c_str()) < 0)
			return AM_ERROR_LOAD_READ_MODEL_ERR;
		if (allow_rep_adjustments)
			ma.fit_model_to_rep();
		if (ma.model_check_required_signals() < 0)
			return AM_ERROR_LOAD_MISSING_REQ_SIGS;
		// if (ma.init_model_for_apply() < 0)
		//	return AM_ERROR_LOAD_READ_MODEL_ERR;
	}
	catch (...)
	{
		return AM_ERROR_LOAD_READ_MODEL_ERR;
	}

	ma.data_load_init();
	// That's it. All is ready for data insert and prediction cycles
	is_loaded = true;
	string vers_info = ma.model_version_info();
	if (vers_info.empty())
		vers_info = "Old model without documented version!";
	MLOG("################ LOADED MODEL VERSION INFO: ##############################\n");
	MLOG("%s\n", vers_info.c_str());
	MLOG("##########################################################################\n");
	return AM_OK_RC;
}

//------------------------------------------------------------------------------------------------
// UnLoad() - clears all data, repository and model, making object ready to be deleted and freed
//------------------------------------------------------------------------------------------------
int MedialInfraAlgoMarker::Unload()
{
	ClearData();
	ma.clear();
	is_loaded = false;
	return AM_OK_RC;
}

//-----------------------------------------------------------------------------------
// ClearData() - clearing current data inserted inside.
//-----------------------------------------------------------------------------------
int MedialInfraAlgoMarker::ClearData()
{
	ma.clear_data();
	return AM_OK_RC;
}

//-----------------------------------------------------------------------------------
// AddData() - adding data for a signal with values and timestamps
//-----------------------------------------------------------------------------------
int MedialInfraAlgoMarker::AddData_data(int patient_id, const char *signalName, int TimeStamps_len, long long *TimeStamps, int Values_len, float *Values,
										map<pair<int, int>, pair<int, vector<char>>> *data)
{
	// At the moment MedialInfraAlgoMarker only loads timestamps given as ints.
	// This may change in the future as needed.
	int *i_times = NULL;
	vector<int> times_int;

	int tu = get_time_unit();
	if (TimeStamps_len > 0)
	{
		times_int.resize(TimeStamps_len);

		// currently assuming we only work with dates ... will have to change this when we'll move to other units
		for (int i = 0; i < TimeStamps_len; i++)
		{
			times_int[i] = AMPoint::auto_time_convert(TimeStamps[i], tu);
			if (times_int[i] < 0)
			{
				MERR("Error in AddData :: patient %d, signals %s, timestamp %lld is ilegal\n",
					 patient_id, signalName, TimeStamps[i]);
				return AM_ERROR_ADD_DATA_FAILED;
			}
			// MLOG("time convert %ld to %d\n", TimeStamps[i], times_int[i]);
		}
		i_times = &times_int[0];
	}

	if (ma.data_load_pid_sig(patient_id, signalName, i_times, TimeStamps_len, Values, Values_len, data) < 0)
		return AM_ERROR_ADD_DATA_FAILED;

	return AM_OK_RC;
}

int MedialInfraAlgoMarker::AddData(int patient_id, const char *signalName, int TimeStamps_len, long long *TimeStamps, int Values_len, float *Values)
{
	MedRepository &rep = ma.get_rep();
	map<pair<int, int>, pair<int, vector<char>>> *data = &rep.in_mem_rep.data;
	return AddData_data(patient_id, signalName, TimeStamps_len, TimeStamps, Values_len, Values, data);
}

//-----------------------------------------------------------------------------------
// AddDatStr() - adding data for a signal with values and timestamps
//-----------------------------------------------------------------------------------
int MedialInfraAlgoMarker::AddDataStr_data(int patient_id, const char *signalName, int TimeStamps_len, long long *TimeStamps, int Values_len, char **Values,
										   map<pair<int, int>, pair<int, vector<char>>> *data)
{
	vector<float> converted_Values;
	vector<long long> final_tm;
	final_tm.reserve(TimeStamps_len);
	converted_Values.reserve(Values_len);
	MedRepository &rep = ma.get_rep();

	try
	{
		string sig = signalName;
		int section_id = rep.dict.section_id(sig);
		int sid = rep.sigs.Name2Sid[sig];

		int Values_i = 0;
		int Time_i = 0;
		const auto &category_map = rep.dict.dict(section_id)->Name2Id;
		int n_elem = 0;
		if (rep.sigs.Sid2Info[sid].n_val_channels > 0)
			n_elem = (int)(Values_len / rep.sigs.Sid2Info[sid].n_val_channels);
		else
			n_elem = (int)(TimeStamps_len / rep.sigs.Sid2Info[sid].n_time_channels);
		for (int i = 0; i < n_elem; i++)
		{
			bool skip_val = false;
			int val_start = Values_i;
			for (int j = 0; j < rep.sigs.Sid2Info[sid].n_val_channels; j++)
			{
				if (rep.sigs.is_categorical_channel(sid, j))
				{
					if (category_map.find(Values[Values_i]) == category_map.end())
					{
						// MWARN("Found undefined code for signal \"%s\" and value \"%s\"\n",
						//	sig.c_str(), Values[Values_i]);
						(*ma.get_unknown_codes(patient_id))[sig].insert(Values[Values_i]);
						skip_val = true;
					}
					++Values_i;
				}
			}
			if (skip_val)
			{
				// remove element!
				Values_len -= rep.sigs.Sid2Info[sid].n_val_channels;
				TimeStamps_len -= rep.sigs.Sid2Info[sid].n_time_channels;
			}
			else
			{
				// All done
				for (int j = 0; j < rep.sigs.Sid2Info[sid].n_time_channels; j++)
				{
					final_tm.push_back(TimeStamps[Time_i]);
					++Time_i;
				}
				for (int j = 0; j < rep.sigs.Sid2Info[sid].n_val_channels; j++)
				{
					float val;
					if (!rep.sigs.is_categorical_channel(sid, j))
						val = stof(Values[Values_i++]);
					else
						val = category_map.at(Values[val_start + j]);
					converted_Values.push_back(val);
				}
			}
		}
	}
	catch (...)
	{
		MERR("Catched Error MedialInfraAlgoMarker::AddDataStr!!\n");
		return AM_FAIL_RC;
	}

	if (TimeStamps_len > 0 || Values_len > 0)
		return AddData_data(patient_id, signalName, TimeStamps_len, final_tm.data(), Values_len, converted_Values.data(), data);
	return AM_OK_RC;
}

int MedialInfraAlgoMarker::AddDataStr(int patient_id, const char *signalName, int TimeStamps_len, long long *TimeStamps, int Values_len, char **Values)
{
	MedRepository &rep = ma.get_rep();
	map<pair<int, int>, pair<int, vector<char>>> *data = &rep.in_mem_rep.data;
	return AddDataStr_data(patient_id, signalName, TimeStamps_len, TimeStamps, Values_len, Values, data);
}

//-----------------------------------------------------------------------------------
// AddDataByType() :
// Supporting loading data directly from a json
//-----------------------------------------------------------------------------------
int MedialInfraAlgoMarker::rec_AddDataByType(int DataType, const char *data, vector<string> &messages)
{
	if (DataType != DATA_JSON_FORMAT && DataType != DATA_BATCH_JSON_FORMAT)
		return AM_ERROR_DATA_UNKNOWN_ADD_DATA_TYPE;

	int ret_code = 0;
	if (DataType == DATA_BATCH_JSON_FORMAT)
	{
		vector<size_t> j_start, j_len;
		vector<char> cdata;
		get_jsons_locations(data, j_start, j_len);
		MedProgress progress("AddDataByType", (int)j_start.size(), 60, 10);
		for (int j = 0; j < j_start.size(); j++)
		{
			if (cdata.size() < j_len[j] + 10)
				cdata.resize(j_len[j] + 10);
			cdata[j_len[j]] = 0;
			strncpy(&cdata[0], &data[j_start[j]], j_len[j]);
			int rc = rec_AddDataByType(DATA_JSON_FORMAT, &cdata[0], messages);
			if (rc != 0)
				ret_code = rc;
			progress.update();
		}
		return ret_code;
	}

	json jsdata;

	try
	{
		jsdata = json::parse(data);
	}
	catch (...)
	{
		char buf[5000];
		snprintf(buf, sizeof(buf), "(%d)Could not parse json data",
				 AM_DATA_BAD_FORMAT_FATAL);
		messages.push_back(string(buf));
		return AM_ERROR_DATA_JSON_PARSE;
	}

	int rc = AddJsonData(-1, jsdata, messages);
	if (rc != 0)
		ret_code = rc;
	return ret_code;
}

int MedialInfraAlgoMarker::AddDataByType(const char *data, char **messages)
{
	vector<string> messages_vec;
	int rc = rec_AddDataByType(DATA_BATCH_JSON_FORMAT, data, messages_vec);
	// transfer messages_vec to messages - each in new line

	if (!messages_vec.empty())
	{
		stringstream ss;
		ss << messages_vec[0];
		for (size_t i = 1; i < messages_vec.size(); ++i)
			ss << "\n"
			   << messages_vec[i];
		string str_ = ss.str();
		*messages = new char[str_.length() + 1];
		(*messages)[str_.length()] = 0;
		strncpy(*messages, str_.data(), str_.length());
	}
	else
		*messages = NULL;

	return rc;
}

string fetch_desription(MedPidRepository &rep, const string &name)
{
	vector<string> tokens;
	boost::split(tokens, name, boost::is_any_of("."));
	if (tokens.size() < 2)
		return "";
	int section_id = rep.dict.section_id(tokens[0]);
	stringstream ss_val;
	ss_val << tokens[1];
	for (size_t i = 2; i < tokens.size(); ++i) // if has more than 1 dot, as part of the code
		ss_val << "." << tokens[i];
	string curr_name = ss_val.str();

	int sig_val = rep.dict.id(section_id, curr_name);
	if (rep.dict.name(section_id, sig_val).empty())
		return "";
	const vector<string> &names = rep.dict.dicts[section_id].Id2Names.at(sig_val);
	if (names.size() <= 1) // has only this "name", no other aliases
		return "";
	stringstream ss;
	bool is_empty = true;
	for (size_t i = 0; i < names.size(); ++i)
	{
		if (names[i] == curr_name)
			continue;
		if (!is_empty)
			ss << "|";
		ss << names[i];
		is_empty = false;
	}
	return ss.str();
}

void process_signal_json(MedPidRepository &rep, Explainer_record_config &e_cfg,
						 int pid, int time, nlohmann::ordered_json &jres)
{
	UniversalSigVec usv;
	int sid = rep.sigs.sid(e_cfg.signal_name);
	if (sid < 0)
		MTHROW_AND_ERR("Error unknown signal %s\n", e_cfg.signal_name.c_str());
	rep.uget(pid, sid, usv);

	bool test_sets = false;
	int section_id = rep.dict.section_id(e_cfg.signal_name);
	if (rep.sigs.is_categorical_channel(sid, e_cfg.val_channel) && !e_cfg.sets.empty())
	{
		// init lut if needed in first time
		if (e_cfg.lut.empty())
			rep.dict.prep_sets_lookup_table(section_id, e_cfg.sets, e_cfg.lut);
		test_sets = true;
	}

	int min_time = med_time_converter.add_subtract_time(time, usv.time_unit(),
														-e_cfg.max_time_window, e_cfg.time_unit);
	int cnt = 0;
	for (int i = usv.len - 1; i >= 0; --i)
	{
		if (usv.Time(i, e_cfg.time_channel) > time)
			continue;
		if (usv.Time(i, e_cfg.time_channel) < min_time)
			break;
		// test categorical values if needed
		if (test_sets && !e_cfg.lut[usv.Val(i, e_cfg.val_channel)])
			continue;

		// In time window - take at most e_cfg.max_count records, most recent.
		nlohmann::ordered_json ele;
		ele["signal"] = e_cfg.signal_name;
		string full_unit = rep.sigs.unit_of_measurement(sid, 0);
		ele["unit"] = nlohmann::ordered_json::array();
		for (int j = 0; j < usv.n_val_channels(); ++j)
			ele["unit"].push_back(rep.sigs.unit_of_measurement(sid, j));

		ele["timestamp"] = nlohmann::ordered_json::array();
		ele["value"] = nlohmann::ordered_json::array();
		// Print time values:
		for (int j = 0; j < usv.n_time_channels(); ++j)
			ele["timestamp"].push_back(usv.Time(i, j));
		for (int j = 0; j < usv.n_val_channels(); ++j)
		{

			if (!rep.sigs.is_categorical_channel(sid, j))
				ele["value"].push_back(to_string(usv.Val(i, j)));
			else
			{
				int code_val = usv.Val<int>(i, j);
				string code_str = rep.dict.name(section_id, code_val);
				ele["value"].push_back(code_str);
			}
		}

		jres.push_back(ele);
		++cnt;
		if (cnt >= e_cfg.max_count)
			break;
	}
}

void process_explainability(nlohmann::ordered_json &jattr,
							Explainer_parameters &ex_params, MedPidRepository &rep,
							int pid, int time)
{
	if (jattr.find("explainer_output") != jattr.end())
	{
		auto final_res = nlohmann::ordered_json::array();
		int total_reasons = 0, tot_pos = 0, tot_neg = 0;
		for (auto &e : jattr["explainer_output"])
		{

			if (e.find("contributor_name") == e.end())
				MTHROW_AND_ERR("Error - expecting to see contributor_name\n");
			string contrib_name = e["contributor_name"].get<string>();
			if (ex_params.ignore_groups_list.find(contrib_name) != ex_params.ignore_groups_list.end())
			{
				// remove group - in ignore, remove "e" from jattr["explainer_output"]
				continue;
			}
			// Add "contributor_description" from name:
			e["contributor_description"] = fetch_desription(rep, contrib_name);

			// test filters:
			bool contrib_positive = true;
			if (e.find("contributor_value") != e.end())
			{
				if (ex_params.threshold_abs > 0)
				{
					if (abs(e["contributor_value"].get<float>()) < ex_params.threshold_abs)
						continue; // sorted, can change to break
				}
				if (e["contributor_value"].get<float>() < 0)
					contrib_positive = false;
				if (contrib_positive && ex_params.total_max_pos_reasons >= 0 && tot_pos >= ex_params.total_max_pos_reasons)
					continue;
				if (!contrib_positive && ex_params.total_max_neg_reasons >= 0 && tot_neg >= ex_params.total_max_pos_reasons)
					continue;
			}
			if (e.find("contributor_percentage") != e.end())
			{
				if (ex_params.threshold_percentage > 0)
				{
					if (abs(e["contributor_percentage"].get<float>()) < ex_params.threshold_percentage)
						continue; // sorted, can change to break
				}
			}
			if (ex_params.total_max_reasons >= 0 && total_reasons >= ex_params.total_max_reasons)
				break;
			// after all filters:
			++total_reasons;
			if (e.find("contributor_value") != e.end())
			{
				if (contrib_positive)
					++tot_pos;
				else
					++tot_neg;
			}

			if (e.find("contributor_value") != e.end() && ex_params.max_threshold > 0 && ex_params.num_groups > 0)
			{
				float level_bin;
				if (ex_params.use_perc)
					level_bin = (abs(e["contributor_percentage"].get<float>()) / ex_params.max_threshold);
				else
					level_bin = (abs(e["contributor_value"].get<float>()) / ex_params.max_threshold);
				if (level_bin > 1)
					level_bin = 1;

				level_bin *= ex_params.num_groups;
				if (level_bin > 0)
					level_bin = (int)(level_bin) + 1;
				else
					level_bin = int(level_bin);

				if (level_bin > ex_params.num_groups)
					level_bin = ex_params.num_groups;

				e["contributor_level"] = (int)level_bin;
				e["contributor_level_max"] = ex_params.num_groups;
			}

			string contib_info = "";
			e["contributor_records"] = nlohmann::ordered_json::array();
			if (ex_params.cfg.records.find(contrib_name) != ex_params.cfg.records.end())
			{
				Explainer_record_config &e_cfg = ex_params.cfg.records.at(contrib_name);
				e["contributor_records_info"] = nlohmann::ordered_json();
				e["contributor_records_info"]["contributor_max_time"] = e_cfg.max_time_window;
				e["contributor_records_info"]["contributor_max_time_unit"] = med_time_converter.type_to_string(e_cfg.time_unit);
				e["contributor_records_info"]["contributor_max_count"] = e_cfg.max_count;
				process_signal_json(rep, e_cfg, pid, time, e["contributor_records"]);
			}
			else
			{
				// check if Age, or if has 1 feature, so present it:
				if (e.find("contributor_elements") != e.end() && e["contributor_elements"].size() == 1)
				{
					string fname = e["contributor_elements"].begin()->at("feature_name").get<string>();
					float fval = e["contributor_elements"].begin()->at("feature_value").get<float>();

					nlohmann::ordered_json element_single;
					element_single["signal"] = fname;
					element_single["unit"] = nlohmann::ordered_json::array();
					if (fname == "Age")
						element_single["unit"].push_back("Year");
					else
					{
						element_single["unit"].push_back("");
					}
					element_single["timestamp"] = nlohmann::ordered_json::array();
					element_single["value"] = nlohmann::ordered_json::array();
					element_single["value"].push_back(to_string(fval));
					e["contributor_records"].push_back(element_single);
				}
			}

			final_res.push_back(e);
		}

		// add static info:
		if (!ex_params.static_features_info.empty())
		{
			// special use cage Age, other names is to fetch signal:
			jattr["static_info"] = nlohmann::ordered_json::array();
			for (const string &feat : ex_params.static_features_info)
			{
				nlohmann::ordered_json feat_js;
				feat_js["signal"] = feat;
				if (boost::to_upper_copy(feat) == "AGE")
				{
					feat_js["unit"] = "Year";
					int sid = rep.sigs.sid("BDATE");
					bool using_byear = false;
					if (sid < 0)
					{
						sid = rep.sigs.sid("BYEAR");
						using_byear = true;
						if (sid < 0)
						{
							MWARN("Error unknown signal BDATE/BYEAR for Age static fetch\n");
							feat_js["value"] = "Missing";
							jattr["static_info"].push_back(feat_js);
							continue;
						}
					}
					int bdate = medial::repository::get_value(rep, pid, sid);
					int byear = bdate;
					if (!using_byear)
						byear = int(bdate / 10000);
					if (bdate < 0)
						feat_js["value"] = "Missing";
					else
					{
						int age = int(time / 10000) - byear;
						feat_js["value"] = to_string(age);
					}
				}
				else
				{
					UniversalSigVec usv;
					int sid = rep.sigs.sid(feat);
					int section_id = rep.dict.section_id(feat);
					if (sid < 0)
						MTHROW_AND_ERR("Error unknown signal %s for static fetch\n", feat.c_str());
					rep.uget(pid, sid, usv);
					feat_js["unit"] = rep.sigs.unit_of_measurement(sid, 0);
					if (usv.len == 0)
						feat_js["value"] = "Missing";
					else
					{
						if (!rep.sigs.is_categorical_channel(sid, 0))
							feat_js["value"] = to_string(usv.Val(usv.len - 1, 0));
						else
						{
							int code_val = usv.Val<int>(usv.len - 1, 0);
							string code_str = rep.dict.name(section_id, code_val);
							feat_js["value"] = code_str;
						}
					}
				}

				jattr["static_info"].push_back(feat_js);
			}
		}

		jattr.erase("explainer_output");
		jattr["explainer_output"] = final_res;
	}
}

//------------------------------------------------------------------------------------------
// Calculate() - after data loading : get a request, get predictions, and pack as responses
//------------------------------------------------------------------------------------------
int MedialInfraAlgoMarker::Calculate(AMRequest *request, AMResponses *responses)
{
	MWARN("Warning : Calculate is deprecated and will not be supported in the future, please use CalculateByType\n");
#ifdef AM_TIMING_LOGS
	MedTimer timer;
	timer.start();
#endif
	if (sort_needed)
	{
		if (ma.data_load_end() < 0)
			return AM_FAIL_RC;
	}
#ifdef AM_TIMING_LOGS
	timer.take_curr_time();
	MLOG("INFO:: MedialInfraAlgoMarker::Calculate :: data_load_end %2.1f milisecond\n", timer.diff_milisec());
#endif

	if (!ma.model_initiated())
	{
#ifdef AM_TIMING_LOGS
		timer.start();
#endif
		if (ma.init_model_for_apply() < 0)
			return AM_FAIL_RC;
#ifdef AM_TIMING_LOGS
		timer.take_curr_time();
		MLOG("INFO:: MedialInfraAlgoMarker::Calculate :: init_model_for_apply %2.1f milisecond\n", timer.diff_milisec());
#endif
	}

	if (responses == NULL)
		return AM_FAIL_RC;

#ifdef AM_TIMING_LOGS
	timer.start();
#endif
	AMMessages *shared_msgs = responses->get_shared_messages();
#ifdef AM_TIMING_LOGS
	timer.take_curr_time();
	MLOG("INFO:: MedialInfraAlgoMarker::Calculate :: get_shared_messages %2.1f milisecond\n", timer.diff_milisec());
#endif

	if (request == NULL)
	{
		string msg = "Error :: (" + to_string(AM_MSG_NULL_REQUEST) + " ) NULL request in Calculate()";
		shared_msgs->insert_message(AM_GENERAL_FATAL, msg.c_str());
		return AM_FAIL_RC;
	}

	// string msg_prefix = "reqId: " + string(request->get_request_id()) + " :: ";
	string msg_prefix = ""; // asked not to put reqId in messages.... (not sure it's a good idea, prev code above in comment)
	responses->set_request_id(request->get_request_id());

#ifdef AM_TIMING_LOGS
	timer.start();
#endif
	for (int i = 0; i < request->get_n_score_types(); i++)
	{
		char *stype = request->get_score_type(i);
		responses->insert_score_types(&stype, 1);
	}
#ifdef AM_TIMING_LOGS
	timer.take_curr_time();
	MLOG("INFO:: MedialInfraAlgoMarker::Calculate :: get_score_type %2.1f milisecond\n", timer.diff_milisec());
#endif

	// We now have to prepare samples for the requested points
	// again - we only deal with int times in this class, so we convert the long long stamps to int
#ifdef AM_TIMING_LOGS
	timer.start();
#endif
	ma.clear_samples();
	int n_points = request->get_n_points();
	int tu = get_time_unit();
	vector<int> conv_times;

	for (int i = 0; i < n_points; i++)
	{
		conv_times.push_back(AMPoint::auto_time_convert(request->get_timestamp(i), tu));
		if ((ma.insert_sample(request->get_pid(i), conv_times.back()) < 0) || (conv_times.back() <= 0))
		{
			string msg = msg_prefix + "(" + to_string(AM_MSG_BAD_PREDICTION_POINT) + ") Failed insert prediction point " + to_string(i) + " pid: " + to_string(request->get_pid(i)) + " ts: " + to_string(request->get_timestamp(i));
			shared_msgs->insert_message(AM_GENERAL_FATAL, msg.c_str());
			return AM_FAIL_RC;
		}
	}

	ma.normalize_samples();
#ifdef AM_TIMING_LOGS
	timer.take_curr_time();
	MLOG("INFO:: MedialInfraAlgoMarker::Calculate :: prepared_samples %2.1f milisecond\n", timer.diff_milisec());

	timer.start();
#endif
	// Checking score types and verify they are supported
	int n_score_types = request->get_n_score_types();
	for (int i = 0; i < n_score_types; i++)
	{
		if (!IsScoreTypeSupported(request->get_score_type(i)))
		{
			// string msg = msg_prefix + "(" + to_string(AM_MSG_BAD_SCORE_TYPE) + ") AlgoMarker of type " + string(get_name()) + " does not support score type " + string(request->get_score_type(i));
			string msg = msg_prefix + "AlgoMarker of type " + string(get_name()) + " does not support score type " + string(request->get_score_type(i));
			shared_msgs->insert_message(AM_GENERAL_FATAL, msg.c_str());
			return AM_FAIL_RC;
		}
	}

	// At this stage we need to create a response entry for each of the requested points
	// Then we have to test for eligibility - err the ones that are not eligible
	// And then score all the eligible ones in a single batch.
	vector<int> eligible_pids, eligible_timepoints;
	vector<long long> eligible_ts;
	MedPidRepository &rep = ma.get_rep();
	unordered_map<unsigned long long, vector<long long>> sample2ts; // conversion of each sample to all the ts that were mapped to it.

	int n_bad_scores = 0;
	for (int i = 0; i < n_points; i++)
	{
		int _pid = request->get_pid(i);
		long long _ts = request->get_timestamp(i);

		// create a response
		AMResponse *res = responses->create_point_response(_pid, _ts);
		// AMResponse *res = responses->get_response_by_point(_pid, (long long)conv_times[i]);
		//		if (res == NULL)
		//			res = responses->create_point_response(_pid, _ts);
		//		res = responses->create_point_response(_pid, (long long)conv_times[i]);

		// test this point for eligibility and add errors if needed
		vector<InputSanityTesterResult> test_res;
		int test_rc = ist.test_if_ok(_pid, (long long)conv_times[i], *ma.get_unknown_codes(_pid), test_res);
		int test_rc2 = ist.test_if_ok(rep, _pid, (long long)conv_times[i], test_res);
		if (test_rc2 < 1)
			test_rc = test_rc2;

		// push messages if there are any
		AMMessages *msgs = res->get_msgs();
		for (auto &tres : test_res)
		{
			// string msg = msg_prefix + tres.err_msg + " Internal Code: " + to_string(tres.internal_rc);
			string msg = msg_prefix + tres.err_msg; // messages without Internal codes...(prev code in comment above).
			msgs->insert_message(tres.external_rc, msg.c_str());
		}

		if (test_rc <= 0)
		{
			n_bad_scores++;
		}
		else
		{
			// MLOG("DEBUG ===> i %d _pid %d conv %d _ts %lld size %d\n", i, _pid, conv_times[i], _ts, eligible_pids.size());
			eligible_pids.push_back(_pid);
			eligible_timepoints.push_back(conv_times[i]);
			eligible_ts.push_back(_ts);
			unsigned long long p = ((unsigned long long)_pid << 32) | (conv_times[i]);
			if (sample2ts.find(p) == sample2ts.end())
				sample2ts[p] = vector<long long>();
			sample2ts[p].push_back(_ts);
		}
	}

	int _n_points = (int)eligible_pids.size();

#ifdef AM_TIMING_LOGS
	timer.take_curr_time();
	MLOG("INFO:: MedialInfraAlgoMarker::Calculate :: tested_eligibility %2.1f milisecond. Has %d samples\n", timer.diff_milisec(), _n_points);

	ma.model_apply_verbose(true);
	timer.start();
#endif

	// Calculating raw scores for eligble points
	vector<float> raw_scores(_n_points, (float)AM_UNDEFINED_VALUE);
	int get_preds_rc = -1;
	try
	{
		if ((get_preds_rc = ma.get_preds(&eligible_pids[0], &eligible_timepoints[0], &raw_scores[0], _n_points)) < 0)
		{
			string msg = msg_prefix + "(" + to_string(AM_MSG_RAW_SCORES_ERROR) + ") Failed getting RAW scores in AlgoMarker " + string(get_name()) + " With return code " + to_string(get_preds_rc);
			shared_msgs->insert_message(AM_GENERAL_FATAL, msg.c_str());
			return AM_FAIL_RC;
		}
	}
	catch (...)
	{
		string msg = msg_prefix + "(" + to_string(AM_MSG_RAW_SCORES_ERROR) + ") Failed getting RAW scores in AlgoMarker " + string(get_name()) + " caught a crash. ";
		shared_msgs->insert_message(AM_GENERAL_FATAL, msg.c_str());
		return AM_FAIL_RC;
	}

	if (am_matrix != "" && _n_points > 0)
	{ // debug only
		if (first_write)
			ma.write_features_mat(am_matrix);
		else
			ma.add_features_mat(am_matrix);
		first_write = false;
	}

#ifdef AM_TIMING_LOGS
	timer.take_curr_time();
	MLOG("INFO:: MedialInfraAlgoMarker::Calculate :: get_preds %2.1f milisecond\n", timer.diff_milisec());
	timer.start();
#endif
	// going over scores, and adding them to the right responses
	char **_score_types;
	int _n_score_types;
	responses->get_score_types(&_n_score_types, &_score_types);

	MedSamples *meds = ma.get_samples_ptr();
	Explainer_parameters ex_params;
	ma.get_explainer_params(ex_params);

	for (auto &id_s : meds->idSamples)
	{
		for (auto &s : id_s.samples)
		{

			// basic info for current sample
			int c_pid = s.id;
			int c_ts = s.time;
			float c_scr = s.prediction.size() > 0 ? s.prediction[0] : (float)AM_UNDEFINED_VALUE;
			string c_ext_scr = "";
			if (s.str_attributes.size() > 0)
			{
				nlohmann::ordered_json c_ext_scr_json({});

				for (auto &ex_res_field_name : extended_result_fields)
				{
					if (s.str_attributes.count(ex_res_field_name))
					{
						c_ext_scr_json[ex_res_field_name] = nlohmann::ordered_json::parse(s.str_attributes[ex_res_field_name]);
						process_explainability(c_ext_scr_json[ex_res_field_name], ex_params,
											   rep, c_pid, c_ts);
					}
				}
				c_ext_scr = c_ext_scr_json.dump();
			}

			unsigned long long p = ((unsigned long long)c_pid << 32) | (c_ts);

			for (auto ts : sample2ts[p])
			{

				// DEBUG
				// for (auto &attr : s.attributes) MLOG("pid %d time %d score %f attr %s %f\n", c_pid, c_ts, c_scr, attr.first.c_str(), attr.second);

				// get the matching response (should be prepared already)
				AMResponse *res = responses->get_response_by_point(c_pid, ts);

				if (res != NULL)
				{

					// we now test the attribute tests
					vector<InputSanityTesterResult> test_res;
					int test_rc = ist.test_if_ok(s, test_res);

					AMMessages *msgs = res->get_msgs();
					for (auto &tres : test_res)
					{
						// string msg = msg_prefix + tres.err_msg + " Internal Code: " + to_string(tres.internal_rc);
						string msg = msg_prefix + tres.err_msg; // no Internal Code message (prev code in comment above).
						// MLOG("Inserting attr error to pid %d ts %d : %d : %s\n", c_pid, ts, tres.external_rc, msg.c_str());
						msgs->insert_message(tres.external_rc, msg.c_str());
					}

					if (test_rc <= 0)
						n_bad_scores++;
					else
					{

						// all is fine, we insert the score into its place
						res->init_scores(_n_score_types);

						for (int j = 0; j < _n_score_types; j++)
						{

							if (strcmp(_score_types[j], "Raw") == 0)
							{
								res->set_score(j, c_scr, _score_types[j], c_ext_scr);
							}
							else
							{
								res->set_score(j, (float)AM_UNDEFINED_VALUE, _score_types[j], "");
								AMScore *am_scr = res->get_am_score(j);
								AMMessages *msgs = am_scr->get_msgs();
								string msg = msg_prefix + "Undefined Score Type: " + string(_score_types[j]);
								msgs->insert_message(AM_GENERAL_FATAL, msg.c_str());
							}
						}
					}
				}
			}
		}
	}

#ifdef AM_TIMING_LOGS
	timer.take_curr_time();
	MLOG("INFO:: MedialInfraAlgoMarker::Calculate :: finished_response %2.1f milisecond\n", timer.diff_milisec());
#endif

	if (n_bad_scores > 0)
	{
		string msg = msg_prefix + "Failed input tests for " + to_string(n_bad_scores) + " out of " + to_string(n_points) + " scores";
		if (n_bad_scores < n_points)
		{
			shared_msgs->insert_message(AM_RESPONSES_ELIGIBILITY_ERROR, msg.c_str());
			return AM_OK_RC;
		}

		shared_msgs->insert_message(AM_RESPONSES_ELIGIBILITY_ERROR, msg.c_str());
		return AM_FAIL_RC;
	}

	return AM_OK_RC;
}

void MedialInfraAlgoMarker::clear_patients_data(const vector<int> &pids)
{
	if (pids.empty())
		return;

	try
	{
		for (int pid : pids)
			ma.get_rep().in_mem_rep.erase_pid_data(pid);
	}
	catch (...)
	{
		MERR("Error in clear_patients_data for %zu patients\n", pids.size());
	}
}

int get_msg_status(const string &msg)
{
	int code = 0;
	int end_idx = -1;
	for (size_t i = 1; i < msg.length(); ++i)
	{
		if (msg[i] == ')')
		{
			end_idx = i;
			break;
		}
	}
	if (end_idx > 0)
	{
		try
		{
			code = stoi(msg.substr(1, end_idx - 1));
		}
		catch (...)
		{
			code = AM_DATA_BAD_FORMAT_FATAL;
		}
	}

	return code;
}

bool is_fatal_load_message(const vector<string> &messages)
{
	bool fatal = false;
	for (const string &msg : messages)
	{
		int message_status = get_msg_status(msg);
		if (message_status == AM_DATA_BAD_FORMAT_FATAL)
		{
			fatal = true;
			break;
		}
	}

	return fatal;
}

//------------------------------------------------------------------------------------------
// CalculateByType : alllows for a general json in -> json out API with many more options
//------------------------------------------------------------------------------------------
int MedialInfraAlgoMarker::CalculateByType(int CalculateType, char *request, char **response)
{
	if (CalculateType != JSON_REQ_JSON_RESP)
		return AM_FAIL_RC;

#ifdef AM_TIMING_LOGS
	ma.model_apply_verbose(true);
	MedTimer timer;
	timer.start();
#endif
	if (sort_needed)
	{
		if (ma.data_load_end() < 0)
			return AM_FAIL_RC;
	}
#ifdef AM_TIMING_LOGS
	timer.take_curr_time();
	MLOG("INFO:: MedialInfraAlgoMarker::CalculateByType :: data_load_end %2.1f milisecond\n", timer.diff_milisec());
	timer.start();
#endif

	if (!ma.model_initiated())
	{
		if (ma.init_model_for_apply() < 0)
			return AM_FAIL_RC;
	}
#ifdef AM_TIMING_LOGS
	timer.take_curr_time();
	MLOG("INFO:: MedialInfraAlgoMarker::CalculateByType :: init_model_for_apply %2.1f milisecond\n", timer.diff_milisec());
	timer.start();
#endif

	json jreq;
	nlohmann::ordered_json jresp;

	jresp = nlohmann::ordered_json({{"type", "response"}});
	try
	{
		jreq = json::parse(request);
	}
	catch (...)
	{
		add_to_json_array(jresp, "errors", "ERROR: Could not parse request as a valid json");
		json_to_char_ptr(jresp, response);
		return AM_FAIL_RC;
	}

	// verify the "type" : "request" , and the "request_id" : something fields
	if (!json_verify_key(jreq, "type", 1, "request"))
		add_to_json_array(jresp, "errors", "ERROR: missing type request");
	string request_id;
	if (!json_verify_key(jreq, "request_id", 0, ""))
		add_to_json_array(jresp, "errors", "ERROR: no request_id provided");
	else
	{
		request_id = jreq["request_id"].get<string>();
		jresp.push_back({"request_id", request_id});
	}

	if (!json_verify_key(jreq, "requests", 0, ""))
		add_to_json_array(jresp, "errors", "ERROR: missing actual requests in request " + request_id);

	if (json_verify_key(jresp, "errors", 0, ""))
	{
		json_to_char_ptr(jresp, response);
		return AM_FAIL_RC;
	} // Leave now if there are errors

	// default parameters
	json_req_info defaults;

	vector<json_req_info> sample_reqs;

	//	try {
	string error_message_in_req;
	if (json_parse_request(jreq, defaults, defaults, error_message_in_req) != 0)
	{
		add_to_json_array(jresp, "errors", "ERROR: general json error in parsing the default request fields in request id " + request_id + " " + error_message_in_req);
		json_to_char_ptr(jresp, response);
		return AM_FAIL_RC;
	}

#ifdef AM_TIMING_LOGS
	timer.take_curr_time();
	MLOG("INFO:: MedialInfraAlgoMarker::CalculateByType :: json_parse_request %2.1f milisecond\n", timer.diff_milisec());
	timer.start();
#endif

	bool has_load_data = false;
	MedPidRepository *rep = &ma.get_rep();
	map<pair<int, int>, pair<int, vector<char>>> new_data;
	vector<vector<string>> all_load_data_msgs(jreq["requests"].size());
	int req_id = 0;
	for (auto &jreq_i : jreq["requests"])
	{
		json_req_info j_i;
		if (json_parse_request(jreq_i, defaults, j_i, error_message_in_req) != 0)
		{
			add_to_json_array(jresp, "errors", "ERROR: general json error in parsing the request fields in request id " + request_id + " " + error_message_in_req);
			json_to_char_ptr(jresp, response);
			return AM_FAIL_RC;
		}
		sample_reqs.push_back(j_i);
		vector<string> &vec_messages = all_load_data_msgs[req_id];
		if (j_i.load_data && json_verify_key(jreq_i, "data", 0, ""))
		{
			has_load_data = true;
			int load_status = LOAD_DATA_STATUS_SUCC;
			if (AddJsonData(j_i.sample_pid, jreq_i["data"], vec_messages, &new_data) != AM_OK_RC)
			{
				load_status = LOAD_DATA_STATUS_NON_FATAL;
				if (is_fatal_load_message(vec_messages))
					load_status = LOAD_DATA_STATUS_FATAL;
			}
			if (load_status == LOAD_DATA_STATUS_FATAL)
			{
				add_to_json_array(jresp, "errors", "ERROR: error when loading data for patient id " + to_string(j_i.sample_pid));
				// Store error message anyway if has (also if successed)
				for (size_t i = 0; i < vec_messages.size(); ++i)
					add_to_json_array(jresp, "errors", vec_messages[i]);
			}
		}
		++req_id;
	}
	MedPidRepository rep_copy;
	if (has_load_data)
	{
		rep_copy.sigs = rep->sigs; // Created a copy - TODO - smarter copy of all data by pointers, except data
		rep_copy.dict = rep->dict;
		rep_copy.switch_to_in_mem_mode();
		rep_copy.in_mem_rep.data = move(new_data);
		rep_copy.in_mem_rep.sortData();
		rep = &rep_copy;
		// Debug print new_data:
		// for (auto &it : rep->in_mem_rep.data) {
		//	MLOG("pid: %d, signal: %s, size: %zu\n", it.first.first, rep->sigs.name(it.first.second).c_str(),
		//	it.second.second.size());
		// }

		// Test get BDATE
		// UniversalSigVec usv_gender;
		// rep->uget(1, "BDATE", usv_gender);
		// MLOG("Will print Gender, size: %d, in_mem=%d\n", usv_gender.len, int(rep->in_mem_mode_active()));
		// for (size_t ii=0; ii< usv_gender.len; ++ii) {
		//	MLOG("Gender = %f\n", usv_gender.Val(ii));
		// }
	}
	//	}

	//	catch (...) {
	//		add_to_json_array(jresp, "errors", "ERROR: error when parsing requests (or loading)");
	//	}
	if (json_verify_key(jresp, "errors", 0, ""))
	{
		json_to_char_ptr(jresp, response);
		return AM_FAIL_RC;
	} // Leave now if there are errors

#ifdef AM_TIMING_LOGS
	timer.take_curr_time();
	MLOG("INFO:: MedialInfraAlgoMarker::CalculateByType :: end load_data %2.1f milisecond\n", timer.diff_milisec());
	timer.start();
#endif

	// We now convert times and do an initial sanity checks
	// again - we only deal with int times in this class, so we convert the long long stamps to int
	// we also run the eligibility tests, keep the results, and make lists of all eligible points for scoring.
	int n_points = (int)sample_reqs.size();
	int tu = get_time_unit();

	vector<int> eligible_pids, eligible_timepoints;
	unordered_map<unsigned long long, vector<int>> sample2ind; // conversion of each sample to all the ts that were mapped to it.
	int n_failed = 0, n_bad = 0;
	// #pragma omp parallel for if (n_points > 10)
	for (int i = 0; i < n_points; i++)
	{
		json_req_info &req_i = sample_reqs[i];
		req_i.conv_time = AMPoint::auto_time_convert(req_i.sample_time, tu);
		int ok_time = 1;
		if (tu == MedTime::Date && (req_i.conv_time < 19500000 || req_i.conv_time > 30000000))
			ok_time = 0;
		if ((req_i.sample_pid <= 0) || (req_i.conv_time <= 0) || ok_time == 0)
		{
#pragma omp critical
			{
				add_to_json_array(jresp, "errors", "ERROR: BAD request patient id or time : failed in inserting pid: " + to_string(req_i.sample_pid) + " ts: " + to_string(req_i.sample_time));
				req_i.sanity_test_rc = -2;
				n_failed++;
			}
		}

		else
		{
			try
			{
				req_i.sanity_test_rc = ist.test_if_ok(req_i.sample_pid, (long long)req_i.conv_time, *ma.get_unknown_codes(req_i.sample_pid), req_i.sanity_res);
				int rc_res = ist.test_if_ok(*rep, req_i.sample_pid, (long long)req_i.conv_time, req_i.sanity_res);
				if (rc_res < 1)
					req_i.sanity_test_rc = rc_res;
			}
			catch (...)
			{
				req_i.sanity_caught_err = 1;
			}
		}

#pragma omp critical
		if (req_i.sanity_caught_err == 0 && req_i.sanity_test_rc > 0)
		{
			unsigned long long p = ((unsigned long long)req_i.sample_pid << 32) | req_i.conv_time;
			if (sample2ind.find(p) == sample2ind.end())
			{
				eligible_pids.push_back(req_i.sample_pid);
				eligible_timepoints.push_back(req_i.conv_time);
				sample2ind[p] = vector<int>();
			}
			sample2ind[p].push_back(i);
		}
		else
			n_bad++;
	}

	if (n_failed > 0)
	{
		json_to_char_ptr(jresp, response);
		return AM_FAIL_RC;
	}

#ifdef AM_TIMING_LOGS
	timer.take_curr_time();
	MLOG("INFO:: MedialInfraAlgoMarker::CalculateByType :: end calc_eligiblilty %2.1f milisecond\n", timer.diff_milisec());
	timer.start();
#endif

	// Fetch what's needed:
	unordered_set<json_req_export, json_req_export::HashFunction> set_rep_f;
	for (int i = 0; i < sample_reqs.size(); i++)
	{
		json_req_info &req_i = sample_reqs[i];
		// set req_i.flag_threshold if empty:
		if (req_i.flag_threshold.empty())
			req_i.flag_threshold = ma.get_default_threshold();
		unordered_map<string, json_req_export> &req_fields = req_i.exports;
		for (auto &it : req_fields)
			set_rep_f.insert(it.second);
	}
	vector<json_req_export> uniq_req_fields(set_rep_f.begin(), set_rep_f.end());
	// convert to requested_fields
	unordered_set<Effected_Field, Effected_Field::HashFunction> set_requested_fields;
	for (size_t i = 0; i < uniq_req_fields.size(); ++i)
	{
		const json_req_export &jinp = uniq_req_fields[i];
		Effected_Field f;
		switch (jinp.type)
		{
		case PREDICTION_SOURCE_ATTRIBUTE:
			f.field = Field_Type::NUMERIC_ATTRIBUTE;
			f.value_name = jinp.field;
			// add also string attributes
			set_requested_fields.insert(Effected_Field(Field_Type::STRING_ATTRIBUTE, jinp.field));
			break;
		case PREDICTION_SOURCE_ATTRIBUTE_AS_JSON:
			f.field = Field_Type::STRING_ATTRIBUTE;
			f.value_name = jinp.field;
			break;
		case PREDICTION_SOURCE_JSON:
			f.field = Field_Type::JSON_DATA;
			f.value_name = jinp.field;
			break;
		case PREDICTION_SOURCE_PREDICTIONS:
			f.field = Field_Type::PREDICTION;
			f.value_name = to_string(jinp.pred_channel);
			break;
		default:
			MLOG("WARN unknown request field %s\n", jinp.field.c_str());
			break;
		}
		set_requested_fields.insert(f);
	}
	vector<Effected_Field> requested_fields(set_requested_fields.begin(), set_requested_fields.end());

	// at this point in time we are ready to score eligible_pids,eligible_timepoints. We will do that, and later wrap it all up into a single json back.
	int _n_points = (int)eligible_pids.size();
	int get_preds_rc = -1;
	try
	{
		if ((get_preds_rc = ma.get_preds(&eligible_pids[0], &eligible_timepoints[0], NULL, _n_points, requested_fields, rep)) < 0)
		{
			add_to_json_array(jresp, "errors", "ERROR: (" + to_string(AM_MSG_RAW_SCORES_ERROR) + ") Failed getting scores in AlgoMarker " + string(get_name()) + " With return code " + to_string(get_preds_rc));
			json_to_char_ptr(jresp, response);
			return AM_FAIL_RC;
		}
	}
	catch (...)
	{
		add_to_json_array(jresp, "errors", "ERROR: (" + to_string(AM_MSG_RAW_SCORES_ERROR) + ") Failed getting scores in AlgoMarker " + string(get_name()) + " caught a crash");
		json_to_char_ptr(jresp, response);
		return AM_FAIL_RC;
	}

#ifdef AM_TIMING_LOGS
	timer.take_curr_time();
	MLOG("INFO:: MedialInfraAlgoMarker::CalculateByType :: end get_preds %2.1f milisecond\n", timer.diff_milisec());
	timer.start();
#endif

	// now we are ready ... we have the results, and we need to put it all into the response json one by one.

	// adding result samples pointers to sample_reqs
	MedSamples *samps = ma.get_samples_ptr();
	for (auto &ids : samps->idSamples)
		for (auto &s : ids.samples)
		{
			unsigned long long p = ((unsigned long long)s.id << 32) | s.time;
			for (auto j : sample2ind[p])
				sample_reqs[j].res = &s;
		}

	jresp.push_back({"responses", json::array()});
	for (int i = 0; i < sample_reqs.size(); i++)
	{
		json_req_info &req_i = sample_reqs[i];
		if (req_i.res != NULL)
		{

			try
			{
				int test_rc = ist.test_if_ok(*(req_i.res), req_i.sanity_res);
				if (test_rc == 0)
					req_i.sanity_test_rc = 0;
			}
			catch (...)
			{
				req_i.sanity_caught_err = 1;
			}
		}
		Explainer_parameters ex_params;
		ma.get_explainer_params(ex_params);
		// MLOG("=====> Working on i %d pid %d time %d sanity_test_rc %d sanity_caught_err %d\n", i, req_i.sample_pid, req_i.sample_time, req_i.sanity_test_rc, req_i.sanity_caught_err);
		nlohmann::ordered_json js = nlohmann::ordered_json({});

		js.push_back({"patient_id", to_string(req_i.sample_pid)});
		js.push_back({"time", to_string(req_i.sample_time)});
		// Store error message from loading if has
		if (i < all_load_data_msgs.size())
			for (size_t j = 0; j < all_load_data_msgs[i].size(); ++j)
				add_to_json_array(js, "messages", all_load_data_msgs[i][j]);

		if (req_i.sanity_caught_err)
			add_to_json_array(js, "messages", "ERROR: sanity tests crashed");
		if (req_i.sanity_res.size() > 0)
			for (auto &its : req_i.sanity_res)
			{
				add_to_json_array(js, "messages", "(" + to_string(its.external_rc) + ")" + its.err_msg);
			}
		for (auto &e : req_i.exports)
		{
			if (e.second.type == PREDICTION_SOURCE_PREDICTIONS)
			{
				if (req_i.res != NULL && req_i.res->prediction.size() > e.second.pred_channel && req_i.sanity_caught_err == 0 && req_i.sanity_test_rc > 0)
				{
					js.push_back({e.first, to_string(req_i.res->prediction[e.second.pred_channel])});
					// Add Flag if configured:
					add_flag_response(js, req_i.res->prediction[e.second.pred_channel], ma, req_i.flag_threshold, req_i.flag_threshold_numeric);
				}
				else
					js.push_back({e.first, to_string(AM_UNDEFINED_VALUE)});
				if (req_i.res == NULL)
					add_to_json_array(js, "messages", "ERROR: did not get result for field " + e.first + " : " + e.second.field);
				else if (req_i.res->prediction.size() <= e.second.pred_channel || e.second.pred_channel < 0)
					add_to_json_array(js, "messages", "ERROR: prediction channel " + to_string(e.second.pred_channel) + " is illegal");
			}
			else if (e.second.type == PREDICTION_SOURCE_ATTRIBUTE && req_i.res != NULL)
			{
				if (req_i.res->attributes.find(e.second.field) != req_i.res->attributes.end())
					js.push_back({e.first, to_string(req_i.res->attributes[e.second.field])});
				else if (req_i.res->str_attributes.find(e.second.field) != req_i.res->str_attributes.end())
					js.push_back({e.first, req_i.res->str_attributes[e.second.field]});
			}
			else if (e.second.type == PREDICTION_SOURCE_ATTRIBUTE_AS_JSON && req_i.res != NULL && req_i.sanity_caught_err == 0 && req_i.sanity_test_rc > 0)
			{
				if (req_i.res->str_attributes.find(e.second.field) != req_i.res->str_attributes.end())
				{
					nlohmann::ordered_json jattr;
					try
					{
						jattr = nlohmann::ordered_json::parse(req_i.res->str_attributes[e.second.field]);
						process_explainability(jattr, ex_params, *rep, req_i.sample_pid, req_i.sample_time);
						js.push_back({e.first, jattr});
					}
					catch (...)
					{
						add_to_json_array(jresp, "messages", "ERROR: could not parse attribute " + e.second.field + " as a valid json");
					}
				}
			}
			else if (e.second.type == PREDICTION_SOURCE_JSON && req_i.res != NULL && req_i.sanity_caught_err == 0 && req_i.sanity_test_rc > 0)
			{
				if (req_i.res->jrec.find(e.second.field) != req_i.res->jrec.end())
				{
					auto &jj = req_i.res->jrec[e.second.field];
					js.push_back({e.first, jj});
				}
			}
		}
		jresp["responses"].push_back(js);
	}

	json_to_char_ptr(jresp, response);

#ifdef AM_TIMING_LOGS
	timer.take_curr_time();
	MLOG("INFO:: MedialInfraAlgoMarker::CalculateByType :: end create_response %2.1f milisecond\n", timer.diff_milisec());
	timer.start();
#endif

	if (am_matrix != "" && _n_points > 0)
	{ // debug only
		if (first_write)
			ma.write_features_mat(am_matrix);
		else
			ma.add_features_mat(am_matrix);
		first_write = false;
	}

	return AM_OK_RC;
}

//-----------------------------------------------------------------------------------
int MedialInfraAlgoMarker::AdditionalLoad(const int LoadType, const char *load)
{
	if (!is_loaded)
		return AM_ERROR_MUST_BE_LOADED;
	if (LoadType != LOAD_DICT_FROM_FILE && LoadType != LOAD_DICT_FROM_JSON)
		return AM_ERROR_UNKNOWN_LOAD_TYPE;

	json js;

	if (LoadType == LOAD_DICT_FROM_FILE)
	{
		string sload;
		string f_in(load);
		if (read_file_into_string(f_in, sload) < 0)
			return AM_ERROR_READING_DICT_FILE;
		js = json::parse(sload.c_str());
	}
	else
		js = json::parse(load);

	try
	{
		ma.add_json_dict(js);
	}
	catch (...)
	{
		return AM_ERROR_PARSING_JSON_DICT;
	}

	// now that we added the json dictionary, we need to reinitialize the model ! as it needs to prepare potential tables using these new definitions and sets
	// ma.init_model_for_apply();

	return AM_OK_RC;
}

//-----------------------------------------------------------------------------------
// private internals for class MedialInfraAlgoMarker
//-----------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------
int MedialInfraAlgoMarker::read_config(const string &conf_f)
{
	set_config(conf_f.c_str());

	ifstream inf(conf_f);

	if (!inf)
		return AM_ERROR_LOAD_NO_CONFIG_FILE;

	string dir = conf_f.substr(0, conf_f.find_last_of("/\\"));
	string curr_line;
	while (getline(inf, curr_line))
	{
		if ((curr_line.size() > 1) && (curr_line[0] != '#'))
		{

			if (curr_line[curr_line.size() - 1] == '\r')
				curr_line.erase(curr_line.size() - 1);

			vector<string> fields;
			split(fields, curr_line, boost::is_any_of("\t"));

			if (fields.size() >= 2)
			{
				if (fields[0] == "TYPE")
					type_in_config_file = fields[1];
				else if (fields[0] == "REPOSITORY")
					rep_fname = fields[1];
				else if (fields[0] == "MODEL")
					model_fname = fields[1];
				else if (fields[0] == "ALLOW_REP_ADJUSTMENTS")
					allow_rep_adjustments = stoi(fields[1]) > 0;
				else if (fields[0] == "MODEL_END_STAGE")
					try
					{
						model_end_stage = stoi(fields[1]);
					}
					catch (...)
					{
						MTHROW_AND_ERR("Could not parse given value MODEL_END_STAGE='%s'\n", fields[1].c_str());
					}
				else if (fields[0] == "EXTENDED_RESULT_FIELDS")
					split(extended_result_fields, fields[1], boost::is_any_of(";"));
				else if (fields[0] == "INPUT_TESTER_CONFIG")
					input_tester_config_file = fields[1];
				else if (fields[0] == "NAME")
					set_name(fields[1].c_str());
				else if (fields[0] == "TIME_UNIT")
				{
					set_time_unit(med_time_converter.string_to_type(fields[1].c_str()));
				}
				else if (fields[0] == "DEBUG_MATRIX")
					am_matrix = fields[1];
				else if (fields[0] == "AM_UDI_DI")
					set_am_udi_di(fields[1].c_str());
				else if (fields[0] == "AM_MANUFACTOR_DATE")
					set_manafactur_date(fields[1].c_str());
				else if (fields[0] == "AM_VERSION")
					set_am_version(fields[1].c_str());
				else if (fields[0] == "EXPLAINABILITY_PARAMS")
					ma.set_explainer_params(fields[1], dir);
				else if (fields[0] == "THRESHOLD_LEAFLET")
					ma.set_threshold_leaflet(fields[1], dir);
				else if (fields[0] == "TESTER_NAME")
				{
				}
				else if (fields[0] == "FILTER")
				{
				}
				else
					MWARN("WRAN: unknown parameter \"%s\". Read and ignored\n", fields[0].c_str());
			}
		}
	}

	if (rep_fname != "" && rep_fname[0] != '/' && rep_fname[0] != '\\')
	{
		// relative path
		rep_fname = dir + "/" + rep_fname;
	}

	if (model_fname != "" && model_fname[0] != '/' && model_fname[0] != '\\')
	{
		// relative path
		model_fname = dir + "/" + model_fname;
	}

	if (input_tester_config_file == ".")
	{
		input_tester_config_file = conf_f; // option to use the general config file as the file to config the tester as well.
	}
	else if (input_tester_config_file != "" && input_tester_config_file[0] != '/' && input_tester_config_file[0] != '\\')
	{
		// relative path
		input_tester_config_file = dir + "/" + input_tester_config_file;
	}

	return AM_OK_RC;
}

// maximal input of 32GB
#define MAX_POSSIBLE_STRING_LEN ((size_t)1 << 35)
//-----------------------------------------------------------------------------------
void MedialInfraAlgoMarker::get_jsons_locations(const char *data, vector<size_t> &j_start, vector<size_t> &j_len)
{
	j_start.clear();
	j_len.clear();

	size_t j = 0;
	int counter = 0;
	int in_string = 0;
	size_t start = 0;
	size_t len = 0;
	while (data[j] != 0 || j > MAX_POSSIBLE_STRING_LEN)
	{
		char ch = data[j];
		if (ch == '\"' || ch == '\'')
			in_string = 1 - in_string;
		if ((!in_string) && ch == '{')
		{
			if (counter == 0)
				start = j;
			counter++;
		}
		if (counter > 0)
			len++;
		if ((!in_string) && ch == '}')
			counter--;

		if (counter == 0 && len > 0)
		{
			j_start.push_back(start);
			j_len.push_back(len);
			len = 0;
			if (j_start.size() > 0 && j_start.size() % 1000 == 0)
				MLOG("Found %d jsons so far\n", j_start.size());
		}
		if (counter < 0)
			MTHROW_AND_ERR("Mismatch in {} count in jsons string\n");

		j++;
	}

	// MLOG("Read %d jsons from data string (debug info: counter = %d j = %ld)\n", j_start.size(), counter, j);
}

string construct_message(int code, const string &msg)
{
	return "(" + to_string(code) + ")" + msg;
}

//-----------------------------------------------------------------------------------
int MedialInfraAlgoMarker::AddJsonData(int patient_id, json &j_data, vector<string> &messages, map<pair<int, int>, pair<int, vector<char>>> *data)
{
	string current_time = "";

	MedRepository &rep = ma.get_rep();
	if (data == NULL) // default pointer to global rep
		data = &rep.in_mem_rep.data;

	bool good = true;
	bool mark_succ_ = false;
	try
	{
		json &js = j_data;

		// supporting also older style jsons that were embeded in a "body" section
		if (j_data.find("body") != j_data.end())
			js = j_data["body"];

		if (patient_id <= 0)
		{
			// in this case we take the patient id directly from the json itself
			if (js.find("patient_id") != js.end())
			{
				if (js["patient_id"].is_number_integer())
					patient_id = js["patient_id"].get<long long>();
				else if (js["patient_id"].is_string())
				{
					try
					{
						patient_id = stoll(js["patient_id"].get<string>());
					}
					catch (...)
					{
						messages.push_back(construct_message(AM_DATA_BAD_FORMAT_FATAL, "Bad data json format - couldn't convert patient_id to integer"));
						return AM_FAIL_RC;
					}
				}
				else
				{
					messages.push_back(construct_message(AM_DATA_BAD_FORMAT_FATAL, "Bad data json format - patient_id suppose to be integer"));
					return AM_FAIL_RC;
				}
			}
			else
			{
				if (js.find("pid") != js.end())
				{
					if (js["pid"].is_number_integer())
						patient_id = js["pid"].get<long long>();
					else if (js["pid"].is_string())
					{
						try
						{
							patient_id = stoll(js["pid"].get<string>());
						}
						catch (...)
						{
							messages.push_back(construct_message(AM_DATA_BAD_FORMAT_FATAL, "Bad data json format - couldn't convert pid to integer"));
							return AM_FAIL_RC;
						}
					}
					else
					{
						messages.push_back(construct_message(AM_DATA_BAD_FORMAT_FATAL, "Bad data json format - pid suppose to be integer"));
						return AM_FAIL_RC;
					}
				}
			}
		}
		if (patient_id <= 0)
		{
			messages.push_back(construct_message(AM_DATA_BAD_FORMAT_FATAL, "Bad data json format - no patient_id was given"));
			return AM_FAIL_RC;
		}

		// MLOG("Loading with pid %d\n", patient_id);

		vector<long long> times;
		int s_data_size = 100000;
		vector<char> sdata(s_data_size);
		vector<int> sinds;
		int curr_s = 0;

		// char str_values[MAX_VALS][MAX_VAL_LEN];
		if (js.find("signals") == js.end() || !js["signals"].is_array())
		{
			char buf[5000];
			if (patient_id != 1)
				snprintf(buf, sizeof(buf), "(%d)Bad format in patient %d. Element should contain signals element as array",
						 AM_DATA_BAD_FORMAT_FATAL, patient_id);
			else
				snprintf(buf, sizeof(buf), "(%d)Bad format. Element should contain signals element as array",
						 AM_DATA_BAD_FORMAT_FATAL);
			messages.push_back(string(buf));
			get_current_time(current_time);
			MLOG("%s::%s\n", current_time.c_str(), buf);
			good = false;
		}
		else
		{
			for (auto &s : js["signals"])
			{
				bool good_sig = true;
				int n_time_channels, n_val_channels, *is_categ;
				string sig;
				times.clear();
				sinds.clear();
				curr_s = 0;
				if (s.find("code") == s.end() || !s["code"].is_string())
				{
					char buf[5000];
					if (patient_id != 1)
						snprintf(buf, sizeof(buf), "(%d)Bad format in patient %d. Element should contain code element as signal name",
								 AM_DATA_BAD_FORMAT_FATAL, patient_id);
					else
						snprintf(buf, sizeof(buf), "(%d)Bad format. Element should contain code element as signal name",
								 AM_DATA_BAD_FORMAT_FATAL);
					messages.push_back(string(buf));
					get_current_time(current_time);
					MLOG("%s::%s\n", current_time.c_str(), buf);
					good = false;
					good_sig = false;
				}
				if (good_sig)
				{
					sig = s["code"].get<string>();
					int sid = rep.sigs.Name2Sid[sig];
					get_sig_structure(sig, n_time_channels, n_val_channels, is_categ);
					if (n_time_channels == 0 && n_val_channels == 0)
					{
						char buf[5000];
						if (patient_id != 1)
							snprintf(buf, sizeof(buf), "(%d)An unknown signal was found: %s for patient %d",
									 AM_DATA_UNKNOWN_SIGNAL, sig.c_str(), patient_id);
						else
							snprintf(buf, sizeof(buf), "(%d)An unknown signal was found: %s",
									 AM_DATA_UNKNOWN_SIGNAL, sig.c_str());

						messages.push_back(string(buf));
						get_current_time(current_time);
						MLOG("%s::%s\n", current_time.c_str(), buf);
						good = false;
						good_sig = false;
						// return AM_FAIL_RC;
						continue;
					}
					// MLOG("%s %d %d\n", sig.c_str(), n_time_channels, n_val_channels);
					int n_data = 0;
					if (s.find("data") == s.end() || !s["data"].is_array())
					{
						char buf[5000];
						if (patient_id != 1)
							snprintf(buf, sizeof(buf), "(%d)Bad format for signal: %s in patient %d. No data element or data element is not array",
									 AM_DATA_BAD_FORMAT_FATAL, sig.c_str(), patient_id);
						else
							snprintf(buf, sizeof(buf), "(%d)Bad format for signal: %s. No data element or data element is not array",
									 AM_DATA_BAD_FORMAT_FATAL, sig.c_str());
						messages.push_back(string(buf));
						get_current_time(current_time);
						MLOG("%s::%s\n", current_time.c_str(), buf);
						good = false;
						good_sig = false;
					}
					if (good_sig)
					{
						for (auto &d : s["data"])
						{
							int nt = 0;
							bool good_record = true;
							if (d.find("timestamp") != d.end() && !d["timestamp"].is_array())
							{
								char buf[5000];
								if (patient_id != 1)
									snprintf(buf, sizeof(buf), "(%d)Bad format for signal: %s in patient %d. timestamp should be array of timestamps, each represents a different channel.",
											 AM_DATA_BAD_FORMAT_NON_FATAL, sig.c_str(), patient_id);
								else
									snprintf(buf, sizeof(buf), "(%d)Bad format for signal: %s. timestamp should be array of timestamps, each represents a different channel.",
											 AM_DATA_BAD_FORMAT_NON_FATAL, sig.c_str());
								messages.push_back(string(buf));
								get_current_time(current_time);
								MLOG("%s::%s\n", current_time.c_str(), buf);
								good = false;
								break;
							}
							if (d.find("value") != d.end() && !d["value"].is_array())
							{
								char buf[5000];
								if (patient_id != 1)
									snprintf(buf, sizeof(buf), "(%d)Bad format for signal: %s in patient %d. value should be array of values, each represents a different channel.",
											 AM_DATA_BAD_FORMAT_NON_FATAL, sig.c_str(), patient_id);
								else
									snprintf(buf, sizeof(buf), "(%d)Bad format for signal: %s. value should be array of values, each represents a different channel.",
											 AM_DATA_BAD_FORMAT_NON_FATAL, sig.c_str());
								messages.push_back(string(buf));
								get_current_time(current_time);
								MLOG("%s::%s\n", current_time.c_str(), buf);
								good = false;
								break;
							}
							for (auto &t : d["timestamp"])
							{
								char buf[5000];
								if (t.is_string())
								{
									try
									{
										times.push_back(stoll(t.get<string>()));
										++nt;
										continue;
									}
									catch (...)
									{
										if (patient_id != 1)
											snprintf(buf, sizeof(buf), "(%d)Bad format for signal: %s in patient %d. Couldn't convert timestamp to integer",
													 AM_DATA_BAD_FORMAT_NON_FATAL, sig.c_str(), patient_id);
										else
											snprintf(buf, sizeof(buf), "(%d)Bad format for signal: %s. Couldn't convert timestamp to integer",
													 AM_DATA_BAD_FORMAT_NON_FATAL, sig.c_str());
										messages.push_back(string(buf));
										get_current_time(current_time);
										MLOG("%s::%s\n", current_time.c_str(), buf);
										good = false;
										good_sig = false;
										good_record = false;
										break;
									}
								}
								else if (!t.is_number_integer())
								{

									if (patient_id != 1)
										snprintf(buf, sizeof(buf), "(%d)Bad format for signal: %s in patient %d. timestamp element should be integer.",
												 AM_DATA_BAD_FORMAT_NON_FATAL, sig.c_str(), patient_id);
									else
										snprintf(buf, sizeof(buf), "(%d)Bad format for signal: %s. timestamp element should be integer.",
												 AM_DATA_BAD_FORMAT_NON_FATAL, sig.c_str());
									messages.push_back(string(buf));
									get_current_time(current_time);
									MLOG("%s::%s\n", current_time.c_str(), buf);
									good = false;
									good_sig = false;
									good_record = false;
									break;
								}

								times.push_back(t.get<long long>());
								++nt;
							}
							if (!good_record)
								break;
							// Check size of timestamps:
							if (nt != n_time_channels)
							{
								char buf[5000];
								if (patient_id != 1)
									snprintf(buf, sizeof(buf), "(%d)Bad signal structure for signal: %s in patient %d. Received %d time channels instead of %d",
											 AM_DATA_BAD_STRUCTURE, sig.c_str(), patient_id, nt, n_time_channels);
								else
									snprintf(buf, sizeof(buf), "(%d)Bad signal structure for signal: %s. Received %d time channels instead of %d",
											 AM_DATA_BAD_STRUCTURE, sig.c_str(), nt, n_time_channels);
								messages.push_back(string(buf));
								get_current_time(current_time);
								MLOG("%s::%s\n", current_time.c_str(), buf);
								good = false;
								good_sig = false;
								good_record = false;
								// return AM_FAIL_RC;
							}
							if (!good_record)
								break;
							int nv = 0;
							for (auto &v : d["value"])
							{
								string sv;
								if (v.is_number() && !is_categ[nv])
								{
									sv = to_string(v.get<double>());
								}
								else
								{
									if (!v.is_string())
									{
										char buf[5000];
										if (patient_id != 1)
											snprintf(buf, sizeof(buf), "(%d)Bad format for signal: %s in patient %d. value element should be string.",
													 AM_DATA_BAD_FORMAT_NON_FATAL, sig.c_str(), patient_id);
										else
											snprintf(buf, sizeof(buf), "(%d)Bad format for signal: %s. value element should be string.",
													 AM_DATA_BAD_FORMAT_NON_FATAL, sig.c_str());
										messages.push_back(string(buf));
										get_current_time(current_time);
										MLOG("%s::%s\n", current_time.c_str(), buf);
										good = false;
										good_sig = false;
										good_record = false;
										break;
									}
									else
										sv = v.get<string>().c_str();
								}

								// Check if "Date"
								string unit_m = rep.sigs.unit_of_measurement(sid, nv);
								boost::to_lower(unit_m);
								if (unit_m == "date")
								{
									try
									{
										int full_date = (int)stod(sv);
										// check if valid date?
										if (!med_time_converter.is_valid_date(full_date))
										{
											char buf[5000];
											if (patient_id != 1)
												snprintf(buf, sizeof(buf), "(%d)Bad format for signal: %s in patient %d. value should be date format. Recieved %d.",
														 AM_DATA_BAD_FORMAT_NON_FATAL, sig.c_str(), patient_id, full_date);
											else
												snprintf(buf, sizeof(buf), "(%d)Bad format for signal: %s. value should be date format. Recieved %d.",
														 AM_DATA_BAD_FORMAT_NON_FATAL, sig.c_str(), full_date);
											messages.push_back(string(buf));
											get_current_time(current_time);
											MLOG("%s::%s\n", current_time.c_str(), buf);
											good = false;
											good_sig = false;
											good_record = false;
											break;
										}
									}
									catch (...)
									{
										char buf[5000];
										if (patient_id != 1)
											snprintf(buf, sizeof(buf), "(%d)Bad format for signal: %s in patient %d. value should be date format. Recieved %s.",
													 AM_DATA_BAD_FORMAT_NON_FATAL, sig.c_str(), patient_id, sv.c_str());
										else
											snprintf(buf, sizeof(buf), "(%d)Bad format for signal: %s. value should be date format. Recieved %s.",
													 AM_DATA_BAD_FORMAT_NON_FATAL, sig.c_str(), sv.c_str());
										messages.push_back(string(buf));
										get_current_time(current_time);
										MLOG("%s::%s\n", current_time.c_str(), buf);
										good = false;
										good_sig = false;
										good_record = false;
										break;
									}
								}

								int slen = (int)sv.length();
								// MLOG("val %d : %s len: %d curr_s %d s_data_size %d %d n_val_channels %d\n", nv, sv.c_str(), slen, curr_s, s_data_size, sdata.size(), n_val_channels);
								if (curr_s + 1 + slen > s_data_size)
								{
									s_data_size *= 2;
									sdata.resize(s_data_size);
								}
								sv.copy(&sdata[curr_s], slen);
								sdata[curr_s + slen] = 0;
								sinds.push_back(curr_s);
								curr_s += slen + 1;
								++nv;
								// char *sp = &sdata[sinds.back()];
								// MLOG("val %d %d %s : %s len: %d curr_s %d s_data_size %d %d\n", sinds.size(), sinds.back(), sp, sv.c_str(), slen, curr_s, s_data_size, sdata.size());
								// MLOG("%s ", v.get<string>().c_str());
							}
							if (!good_record)
								break;
							// Check size of value:
							if (nv != n_val_channels)
							{
								char buf[5000];
								if (patient_id != 1)
									snprintf(buf, sizeof(buf), "(%d)Bad signal structure for signal: %s in patient %d. Received %d value instead of %d",
											 AM_DATA_BAD_STRUCTURE, sig.c_str(), patient_id, nv, n_val_channels);
								else
									snprintf(buf, sizeof(buf), "(%d)Bad signal structure for signal: %s. Received %d value channels instead of %d",
											 AM_DATA_BAD_STRUCTURE, sig.c_str(), nv, n_val_channels);
								messages.push_back(string(buf));
								get_current_time(current_time);
								MLOG("%s::%s\n", current_time.c_str(), buf);
								good = false;
								good_sig = false;
								good_record = false;
								// return AM_FAIL_RC;
							}
							// MLOG("\n");
							if (!good_record)
								break;
							n_data++;
						}
					}
				}
				vector<char *> p_str;
				for (auto j : sinds)
					p_str.push_back(&sdata[j]);
				long long *p_times = &times[0];
				int n_times = (int)times.size();
				char **str_values = &p_str[0];
				int n_vals = (int)p_str.size();

				// MLOG("%s n_times %d n_vals %d n_data %d\n", sig.c_str(), n_times, n_vals, n_data);
				// MLOG("times: "); for (int j = 0; j < n_times; j++) MLOG("%d,", p_times[j]); 	MLOG("\nvals: ");
				// for (int j = 0; j < n_vals; j++) MLOG("%s, ", str_values[j]); MLOG("\n");

				if (good_sig)
				{
					if (AddDataStr_data(patient_id, sig.c_str(), n_times, p_times, n_vals, str_values, data) != AM_OK_RC)
					{
						char buf[5000];
						if (patient_id != 1)
							snprintf(buf, sizeof(buf), "(%d)General error in signal: %s for patient %d",
									 AM_DATA_GENERAL_ERROR, sig.c_str(), patient_id);
						else
							snprintf(buf, sizeof(buf), "(%d)General error in signal: %s",
									 AM_DATA_GENERAL_ERROR, sig.c_str());
						messages.push_back(string(buf));
						get_current_time(current_time);
						MLOG("%s::%s\n", current_time.c_str(), buf);
						good_sig = false;
						good = false;
						// return AM_FAIL_RC;
					}
					else
						mark_succ_ = true;
				}
			}
		}
	}
	catch (...)
	{
		char buf[5000];
		snprintf(buf, sizeof(buf), "(%d)Bad data json format",
				 AM_DATA_BAD_FORMAT_FATAL);
		messages.push_back(string(buf));
		get_current_time(current_time);
		MLOG("%s::%s\n", current_time.c_str(), buf);
		good = false;
		// return AM_FAIL_RC;
	}
	if (!good)
	{
		if (mark_succ_) // add message that some was loaded:
			MLOG("AddDataByType() WARN - some of the data signals were loaded for patient %d. Consider calling ClearData if rerun again after fixing.\n",
				 patient_id);
		return AM_FAIL_RC;
	}

	return AM_OK_RC;
}

//-----------------------------------------------------------------------------------
string MedialInfraAlgoMarker::get_lib_code_version()
{
	return medial::get_git_version();
}

//-----------------------------------------------------------------------------------
int MedialInfraAlgoMarker::Discovery(char **response)
{
	nlohmann::ordered_json jresp;
	jresp = nlohmann::ordered_json(
		{{"name", get_name()},
		 {"version", get_am_version()},
		 {"udi", get_am_udi_di()},
		 {"mfg_date", get_manfactor_date()},
		 {"model_code_version", ma.model_version_info()},
		 {"lib_code_version", get_lib_code_version()},
		 {"signals", nlohmann::ordered_json::array()}});

	nlohmann::ordered_json &json_signals = jresp["signals"];

	// Add signals to json_signals
	if (!ma.model_initiated())
		ma.init_model_for_rep();

	unordered_set<string> sigs;
	get_am_rep_signals(sigs);
	vector<string> req_sigs;
	unordered_map<string, vector<string>> sig_categ, sig_categ_final;
	ma.get_model_signals_info(req_sigs, sig_categ);
	unordered_set<string> req_set(req_sigs.begin(), req_sigs.end());
	// rename sig_categ if needed in some cases!
	for (const auto &it : sig_categ)
	{
		if (req_set.find(it.first) != req_set.end())
			sig_categ_final[it.first] = it.second;
		else
		{
			// test by section name!
			int sig_section = ma.get_rep().dict.section_id(it.first);
			// retrieve name that exists in SECTION:
			const unordered_set<string> &sec_names = ma.get_rep().dict.dict(sig_section)->section_name;
			vector<string> candidates;
			for (const string &cand : sec_names)
			{
				if (req_set.find(cand) == req_set.end())
					continue;
				candidates.push_back(cand);
			}

			// Add to all:
			if (candidates.empty())
				MWARN("Warn - has used categorical signal %s without mapping\n",
					  it.first.c_str());
			else if (candidates.size() > 1)
				MWARN("Warn - has used categorical signal %s with multiple mapping\n",
					  it.first.c_str());
			else
			{
				unordered_set<string> sig_list_c(sig_categ_final[candidates[0]].begin(), sig_categ_final[candidates[0]].end());
				sig_list_c.insert(it.second.begin(), it.second.end());
				vector<string> final_list(sig_list_c.begin(), sig_list_c.end());
				sort(final_list.begin(), final_list.end());
				sig_categ_final[candidates[0]] = move(final_list);
			}
		}
	}

	for (const string &sig_name : req_sigs)
	{
		string sig_nm = sig_name;
		string sig_unit = "";
		string sig_type = "";
		vector<int> categorical_ch;
		vector<string> categorical_vals;
		if (sigs.find(sig_name) != sigs.end())
		{
			int sid = ma.get_rep().sigs.sid(sig_name);
			const SignalInfo &si = ma.get_rep().sigs.Sid2Info[sid];
			for (int i = 0; i < si.n_val_channels; ++i)
			{
				if (i > 0)
					sig_unit += ",";
				sig_unit += si.unit_of_measurement_per_val_channel[i];
			}
			UniversalSigVec usv;
			usv.init(si);
			sig_type = usv.get_signal_generic_spec();
			categorical_ch.resize(si.n_val_channels);
			for (size_t i = 0; i < si.n_val_channels; ++i)
				if (si.is_categorical_per_val_channel[i])
					categorical_ch[i] = 1;
		}
		else
		{
			sig_nm = sig_nm + "(virtual)";
		}
		if (sig_categ_final.find(sig_name) != sig_categ_final.end())
			categorical_vals = std::move(sig_categ_final.at(sig_name));

		nlohmann::ordered_json sig_js;
		sig_js = {
			{"code", sig_nm},
			{"unit", sig_unit},
			{"type", sig_type},
			{"categorical_channels", categorical_ch},
			{"categorical_values", categorical_vals}};

		json_signals += sig_js;
	}

	Explainer_parameters ex_params;
	ma.get_explainer_params(ex_params);
	if (ex_params.max_threshold > 0)
	{ // has settings
		MedPidRepository &rep = ma.get_rep();
		jresp["explainability_options"] = nlohmann::ordered_json::array();
		vector<string> model_groups; // fetch from model
		ma.get_explainer_output_options(model_groups);
		// filer ignore:
		for (const string &str_g : model_groups)
			if (ex_params.ignore_groups_list.find(str_g) == ex_params.ignore_groups_list.end())
			{
				nlohmann::ordered_json element_exp;
				element_exp["name"] = str_g;
				element_exp["description"] = fetch_desription(rep, str_g);
				jresp["explainability_options"].push_back(element_exp);
			}
	}
	// ma.get_rep().sigs.Sid2Info[1].

	vector<string> mbr_opts;
	ma.fetch_all_thresholds(mbr_opts);
	if (!mbr_opts.empty())
	{
		jresp["default_threshold"] = ma.get_default_threshold();
		jresp["flag_threshold_options"] = json::array();
		for (const string &opt : mbr_opts)
			jresp["flag_threshold_options"].push_back(opt);
	}

	json_to_char_ptr(jresp, response);
	return 0;
}

//-----------------------------------------------------------------------------------
void add_to_json_array(json &js, const string &key, const string &s_add)
{
	if (js.find(key) == js.end())
		js.push_back({key, json::array()});
	js[key] += s_add;
}
void add_to_json_array(nlohmann::ordered_json &js, const string &key, const string &s_add)
{
	if (js.find(key) == js.end())
		js.push_back({key, json::array()});
	js[key] += s_add;
}

//-----------------------------------------------------------------------------------
void json_to_char_ptr(json &js, char **jarr)
{
	*jarr = NULL;
	string sj = js.dump(1, '\t');

	*jarr = new char[sj.length() + 1];

	if (*jarr != NULL)
	{
		(*jarr)[sj.length()] = 0;
		strncpy(*jarr, sj.c_str(), sj.length());
	}
}

void json_to_char_ptr(nlohmann::ordered_json &js, char **jarr)
{
	*jarr = NULL;
	string sj = js.dump(1, '\t');

	*jarr = new char[sj.length() + 1];

	if (*jarr != NULL)
	{
		(*jarr)[sj.length()] = 0;
		strncpy(*jarr, sj.c_str(), sj.length());
	}
}

//-----------------------------------------------------------------------------------
bool json_verify_key(json &js, const string &key, int verify_val_flag, const string &val)
{
	bool is_in = false;
	if (js.find(key) != js.end())
		is_in = true;

	if (is_in && verify_val_flag)
	{
		if (!js[key].is_string() || (js[key].get<string>() != val))
			is_in = false;
	}

	return is_in;
}
bool json_verify_key(nlohmann::ordered_json &js, const string &key, int verify_val_flag, const string &val)
{
	bool is_in = false;
	if (js.find(key) != js.end())
		is_in = true;

	if (is_in && verify_val_flag)
	{
		if (!js[key].is_string() || (js[key].get<string>() != val))
			is_in = false;
	}

	return is_in;
}

//------------------------------------------------------------------------------------------
int json_parse_request(json &jreq, json_req_info &defaults, json_req_info &req_i, string &error_message)
{
	req_i = defaults;
	error_message = "";
	// read defaults (if exist)
	try
	{
		if (json_verify_key(jreq, "patient_id", 0, "") || json_verify_key(jreq, "pid", 0, ""))
		{
			if (json_verify_key(jreq, "patient_id", 0, ""))
			{
				if (jreq["patient_id"].is_string())
					req_i.sample_pid = stoi(jreq["patient_id"].get<string>());
				else if (jreq["patient_id"].is_number_integer())
					req_i.sample_pid = jreq["patient_id"].get<int>();
				else
				{
					error_message = "Error in patient_id field - unsupported type";
					MTHROW_AND_ERR("Error in patient_id field - unsupported type\n");
				}
			}
			else
			{
				if (jreq["pid"].is_string())
					req_i.sample_pid = stoi(jreq["pid"].get<string>());
				else if (jreq["pid"].is_number_integer())
					req_i.sample_pid = jreq["pid"].get<int>();
				else
				{
					error_message = "Error in patient_id field - unsupported type";
					MTHROW_AND_ERR("Error in pid field - unsupported type\n");
				}
			}
		}

		if (json_verify_key(jreq, "scoreOnDate", 0, "") || json_verify_key(jreq, "time", 0, ""))
		{
			if (json_verify_key(jreq, "scoreOnDate", 0, ""))
			{
				if (jreq["scoreOnDate"].is_string())
					req_i.sample_time = stoll(jreq["scoreOnDate"].get<string>());
				else if (jreq["scoreOnDate"].is_number_integer())
					req_i.sample_time = jreq["scoreOnDate"].get<long long>();
				else
				{
					error_message = "Error in scoreOnDate field - unsupported type";
					MTHROW_AND_ERR("Error in scoreOnDate field - unsupported type\n");
				}
			}
			else
			{
				if (jreq["time"].is_string())
					req_i.sample_time = stoll(jreq["time"].get<string>());
				else if (jreq["time"].is_number_integer())
					req_i.sample_time = jreq["time"].get<long long>();
				else
				{
					error_message = "Error in time field - unsupported type";
					MTHROW_AND_ERR("Error in time field - unsupported type\n");
				}
			}
		}

		if (json_verify_key(jreq, "load", 0, ""))
		{
			if (jreq["load"].is_string())
				req_i.load_data = stoi(jreq["load"].get<string>());
			else if (jreq["load"].is_number_integer())
				req_i.load_data = jreq["load"].get<int>();
			else
			{
				error_message = "Error in load field - unsupported type";
				MTHROW_AND_ERR("Error in load field - unsupported type\n");
			}
		}

		if (json_verify_key(jreq, "export", 0, ""))
		{

			for (auto &jexp : jreq["export"].items())
			{

				string name = jexp.key();
				string field = jexp.value().get<string>();

				// MLOG("Working on %s : %s\n", name.c_str(), field.c_str());
				int type = PREDICTION_SOURCE_UNKNOWN;
				int pred_channel = -1;

				vector<string> f;
				boost::split(f, field, boost::is_any_of(" "));
				if (f.size() == 2)
				{
					if (f[0] == "attr")
					{
						type = PREDICTION_SOURCE_ATTRIBUTE;
						field = f[1];
					}
					else if (f[0] == "json_attr")
					{
						type = PREDICTION_SOURCE_ATTRIBUTE_AS_JSON;
						field = f[1];
					}
					else if (f[0] == "pred")
					{
						type = PREDICTION_SOURCE_PREDICTIONS;
						field = "pred_" + f[1];
					}
					else if (f[0] == "json")
					{
						type = PREDICTION_SOURCE_JSON;
						field = f[1];
					}
				}

				if ((type == PREDICTION_SOURCE_UNKNOWN || type == PREDICTION_SOURCE_PREDICTIONS) && (field.length() > 5) && (field.substr(0, 5) == "pred_"))
				{
					type = PREDICTION_SOURCE_PREDICTIONS;
					pred_channel = stoi(field.substr(5));
				}

				if (type == PREDICTION_SOURCE_UNKNOWN)
					type = PREDICTION_SOURCE_ATTRIBUTE;

				json_req_export jexport;

				jexport.field = field;
				jexport.pred_channel = pred_channel;
				jexport.type = type;

				req_i.exports[name] = jexport;
			}
		}

		if (json_verify_key(jreq, "flag_threshold", 0, ""))
		{
			if (!jreq["flag_threshold"].is_string())
			{
				error_message = "Error in flag_threshold field - unsupported type, expecting string";
				MTHROW_AND_ERR("Error in flag_threshold field - unsupported type, expecting string\n");
			}
			req_i.flag_threshold = jreq["flag_threshold"].get<string>();
		}

		if (json_verify_key(jreq, "flag_threshold_numeric", 0, ""))
		{
			if (jreq["flag_threshold_numeric"].is_number())
				req_i.flag_threshold_numeric = jreq["flag_threshold_numeric"].get<float>();
			else if (jreq["flag_threshold_numeric"].is_string())
			{
				try
				{
					req_i.flag_threshold_numeric = stof(jreq["flag_threshold_numeric"].get<string>());
				}
				catch (...)
				{
					error_message = "Error in flag_threshold_numeric field - unsupported type, expecting float";
					return -1;
				}
			}
			else
			{
				error_message = "Error in flag_threshold_numeric field - unsupported type, expecting float";
				MTHROW_AND_ERR("Error in flag_threshold_numeric field - unsupported type, expecting float\n");
			}
		}
	}
	catch (...)
	{

		return -1;
	}
	return 0;
}

void Explainer_description_config::read_cfg_file(const string &file)
{
	ifstream file_reader(file);
	if (!file_reader.good())
		MTHROW_AND_ERR("Error Explainer_description_config::read_cfg_file - file %s wasn't found\n",
					   file.c_str());
	records.clear();
	unordered_set<string> contrib_seen;

	string line;
	while (getline(file_reader, line))
	{
		mes_trim(line);
		if (line.empty() || line[0] == '#')
			continue;
		vector<string> tokens; // by the order of Explainer_record_config
		boost::split(tokens, line, boost::is_any_of("\t"));
		if (tokens.size() < 2)
			MTHROW_AND_ERR("Error Explainer_description_config::read_cfg_file - in file %s. Expected to have at least 2 tokens.\nGot line: \"%s\"\n",
						   file.c_str(), line.c_str());
		Explainer_record_config r;
		r.contributer_group_name = tokens[0];
		r.signal_name = tokens[1];
		if (tokens.size() > 2)
			r.max_count = med_stoi(tokens[2]);
		if (tokens.size() > 3)
			r.max_time_window = med_stoi(tokens[3]);
		if (tokens.size() > 4)
			r.time_channel = med_stoi(tokens[4]);
		if (tokens.size() > 5)
			r.time_unit = med_time_converter.string_to_type(tokens[5]);
		if (tokens.size() > 6)
			r.val_channel = med_stoi(tokens[6]);
		if (tokens.size() > 7)
			boost::split(r.sets, tokens[7], boost::is_any_of(","));

		if (contrib_seen.find(r.contributer_group_name) != contrib_seen.end())
			MTHROW_AND_ERR("Error, already defined rule for %s\n", r.contributer_group_name.c_str());
		contrib_seen.insert(r.contributer_group_name);

		records[r.contributer_group_name] = r;
	}

	file_reader.close();
}

void add_flag_response(nlohmann::ordered_json &js, float score, const MedAlgoMarkerInternal &ma,
					   const string &flag_threshold, float flag_threshold_numeric)
{
	if (!ma.has_threshold_settings()) // No flag settings - do nothing
		return;

	string err_msg;
	float cutoff = flag_threshold_numeric;
	if (flag_threshold_numeric == MED_MAT_MISSING_VALUE)
	{ // If not given take from flag_threshold
		cutoff = ma.fetch_threshold(flag_threshold, err_msg);
		js.push_back({"flag_threshold", flag_threshold});
	}

	if (!err_msg.empty())
	{
		js.push_back({"flag_result", AM_UNDEFINED_VALUE});
		add_to_json_array(js, "messages", err_msg);
		return;
	}
	js.push_back({"flag_threshold_numeric", cutoff});
	// All OK:
	int flag = int(score >= cutoff);
	js.push_back({"flag_result", flag});
}

void MedialInfraAlgoMarker::show_rep_data(char **response)
{
	MedPidRepository &rep = ma.get_rep();
	nlohmann::ordered_json js;
	js["data"] = nlohmann::ordered_json::array();
	unordered_map<int, unordered_set<string>> pid_to_signals;
	for (auto &it : rep.in_mem_rep.data)
	{
		const pair<int, int> &pid_sig = it.first;
		const string &sig_name = rep.sigs.name(pid_sig.second);
		pid_to_signals[pid_sig.first].insert(sig_name);
	}
	for (auto &it : pid_to_signals)
	{
		vector<string> all_sigs(it.second.begin(), it.second.end());
		js["data"] += {{"pid", it.first},
					   {"signals", all_sigs}};
	}

	json_to_char_ptr(js, response);
}