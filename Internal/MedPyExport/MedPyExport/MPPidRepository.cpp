#include "MPPidRepository.h"
#include "MPSigExporter.h"
#include "MPDictionary.h"

#include <time.h>
#include <string>

#include "InfraMed/InfraMed/MedConvert.h"
#include "InfraMed/InfraMed/InfraMed.h"
#include "InfraMed/InfraMed/Utils.h"
#include "MedUtils/MedUtils/MedUtils.h"
#include "InfraMed/InfraMed/MedPidRepository.h"
#include "MedProcessTools/MedProcessTools/MedModel.h"
#include "MedProcessTools/MedProcessTools/SampleFilter.h"

#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL LOG_DEF_LEVEL

MPPidRepository::MPPidRepository() : o(new MedPidRepository()), dict(this)
{
#ifdef AM_API_FOR_CLIENT
	switch_to_in_mem();
#endif
};
MPPidRepository::~MPPidRepository()
{
	delete o;
	o = nullptr;
};
int MPPidRepository::dict_section_id(const std::string &secName) { return o->dict.section_id(secName); };
std::string MPPidRepository::dict_name(int section_id, int id) { return o->dict.name(section_id, id); };
#ifndef AM_API_FOR_CLIENT
int MPPidRepository::read_all(const std::string &conf_fname) { return val_or_exception(o->read_all(conf_fname), string("Error: read_all() failed")); };

int MPPidRepository::read_all_i(const std::string &conf_fname, const std::vector<int> &pids_to_take, const std::vector<int> &signals_to_take)
{
	return o->read_all(conf_fname, pids_to_take, signals_to_take);
};
/*
int MPPidRepository::read_all_s(const std::string &conf_fname, const std::vector<int> &pids_to_take, const vector<std::string> &signals_to_take)
{
	return o->read_all(conf_fname, pids_to_take, signals_to_take);
};*/

int MPPidRepository::read_all(const std::string &conf_fname, MEDPY_NP_INPUT(int *pids_to_take, unsigned long long num_pids_to_take), const std::vector<std::string> &signals_to_take)
{
	vector<int> pids_tt;
	buf_to_vector(pids_to_take, num_pids_to_take, pids_tt);
	return val_or_exception(o->read_all(conf_fname, pids_tt, signals_to_take), string("Error: read_all() failed"));
}

int MPPidRepository::loadsig(const std::string &signame) { return o->load(signame); };
int MPPidRepository::loadsig_pids(const std::string &signame, MEDPY_NP_INPUT(int *pids_to_take, unsigned long long num_pids_to_take))
{
	vector<int> pids(pids_to_take, pids_to_take + num_pids_to_take);
	return o->load(signame, pids);
};

#endif

int MPPidRepository::init(const std::string &conf_fname) { return o->init(conf_fname); };
const std::vector<int> &MPPidRepository::MEDPY_GET_pids() { return o->index.pids; };
int MPPidRepository::sig_id(const std::string &signame) { return o->dict.id(signame); };
int MPPidRepository::sig_type(const std::string &signame) { return o->sigs.type(signame); };
std::string MPPidRepository::sig_description(const std::string &signame)
{
	int id = o->sigs.sid(signame);
	if (id < 0)
		MTHROW_AND_ERR("Error signal %s wasn't found\n", signame.c_str());
	const SignalInfo &si = o->sigs.Sid2Info[id];
	UniversalSigVec usv;
	usv.init(si);
	string res = usv.get_signal_generic_spec();
	return res;
};
bool MPPidRepository::is_categorical(const std::string &signame, int val_channel)
{
	return o->sigs.is_categorical_channel(signame, val_channel);
}
MPSigVectorAdaptor MPPidRepository::uget(int pid, int sid)
{
	MPSigVectorAdaptor ret;
	o->uget(pid, sid, *((UniversalSigVec *)(ret.o)));
	return ret;
};

void MPPidRepository::finish_load_data()
{
	if (!o->in_mem_mode_active())
	{
		MWARN("WARN: not in mem_mode - doing nothing\n");
		return;
	}
	if (!data_load_sorted)
		o->in_mem_rep.sortData();
	data_load_sorted = true;
}

std::vector<bool> MPPidRepository::dict_prep_sets_lookup_table(int section_id, const std::vector<std::string> &set_names)
{
	vector<char> lut_cvec;
	o->dict.prep_sets_lookup_table(section_id, set_names, lut_cvec);
	vector<bool> lut_bvec;
	lut_bvec.reserve(lut_cvec.size());
	for (int i = 0; i < lut_cvec.size(); i++)
		lut_bvec.push_back(lut_cvec[i] != 0);
	// vector<bool> lut_bvec(lut_cvec.begin(), lut_cvec.end());
	return lut_bvec;
}

std::vector<bool> MPPidRepository::get_lut_from_regex(int section_id, const std::string &regex_s)
{
	vector<std::string> names;
	o->dict.dicts[section_id].get_regex_names(regex_s, names);
	return dict_prep_sets_lookup_table(section_id, names);
}

#ifndef AM_API_FOR_CLIENT
MPSigExporter MPPidRepository::export_to_numpy(string signame, MEDPY_NP_INPUT(int *pids_to_take, unsigned long long num_pids_to_take), int use_all_pids, int translate_flag, int free_sig, string filter_regex_str)
{
	return MPSigExporter(*this, signame, pids_to_take, num_pids_to_take, use_all_pids, translate_flag, free_sig, filter_regex_str);
}

int MPPidRepository::free(string signame)
{
	return o->free(signame);
}
#endif

void MPPidRepository::get_sig_structure(string &sig, int &n_time_channels, int &n_val_channels, int *&is_categ)
{
	int sid = o->sigs.sid(sig);
	if (sid <= 0)
	{
		n_time_channels = 0;
		n_val_channels = 0;
	}
	else
	{
		n_time_channels = o->sigs.Sid2Info[sid].n_time_channels;
		n_val_channels = o->sigs.Sid2Info[sid].n_val_channels;
		is_categ = &(o->sigs.Sid2Info[sid].is_categorical_per_val_channel[0]);
	}
}

int auto_time_convert(long long ts, int to_type)
{
	long long date_t = 0;
	long long hhmm = 0;

	// MLOG("auto time convert: Date is %d , ts %lld , to_type %d\n", MedTime::Date, ts, to_type);

	if ((ts / (long long)1000000000) == 0)
	{
		date_t = ts; // yyyymmdd
		hhmm = 0;
	}
	else if (((ts / (long long)100000000000) == 0))
	{
		date_t = ts / 100; // yyyymmddhh
		hhmm = 60 * (ts % 100);
	}
	else if (((ts / (long long)10000000000000) == 0))
	{
		date_t = ts / 10000; // yyyymmddhhmm
		hhmm = 60 * ((ts % 10000) / 100) + (ts % 100);
	}
	else
	{
		date_t = ts / 1000000; // yyyymmddhhmmss
		hhmm = 60 * ((ts % 1000000) / 10000) + ((ts % 10000) / 100);
	}

	// MLOG("auto_time_converter: ts %lld to_type %d data_t %lld hhmm %lld\n", ts, to_type, date_t, hhmm);

	if (to_type == MedTime::Date)
	{
		// MLOG("auto time convert: date_t %d\n", date_t);
		// Ensure valid date:
		int year = int(date_t / 10000);
		if (year < 1900 || year > 3000)
		{
			// MTHROW_AND_ERR("Error invalid date %lld\n", ts);
			return -1;
		}
		return (int)date_t;
	}

	if (to_type == MedTime::Minutes)
	{
		int year = int(date_t / 10000);
		if (year < 1900 || year > 2100)
		{
			// MTHROW_AND_ERR("Error invalid timestamp %lld\n", ts);
			return -1;
		}
		int minutes = med_time_converter.convert_date(MedTime::Minutes, (int)date_t);
		return minutes + (int)hhmm;
	}

	return 0;
}

void MPPidRepository::AddData(int patient_id, const char *signalName, int TimeStamps_len, long long *TimeStamps, int Values_len, float *Values)
{
	// At the moment MedialInfraAlgoMarker only loads timestamps given as ints.
	// This may change in the future as needed.
	int *i_times = NULL;
	vector<int> times_int;

	int tu = o->time_unit;
	if (TimeStamps_len > 0)
	{
		times_int.resize(TimeStamps_len);

		// currently assuming we only work with dates ... will have to change this when we'll move to other units
		for (int i = 0; i < TimeStamps_len; i++)
		{
			times_int[i] = auto_time_convert(TimeStamps[i], tu);
			if (times_int[i] < 0)
			{
				MERR("Error in AddData :: patient %d, signals %s, timestamp %lld is ilegal\n",
					 patient_id, signalName, TimeStamps[i]);
				return;
			}
			// MLOG("time convert %ld to %d\n", TimeStamps[i], times_int[i]);
		}
		i_times = &times_int[0];
	}

	int sid = o->sigs.Name2Sid[string(signalName)];
	if (sid < 0)
		MTHROW_AND_ERR("No signal %s\n", signalName); // no such signal
	if (o->in_mem_rep.insertData(patient_id, sid, i_times, Values, TimeStamps_len, Values_len) < 0)
		MTHROW_AND_ERR("Failed load data %s\n", signalName);
}

void MPPidRepository::AddDataStr(int patient_id, const char *signalName, int TimeStamps_len, long long *TimeStamps, int Values_len, char **Values)
{
	vector<float> converted_Values;
	MedRepository &rep = *o;

	try
	{
		string sig = signalName;
		int section_id = rep.dict.section_id(sig);
		int sid = rep.sigs.Name2Sid[sig];
		if (rep.sigs.Sid2Info[sid].n_val_channels > 0)
		{
			int Values_i = 0;
			const auto &category_map = rep.dict.dict(section_id)->Name2Id;
			int n_elem = (int)(Values_len / rep.sigs.Sid2Info[sid].n_val_channels);
			for (int i = 0; i < n_elem; i++)
			{
				for (int j = 0; j < rep.sigs.Sid2Info[sid].n_val_channels; j++)
				{
					float val = -1;
					if (!rep.sigs.is_categorical_channel(sid, j))
					{
						val = stof(Values[Values_i++]);
					}
					else
					{
						if (category_map.find(Values[Values_i]) == category_map.end())
						{
							MWARN("Found undefined code for signal \"%s\" and value \"%s\"\n",
								  sig.c_str(), Values[Values_i]);
						}
						else
							val = category_map.at(Values[Values_i]);
						++Values_i;
					}

					converted_Values.push_back(val);
				}
			}
		}
	}
	catch (...)
	{
		MTHROW_AND_ERR("Catched Error MedialInfraAlgoMarker::AddDataStr for signal %s in patient %d!!\n", signalName, patient_id);
		return;
	}

	AddData(patient_id, signalName, TimeStamps_len, TimeStamps, Values_len, converted_Values.data());
}

void MPPidRepository::switch_to_in_mem()
{
	if (!o->in_mem_mode_active())
	{
#ifndef AM_API_FOR_CLIENT
		MLOG("Switch to in mem repository\n");
#endif
		o->switch_to_in_mem_mode();
	}
}

int MPPidRepository::_load_single_json(void *_js)
{

	vector<string> messages;
	string current_time = "";
	bool good = true;
	bool mark_succ_ = false;
	json &j_data = *(json *)_js;
	json &js = *(json *)_js;
	int patient_id = -1;
	if (js.find("patient_id") != js.end())
		patient_id = js["patient_id"].get<long long>();
	else if (js.find("pid") != js.end())
		patient_id = js["pid"].get<long long>();

	if (patient_id < 0)
		MTHROW_AND_ERR("Error patient_id wasn't provided, should be provided with number bigger than zero.\n");

	try
	{
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
						MTHROW_AND_ERR("(330)Bad data json format - couldn't convert patient_id to integer");
					}
				}
				else
				{
					MTHROW_AND_ERR("(330)Bad data json format - patient_id suppose to be integer");
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
							MTHROW_AND_ERR("(330)Bad data json format - couldn't convert pid to integer");
						}
					}
					else
					{
						MTHROW_AND_ERR("(330)Bad data json format - pid suppose to be integer");
					}
				}
			}
		}
		if (patient_id <= 0)
		{
			MTHROW_AND_ERR("(330)Bad data json format - no patient_id was given");
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
				snprintf(buf, sizeof(buf), "Bad format in patient %d. Element should contain signals element as array",
						 patient_id);
			else
				snprintf(buf, sizeof(buf), "Bad format. Element should contain signals element as array");
			MTHROW_AND_ERR("%s", string(buf).c_str());
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
						snprintf(buf, sizeof(buf), "Bad format in patient %d. Element should contain code element as signal name",
								 patient_id);
					else
						snprintf(buf, sizeof(buf), "Bad format. Element should contain code element as signal name");
					messages.push_back(string(buf));
					get_current_time(current_time);
					MLOG("%s::%s\n", current_time.c_str(), buf);
					good = false;
					good_sig = false;
				}
				if (good_sig)
				{
					sig = s["code"].get<string>();
					int sid = o->sigs.sid(sig);
					get_sig_structure(sig, n_time_channels, n_val_channels, is_categ);
					if (n_time_channels == 0 && n_val_channels == 0)
					{
						char buf[5000];
						if (patient_id != 1)
							snprintf(buf, sizeof(buf), "An unknown signal was found: %s for patient %d",
									 sig.c_str(), patient_id);
						else
							snprintf(buf, sizeof(buf), "An unknown signal was found: %s",
									 sig.c_str());

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
							snprintf(buf, sizeof(buf), "Bad format for signal: %s in patient %d. No data element or data element is not array",
									 sig.c_str(), patient_id);
						else
							snprintf(buf, sizeof(buf), "Bad format for signal: %s. No data element or data element is not array",
									 sig.c_str());
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
									snprintf(buf, sizeof(buf), "Bad format for signal: %s in patient %d. timestamp should be array of timestamps, each represents a different channel.",
											 sig.c_str(), patient_id);
								else
									snprintf(buf, sizeof(buf), "Bad format for signal: %s. timestamp should be array of timestamps, each represents a different channel.",
											 sig.c_str());
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
									snprintf(buf, sizeof(buf), "Bad format for signal: %s in patient %d. value should be array of values, each represents a different channel.",
											 sig.c_str(), patient_id);
								else
									snprintf(buf, sizeof(buf), "Bad format for signal: %s. value should be array of values, each represents a different channel.",
											 sig.c_str());
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
											snprintf(buf, sizeof(buf), "Bad format for signal: %s in patient %d. Couldn't convert timestamp to integer",
													 sig.c_str(), patient_id);
										else
											snprintf(buf, sizeof(buf), "Bad format for signal: %s. Couldn't convert timestamp to integer",
													 sig.c_str());
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
										snprintf(buf, sizeof(buf), "Bad format for signal: %s in patient %d. timestamp element should be integer.",
												 sig.c_str(), patient_id);
									else
										snprintf(buf, sizeof(buf), "Bad format for signal: %s. timestamp element should be integer.",
												 sig.c_str());
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
									snprintf(buf, sizeof(buf), "Bad signal structure for signal: %s in patient %d. Received %d time channels instead of %d",
											 sig.c_str(), patient_id, nt, n_time_channels);
								else
									snprintf(buf, sizeof(buf), "Bad signal structure for signal: %s. Received %d time channels instead of %d",
											 sig.c_str(), nt, n_time_channels);
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
											snprintf(buf, sizeof(buf), "Bad format for signal: %s in patient %d. value element should be string.",
													 sig.c_str(), patient_id);
										else
											snprintf(buf, sizeof(buf), "Bad format for signal: %s. value element should be string.",
													 sig.c_str());
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
								string unit_m = o->sigs.unit_of_measurement(sid, nv);
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
												snprintf(buf, sizeof(buf), "Bad format for signal: %s in patient %d. value should be date format. Recieved %d.",
														 sig.c_str(), patient_id, full_date);
											else
												snprintf(buf, sizeof(buf), "Bad format for signal: %s. value should be date format. Recieved %d.",
														 sig.c_str(), full_date);
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
											snprintf(buf, sizeof(buf), "Bad format for signal: %s in patient %d. value should be date format. Recieved %s.",
													 sig.c_str(), patient_id, sv.c_str());
										else
											snprintf(buf, sizeof(buf), "Bad format for signal: %s. value should be date format. Recieved %s.",
													 sig.c_str(), sv.c_str());
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
									snprintf(buf, sizeof(buf), "Bad signal structure for signal: %s in patient %d. Received %d value instead of %d",
											 sig.c_str(), patient_id, nv, n_val_channels);
								else
									snprintf(buf, sizeof(buf), "Bad signal structure for signal: %s. Received %d value channels instead of %d",
											 sig.c_str(), nv, n_val_channels);
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
					try
					{
						AddDataStr(patient_id, sig.c_str(), n_times, p_times, n_vals, str_values);
						mark_succ_ = true;
					}
					catch (...)
					{
						char buf[5000];
						if (patient_id != 1)
							snprintf(buf, sizeof(buf), "General error in signal: %s for patient %d",
									 sig.c_str(), patient_id);
						else
							snprintf(buf, sizeof(buf), "General error in signal: %s",
									 sig.c_str());
						messages.push_back(string(buf));
						get_current_time(current_time);
						MLOG("%s::%s\n", current_time.c_str(), buf);
						good_sig = false;
						good = false;
						// return AM_FAIL_RC;
					}
				}
			}
		}
	}
	catch (...)
	{
		char buf[5000];
		snprintf(buf, sizeof(buf), "Bad data json format");
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
		MTHROW_AND_ERR("Please see errror messages\n");
	}

	// MLOG("Loading pid %d\n", patient_id);
	// Add to patient list
	unordered_set<int> s_pids(o->all_pids_list.begin(), o->all_pids_list.end());
	if (s_pids.find(patient_id) == s_pids.end())
	{
		o->all_pids_list.push_back(patient_id);
		o->pids.push_back(patient_id);
		o->index.pids.push_back(patient_id);
	}

	return patient_id;
}

void MPPidRepository::load_from_json(const std::string &json_file_path)
{
	switch_to_in_mem();
	if (o->sigs.Sid2Info.empty())
		MTHROW_AND_ERR("Error - please call init with repository config before\n");

	ifstream inf(json_file_path);
	if (!inf)
		MTHROW_AND_ERR("can't open json file [%s] for read\n", json_file_path.c_str());
	stringstream sstr;
	sstr << inf.rdbuf();
	inf.close();

	json j_data;
	try
	{
		j_data = json::parse(sstr);
	}
	catch (json::parse_error &err)
	{
		MTHROW_AND_ERR("Parsing error:\n%s", err.what());
	}
	catch (...)
	{
		MTHROW_AND_ERR("Error bad json format\n");
	}
	json *js = &j_data;
	data_load_sorted = false;

	if (j_data.find("multiple") != j_data.end())
	{
		for (auto &p_js : j_data["multiple"])
		{
			js = &p_js;
			if (p_js.find("body") != p_js.end())
				js = &p_js["body"];
			_load_single_json(js);
		}
	}
	else
	{ // single load
		if (j_data.find("body") != j_data.end())
			js = &j_data["body"];

		_load_single_json(js);
	}
}

void MPPidRepository::load_from_json_str(const std::string &json_content)
{
	switch_to_in_mem();
	if (o->sigs.Sid2Info.empty())
		MTHROW_AND_ERR("Error - please call init with repository config before\n");
	json j_data;
	try
	{
		j_data = json::parse(json_content);
	}
	catch (json::parse_error &err)
	{
		MTHROW_AND_ERR("Parsing error:\n%s", err.what());
	}
	catch (...)
	{
		MTHROW_AND_ERR("Error bad json format\n");
	}
	json *js = &j_data;
	data_load_sorted = false;

	if (j_data.find("multiple") != j_data.end())
	{
		for (auto &p_js : j_data["multiple"])
		{
			js = &p_js;
			if (p_js.find("body") != p_js.end())
				js = &p_js["body"];
			_load_single_json(js);
		}
	}
	else
	{ // single load
		if (j_data.find("body") != j_data.end())
			js = &j_data["body"];

		_load_single_json(js);
	}
}

void MPPidRepository::clear()
{
	if (!o->in_mem_mode_active())
		o->clear();
	else
	{
		// keep "init" and in_mem, just clear loaded data:
		o->in_mem_rep.clear();
		o->all_pids_list.clear();
		o->pids.clear();
		o->index.pids.clear();
	}
}

// ****************************** MPSig      *********************************
MPSig::MPSig(void *_o, int index) : o(_o), idx(index) {};
MPSig::MPSig(const MPSig &other)
{
	o = other.o;
	idx = other.idx;
};

std::vector<std::string> MPPidRepository::list_signals()
{
	std::vector<std::string> res;
	res.reserve(o->sigs.Name2Sid.size());
	for (const auto &it : o->sigs.Name2Sid)
		res.push_back(it.first);
	return res;
}

int MPSig::time(int chan) { return ((UniversalSigVec *)o)->Time(idx, chan); }
float MPSig::val(int chan) { return ((UniversalSigVec *)o)->Val(idx, chan); }
int MPSig::timeU(int to_time_unit) { return ((UniversalSigVec *)o)->TimeU(idx, to_time_unit); }
int MPSig::date(int chan) { return ((UniversalSigVec *)o)->Date(idx, chan); }
int MPSig::years(int chan) { return ((UniversalSigVec *)o)->Years(idx, chan); }
int MPSig::months(int chan) { return ((UniversalSigVec *)o)->Months(idx, chan); }
int MPSig::days(int chan) { return ((UniversalSigVec *)o)->Days(idx, chan); }
int MPSig::hours(int chan) { return ((UniversalSigVec *)o)->Hours(idx, chan); }
int MPSig::minutes(int chan) { return ((UniversalSigVec *)o)->Minutes(idx, chan); }

// ****************************** MPSigVectorAdaptor      *********************

MPSigVectorAdaptor::MPSigVectorAdaptor() { o = new UniversalSigVec(); };
MPSigVectorAdaptor::MPSigVectorAdaptor(const MPSigVectorAdaptor &other)
{
	o = new UniversalSigVec();
	*((UniversalSigVec *)(this->o)) = *((UniversalSigVec *)(other.o));
};
MPSigVectorAdaptor::~MPSigVectorAdaptor() { delete ((UniversalSigVec *)o); };

int MPSigVectorAdaptor::__len__() { return ((UniversalSigVec *)o)->len; };
MPSig MPSigVectorAdaptor::__getitem__(int i) { return MPSig(o, i); };

int MPSigVectorAdaptor::MEDPY_GET_type() { return ((UniversalSigVec *)o)->get_type(); }

int MPSigVectorAdaptor::MEDPY_GET_n_time_channels() { return ((UniversalSigVec *)o)->n_time_channels(); }
int MPSigVectorAdaptor::MEDPY_GET_n_val_channels() { return ((UniversalSigVec *)o)->n_val_channels(); }
int MPSigVectorAdaptor::MEDPY_GET_time_unit() { return ((UniversalSigVec *)o)->time_unit(); }
int MPSigVectorAdaptor::MEDPY_GET_size() { return (int)(((UniversalSigVec *)o)->size()); }

MPMedConvert::MPMedConvert(): o(new MedConvert()) {}
MPMedConvert::~MPMedConvert()
{
	delete o;
	o = nullptr;
}

void MPMedConvert::init_load_params(const std::string &load_args)
{
	o->init_load_params(load_args);
}

void MPMedConvert::create_rep(const std::string &conf_fname)
{
	o->create_rep(conf_fname);
}

int MPMedConvert::create_index(std::string &conf_fname)
{
	MedPidRepository mpr;

	if (mpr.read_config(conf_fname) < 0)
		return -1;
	if (mpr.read_pid_list() < 0)
		return -1;

	if (mpr.all_pids_list.size() == 0)
	{
		MERR("Flow: rep_create_pids: ERROR: got 0 pids in repository...., maybe pids list file was not created?...\n");
		return -1;
	}

	int jump = 1000000;
	int first = mpr.all_pids_list[0] - 1;
	int last = mpr.all_pids_list.back() + 1;
	if (((last - first) / 10) > jump)
		jump = (last - first) / 10;

	if (mpr.create(conf_fname, first, last, jump) < 0)
	{
		MERR("Flow: rep_create_pids: ERROR: Failed creating pid transpose rep.\n");
		return -1;
	}

	MLOG("Flow: rep_create_pids: Succeeded\n");
	return 0;
}