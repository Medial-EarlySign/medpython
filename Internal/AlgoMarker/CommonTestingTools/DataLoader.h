#ifndef __DATALOADER_H
#define __DATALOADER_H

#define AM_DLL_IMPORT

#include <AlgoMarker/AlgoMarker/AlgoMarker.h>
#include <AlgoMarker/DynAMWrapper/DynAMWrapper.h>
#include <AlgoMarker/CommonTestingTools/CommonTestingTools.h>
#include <Logger/Logger/Logger.h>
#include <MedProcessTools/MedProcessTools/MedModel.h>
#include <MedProcessTools/MedProcessTools/MedSamples.h>
#include <MedIO/MedIO/MedIO.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string.hpp>
//#include "internal_am.h"
#include <algorithm>

#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL	LOG_DEF_LEVEL
namespace CommonTestingTools {

	class DataLoader {
	public:
		static const int base_pid;

		MedModel model;
		MedSamples samples;
		MedPidRepository rep;
		vector<int> pids;
		vector<string> sigs;
		map<int, MedIdSamples* > pid2samples;
		map<string, vector<map<int, string> > > sig_dict_cached;

		void load(const string& rep_fname, const string& model_fname, const string& samples_fname = "", bool read_signals = true) {
			// read model file
			if (model.read_from_file(model_fname) < 0) {
				MERR("FAILED reading model file %s\n", model_fname.c_str());
				throw runtime_error(string("FAILED reading model file ") + model_fname);
			}

			unordered_set<string> sigs_set;
			model.get_required_signal_names(sigs_set);

			MLOG("Required signals:");
			for (auto &sig : sigs_set) {
				MLOG(" %s", sig.c_str());
				sigs.push_back(sig);
			}
			MLOG("\n");
			if (samples_fname != "") {
				if (samples.read_from_file(samples_fname)) {
					MERR("FAILED reading samples file %s\n", samples_fname.c_str());
					throw runtime_error(string("FAILED reading samples file ") + samples_fname);
				}
			}
			MLOG("\n");
			samples.get_ids(pids);
			if (read_signals) {
				if (rep.read_all(rep_fname, pids, sigs) < 0) {
					MERR("FAILED loading pids and signals from repository %s\n", rep_fname.c_str());
					throw runtime_error(string("FAILED loading pids and signals from repository"));
				}
			}
			else {
				if (rep.MedRepository::init(rep_fname) < 0) {
					MERR("Could not read repository definitions from %s\n", rep_fname.c_str());
					throw runtime_error(string("FAILED MedRepository::init(") + rep_fname + "\")");
				}
			}
			for (auto &id : samples.idSamples)
				pid2samples[id.id] = &id;
		}

		void get_sig_dict_cached(const string& cat_prefix = "", bool force_cat_prefix = false) {
			sig_dict_cached = get_sig_dict(cat_prefix, force_cat_prefix);
		}

		map<string, vector<map<int, string> > > get_sig_dict(const string& cat_prefix = "", bool force_cat_prefix = false) {
			map<string, vector<map<int, string> > > sig_dict;
			for (auto& sig : sigs) {
				vector<map<int, string > > chan_dict;
				int section_id = rep.dict.section_id(sig);
				int sid = rep.sigs.Name2Sid[sig];
				int n_vchan = rep.sigs.Sid2Info[sid].n_val_channels;
				for (int vchan = 0; vchan < n_vchan; ++vchan) {
					if (rep.sigs.is_categorical_channel(sig, vchan)) {
						map<int, string> new_dict;
						const auto& Id2Names = rep.dict.dict(section_id)->Id2Names;
						const auto& Member2Sets = rep.dict.dict(section_id)->Member2Sets;
						for (const auto& entry : Id2Names) {
							if (boost::starts_with(entry.second[0], cat_prefix)) {
								new_dict[entry.first] = entry.second[0];
								continue;
							}
							string new_ent = entry.second[0];
							if (Member2Sets.count(entry.first) != 0)
								for (const auto& setid : Member2Sets.at(entry.first)) {
									if (Id2Names.count(setid) != 0 && boost::starts_with(Id2Names.at(setid)[0], cat_prefix)) {
										if (!boost::starts_with(new_ent, cat_prefix) || new_ent.length() > Id2Names.at(setid)[0].length())
											new_ent = Id2Names.at(setid)[0];
									}
								}
							if (!force_cat_prefix || boost::starts_with(new_ent, cat_prefix))
								new_dict[entry.first] = new_ent;
						}

						chan_dict.push_back(new_dict);
					}
					else chan_dict.push_back(map<int, string>());
				}
				sig_dict[sig] = chan_dict;
			}
			return sig_dict;
		}

		map<string, vector<map<string, int>* > > get_sig_reverse_dict() {
			map<string, vector<map<string, int >* > > sig_dict;
			MLOG("(II)   Preparing signal reverse dictionary for signals\n");
			for (auto& sig : sigs) {
				//MLOG("(II)   Preparing signal dictionary for signal '%s'\n", sig.c_str());
				vector<map<string, int >* > chan_dict;
				if (rep.sigs.Name2Sid.count(sig) == 0) {
					MERR("no Name2Sid entry for signal '%s'\n", sig.c_str());
					exit(-1);
				}
				int section_id = rep.dict.section_id(sig);
				int sid = rep.sigs.Name2Sid[sig];
				int n_vchan = rep.sigs.Sid2Info[sid].n_val_channels;
				for (int vchan = 0; vchan < n_vchan; ++vchan) {
					if (rep.sigs.is_categorical_channel(sig, vchan))
					{
						chan_dict.push_back(&(rep.dict.dict(section_id)->Name2Id));
					}
					else {
						chan_dict.push_back(nullptr);
					}
				}
				sig_dict[sig] = chan_dict;
			}
			return sig_dict;
		}

		void export_required_data(const string& fname, const string& cat_prefix, bool force_cat_prefix) {
			ofstream outfile(fname, ios::binary | ios::out);

			MLOG("(II) Preparing dictinaries to export\n", fname.c_str());

			auto sig_dict = get_sig_dict(cat_prefix, force_cat_prefix);

			MLOG("(II) Exporting required data to %s\n", fname.c_str());

			UniversalSigVec usv;

			for (int pid : pids) {
				for (auto &sig : sigs) {
					rep.uget(pid, sig, usv);
					for (int i = 0; i < usv.len; ++i) {
						stringstream outss;
						outss << pid << '\t';
						outss << sig;
						for (int tchan = 0, n_tchan = usv.n_time_channels(); tchan < n_tchan; ++tchan) {
							outss << '\t' << usv.Time(i, tchan);
						}
						bool ignore_line = false;
						for (int vchan = 0, n_vchan = usv.n_val_channels(); vchan < n_vchan; ++vchan) {
							if (sig_dict.at(sig)[vchan].size() == 0)
								outss << '\t' << setprecision(10) << usv.Val(i, vchan);
							else {
								if (sig_dict.at(sig)[vchan].count((int)(usv.Val(i, vchan))) != 0) {
									outss << '\t' << sig_dict.at(sig)[vchan].at((int)(usv.Val(i, vchan)));
								}
								else {
									ignore_line = true;
								}
							}
						}
						if (!ignore_line)
							outfile << outss.str() << '\n';
					}
				}
			}
			outfile.close();
		}

		static void convert_reqfile_to_data(const string& input_json_fname, const string& output_data_fname) {
			ofstream outfile(output_data_fname, ios::binary | ios::out);
			ifstream infile(input_json_fname, ios::binary | ios::in);

			MLOG("(II) Exporting required data to %s\n", output_data_fname.c_str());

			json j;
			infile >> j;

			MLOG("(II) num of requests = %d\n", j.size());

			for (int pid = 0; pid < j.size(); ++pid) {
				json j_req_signals;
				if (j[pid].count("body") != 0)
					j_req_signals = j[pid]["body"]["signals"];
				else if (j[pid].count("signals") != 0)
					j_req_signals = j[pid]["signals"];
				else throw runtime_error("Unrecognized JSON fromat");

				for (const auto& j_sig : j_req_signals)
				{
					string sig = j_sig["code"];
					for (const auto& j_data : j_sig["data"]) {
						outfile << pid + base_pid << '\t';
						outfile << sig;
						for (const auto& j_time : j_data["timestamp"]) {
							outfile << '\t' << j_time;
						}
						for (const auto& j_val : j_data["value"]) {
							if (boost::to_upper_copy(sig) == "GENDER")
								outfile << '\t' << (boost::to_upper_copy(j_val.get<string>()) == "MALE" ? "1" : "2");
							else
								outfile << '\t' << j_val.get<string>();
						}

						outfile << "\n";
					}

				}
			}
			outfile.close();
		}



		void import_required_data(const string& fname) {
			ifstream infile(fname, ios::binary | ios::in);

			auto sig_dict = get_sig_reverse_dict();
			MLOG("(II)   Switching repo to in-mem mode\n");
			string curr_line;
			rep.switch_to_in_mem_mode();

			vector<int> tchan_vec;
			vector<float> vchan_vec;
			tchan_vec.reserve(10);
			vchan_vec.reserve(10);

			MLOG("(II)   reading data in to in-mem repository\n");
			int cur_line = 0;
			while (getline(infile, curr_line)) {
				cur_line++;
				if ((curr_line.size() > 1) && (curr_line[0] != '#')) {
					if (curr_line[curr_line.size() - 1] == '\r')
						curr_line.erase(curr_line.size() - 1);
					vector<string> fields;
					split(fields, curr_line, boost::is_any_of("\t"));
					int fields_i = 0;
					const auto& pid_str = fields[fields_i++];
					int pid = -1;
					try {
						pid = stoi(pid_str);
					}
					catch (...) {
						MERR("failed reading pid, performing stoi(\"%s\") at %s:%d\n", pid_str.c_str(), fname.c_str(), cur_line);
						exit(-1);
					}
					string sig = fields[fields_i++];
					int sid = rep.sigs.Name2Sid[sig];
					int n_vchan = rep.sigs.Sid2Info[sid].n_val_channels;
					int n_tchan = rep.sigs.Sid2Info[sid].n_time_channels;
					tchan_vec.clear();
					vchan_vec.clear();
					for (int tchan = 0; tchan < n_tchan; ++tchan) {
						const auto& field_str = fields[fields_i++];
						try {
							tchan_vec.push_back(stoi(field_str));
						}
						catch (...) {
							MERR("failed reading time channel #%d, performing stoi(\"%s\") at %s:%d\n", tchan, field_str.c_str(), fname.c_str(), cur_line);
							exit(-1);
						}

					}
					for (int vchan = 0; vchan < n_vchan; ++vchan) {
						if (sig_dict[sig][vchan] == nullptr) {
							const auto& field_str = fields[fields_i++];
							try {
								vchan_vec.push_back(stof(field_str));
							}
							catch (...) {
								MERR("failed reading value channel #%d, performing stof(\"%s\") at %s:%d\n", vchan, field_str.c_str(), fname.c_str(), cur_line);
								exit(-1);
							}
						}
						else
						{
							try {
								vchan_vec.push_back((*(sig_dict.at(sig)[vchan])).at(fields[fields_i++]));
							}
							catch (...) {
								MERR("Error converting sig %s, chan %d, '%s' back to code\n", sig.c_str(), vchan, fields[fields_i - 1].c_str());
								exit(-1);
							}
						}
					}
					rep.in_mem_rep.insertData(pid, sid, tchan_vec.data(), vchan_vec.data(), n_tchan, n_vchan);
				}
			}

			rep.in_mem_rep.sortData();

			infile.close();

			////REMOVE THIS
			//export_required_data("/nas1/Work/Users/Shlomi/apply-program/generated/repdata-re-export-after-import.txt", "ATC_", true);

		}

		void import_json_request_data(const string& fname) {
			ifstream infile(fname, ios::binary | ios::in);

			auto sig_dict = get_sig_reverse_dict();
			MLOG("(II)   Switching repo to in-mem mode\n");
			rep.switch_to_in_mem_mode();

			vector<int> tchan_vec;
			vector<float> vchan_vec;
			vector<string> vchan_vec_actual;
			tchan_vec.reserve(10);
			vchan_vec.reserve(10);

			MLOG("(II)   Started reading json data to in-mem repository\n");
			bool context = false;
			int cur_rec_no = 0;
			for (json j = read_json_array_next_chunk(infile, context); j != nullptr; j = read_json_array_next_chunk(infile, context)) {
				int pid = base_pid + cur_rec_no;
				if (j.count("body"))
					j = j["body"];
				if (j.count("signals") == 0 || !j.at("signals").is_array())
					MTHROW_AND_ERR("In file %s, Failed reading req #%d, no signals. request = \n'%s'\n", fname.c_str(), cur_rec_no, j.dump(1).c_str());
				for (auto j_sig : j.at("signals")) {
					string sig = j_sig.at("code");
					int sid = rep.sigs.Name2Sid[sig];
					int n_vchan = rep.sigs.Sid2Info[sid].n_val_channels;
					int n_tchan = rep.sigs.Sid2Info[sid].n_time_channels;
					tchan_vec.clear();
					vchan_vec.clear();
					vchan_vec_actual.clear();

					if (j_sig.count("data") == 0 || !j_sig.at("data").is_array())
						MTHROW_AND_ERR("In file %s, Failed reading data in req #%d . signal = '%s' json = \n'%s'\n\n", fname.c_str(), cur_rec_no, sig.c_str(), j_sig.dump(1).c_str());

					for (auto d_sig : j_sig.at("data")) {

						if (n_tchan > 0 && (d_sig.count("timestamp") == 0 || (!d_sig.at("timestamp").is_array())))
							MTHROW_AND_ERR("In file %s, Failed reading timestamp in req #%d . signal = '%s' json = \n'%s'\n\n", fname.c_str(), cur_rec_no, sig.c_str(), d_sig.dump(1).c_str());
						for (int tchan = 0; tchan < n_tchan; ++tchan) {
							string field_str = to_string(d_sig.at("timestamp")[tchan].get<int>());
							try {
								tchan_vec.push_back(stoi(field_str));
							}
							catch (...) {
								MERR("failed reading time channel #%d, performing stoi(\"%s\") at %s:%d\n", tchan, field_str.c_str(), fname.c_str(), cur_rec_no);
								exit(-1);
							}

						}
						if (n_vchan > 0 && (d_sig.count("value") == 0 || (!d_sig.at("value").is_array())))
							MTHROW_AND_ERR("In file %s, Failed reading value in req #%d . signal = '%s'\n", fname.c_str(), cur_rec_no, sig.c_str());
						for (int vchan = 0; vchan < n_vchan; ++vchan) {
							string field_str = d_sig.at("value")[vchan].get<string>();
							if (boost::to_upper_copy(sig) == "GENDER") {
								if (boost::to_upper_copy(field_str) == "MALE") field_str = "1";
								else if (boost::to_upper_copy(field_str) == "FEMALE") field_str = "2";
							}
							if (sig_dict[sig][vchan] == nullptr) {
								try {
									vchan_vec.push_back(stof(field_str));
									vchan_vec_actual.push_back(field_str);
								}
								catch (...) {
									MERR("failed reading value channel #%d, performing stof(\"%s\") at %s:%d\n", vchan, field_str.c_str(), fname.c_str(), cur_rec_no);
									exit(-1);
								}
							}
							else
							{
								try {
									vchan_vec.push_back((*(sig_dict.at(sig)[vchan])).at(field_str));
									vchan_vec_actual.push_back(field_str);
								}
								catch (...) {
									MERR("Error converting sig %s, chan %d, '%s' back to code in request #%d\n", sig.c_str(), vchan, field_str.c_str(), cur_rec_no);
									exit(-1);
								}
							}
						}
					}
					/* Write .data repo for testing
					int nelem = 0;
					if (tchan_vec.size() != 0)
					nelem = tchan_vec.size() / n_tchan;
					else nelem = vchan_vec.size() / n_vchan;
					int ti = 0;
					int vi = 0;
					for(int j = 0; j < nelem; j++) {
					MLOG("%d\t%s", pid, sig.c_str());
					for (int i = 0; i < n_tchan; i++)
					MLOG("\t%d", tchan_vec[ti++]);
					for (int i = 0; i < n_vchan; i++)
					MLOG("\t%s", vchan_vec_actual[vi++].c_str());
					MLOG("\n");
					}
					*/
					rep.in_mem_rep.insertData(pid, sid, tchan_vec.data(), vchan_vec.data(), (int)tchan_vec.size(), (int)vchan_vec.size());
				}
				cur_rec_no++;
			}

			rep.in_mem_rep.sortData();

			infile.close();
			MLOG("(II)   Finished loading json data to in-mem repository\n");
		}

		int load_samples_from_dates_to_score(const string& fname)
		{
			// read scores file
			vector<vector<string>> raw_scores;
			if (read_text_file_cols(fname, " \t", raw_scores) < 0) {
				MERR("Could not read scores file %s\n", fname.c_str());
				return -1;
			}
			MLOG("(II) Read %d lines from scores file %s\n", raw_scores.size(), fname.c_str());

			// prepare MedSamples
			for (auto &v : raw_scores)
				if (v.size() >= 2) {
					samples.insertRec(stoi(v[0]), stoi(v[1]));
				}
			samples.normalize();
			MLOG("(II) Prepared MedSamples\n");
			for (auto &id : samples.idSamples)
				pid2samples[id.id] = &id;
			return 0;
		}

		void am_add_data(AlgoMarker *am, int pid, int max_date, bool force_add_data, vector<string> ignore_sig, json& json_out) {
			static bool print_once = false;
			UniversalSigVec usv;
			int reserve_capacity = 100000;
			vector<long long> times;
			vector<float> vals;
			vector<bool> take_nelem;
			charpp_adaptor str_vals;
			times.reserve(reserve_capacity);
			vals.reserve(reserve_capacity);
			str_vals.reserve(reserve_capacity);
			if (!print_once) {
				print_once = true;
				MLOG("(INFO) force_add_data=%d\n", ((int)force_add_data));
				MLOG("(INFO) Will use %s API to insert data\n", (DynAM::so->addr_AM_API_AddDataStr == nullptr || force_add_data) ? "AddData()" : "AddDataStr()");
			}
			string reqId = string("req_") + to_string(pid) + "_" + to_string(max_date);
			json_out = json({});
			json_out["body"] = {
				{ "accountId", "A" },
			{ "requestId", reqId.c_str() },
			{ "customerId", "Earlysign" },
			{ "calculator" , "LC" },
			{ "signals",json::array() }
			};
			json_out["header"] = {
				{ "Accept", "application/json" },
			{ "Content-Type", "application/json" }
			};

			for (auto &sig : sigs) {
				if (std::find(ignore_sig.begin(), ignore_sig.end(), sig) != ignore_sig.end())
					continue;
				json json_sig;
				int sid = rep.sigs.Name2Sid[sig];
				//			int section_id = rep.dict.section_id(sig);
				usv.init(rep.sigs.Sid2Info[sid]);
				rep.uget(pid, sig, usv);
				int nelem = usv.len;
				if (nelem == 0)
					continue;
				vals.clear();
				times.clear();
				take_nelem.resize(nelem);

				if (usv.n_time_channels() <= 0) {
					std::fill(take_nelem.begin(), take_nelem.end(), true);
				}
				else {
					std::fill(take_nelem.begin(), take_nelem.end(), false);
					for (int i = 0; i < nelem; i++) {
						bool take_elem = true;
						for (int j = 0; j < usv.n_time_channels(); j++) {
							if (usv.Time(i, j) > max_date) {
								take_elem = false;
								break;
							}
						}
						if (take_elem) {
							for (int j = 0; j < usv.n_time_channels(); j++) {
								times.push_back((long long)usv.Time(i, j));
							}
							take_nelem[i] = true;
						}
						else {
							if (usv.n_time_channels() == 1)
								break;
						}
					}
				}

				if (DynAM::so->addr_AM_API_AddDataStr == nullptr || force_add_data) {
					vals.clear();
					if (usv.n_val_channels() > 0) {
						for (int i = 0; i < nelem; i++) {
							if (!take_nelem[i])
								continue;
							for (int j = 0; j < usv.n_val_channels(); j++) {
								vals.push_back(usv.Val(i, j));
							}

						}
					}

					if ((times.size() > 0) || (vals.size() > 0)) {
						get_volatile_data_adaptor<long long> p_times(times);
						get_volatile_data_adaptor<float> p_vals(vals);
						DynAM::AM_API_AddData(am, pid, sig.c_str(), (int)times.size(), p_times.get_volatile_data(), (int)vals.size(), p_vals.get_volatile_data());
						json_sig = json_AddData(sig.c_str(), (int)times.size(), p_times.get_volatile_data(), (int)vals.size(), p_vals.get_volatile_data(), usv.n_time_channels(), usv.n_val_channels());
					}
				}
				else {
					str_vals.clear();
					if (usv.n_val_channels() > 0) {
						for (int i = 0; i < nelem; i++) {
							if (!take_nelem[i])
								continue;
							for (int j = 0; j < usv.n_val_channels(); j++) {
								string val = "";
								if (rep.sigs.is_categorical_channel(sid, j)) {
									val = sig_dict_cached.at(sig)[j].at((int)(usv.Val(i, j)));
								}
								else {
									val = precision_float_to_string(usv.Val(i, j));
								}
								str_vals.push_back(val);
							}
						}
					}

					if ((times.size() > 0) || (str_vals.size() > 0)) {
						get_volatile_data_adaptor<long long> p_times(times);
						DynAM::AM_API_AddDataStr(am, pid, sig.c_str(), (int)times.size(), p_times.get_volatile_data(), (int)str_vals.size(), str_vals.get_charpp());
						json_sig = json_AddDataStr(sig.c_str(), (int)times.size(), p_times.get_volatile_data(), (int)vals.size(), str_vals.get_charpp(), usv.n_time_channels(), usv.n_val_channels());
					}
				}
				if (!json_sig.is_null())
					json_out["body"]["signals"].push_back(json_sig);
			}
		}

	};

}
#endif // __DATALOADER_H