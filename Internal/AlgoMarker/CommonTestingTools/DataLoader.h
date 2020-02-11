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

		void load(const string& rep_fname, const string& model_fname, const string& samples_fname = "", bool read_signals = true);

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



		void import_required_data(const string& fname);

		void import_json_request_data(const string& fname);

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

		void am_add_data(AlgoMarker *am, int pid, int max_date, bool force_add_data, vector<string> ignore_sig, json& json_out);

	};

}
#endif // __DATALOADER_H