//
// helper routines for converting data into a new repository
//
#define __INFRAMED_DLL

#ifndef _SCL_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS
#endif

#include "Logger/Logger/Logger.h"
#include "MedConvert.h"
#include "Utils.h"
#include "MedUtils/MedUtils/MedUtils.h"
#include <stdio.h>
#include <fstream>
#include <algorithm>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/replace.hpp>


using namespace boost;

#define LOCAL_SECTION LOG_CONVERT
#define LOCAL_LEVEL	LOG_DEF_LEVEL
extern MedLogger global_logger;

//#define MAX_PID_TO_TAKE	1000
#define MAX_PID_TO_TAKE	100000000
//#define MAX_PID_TO_TAKE	5010000

void MedConvert::clear()
{
	config_fname = "";
	repository_config_fname = "";
	path = "";
	out_path = "";
	code_to_signal_fname = "";
	dict_fnames.clear();
	sig_fnames.clear();
	registry_fname = "";
	relative = 0;
	in_data_fnames.clear();
	in_strings_data_fnames.clear();
	prefix_names.clear();
	index_fnames.clear();
	data_fnames.clear();
	codes2names.clear();
	dict.clear();
	sid2fno.clear();
	sid2serial.clear();
	serial2sid.clear();
	forced.clear();
	safe_mode = 0;
}

//------------------------------------------------
int MedConvert::read_config(const string &fname)
{
	ifstream inf(fname);

	if (!inf) {
		MERR("MedConvert: read_config: Can't open file [%s]\n", fname.c_str());
		return -1;
	}

	clear();

	string curr_line;
	while (getline(inf, curr_line)) {
		if ((curr_line.size() > 1) && (curr_line[0] != '#')) {

			vector<string> fields;
			split(fields, curr_line, boost::is_any_of(" \t"));

			if (fields.size() >= 2) {

				if (fields[0].compare("DIR") == 0) {
					path = fields[1];
					if (path.compare(".") == 0) {
						// in this case we fix our path to be from where we were called
						size_t found = fname.find_last_of("/\\");
						path = fname.substr(0, found);
					}
				}
				if (fields[0].compare("OUTDIR") == 0) out_path = fields[1];
				if (fields[0].compare("CONFIG") == 0) repository_config_fname = fields[1];
				if (fields[0].compare("DICTIONARY") == 0) dict_fnames.push_back(fields[1]);
				if (fields[0].compare("SIGNAL") == 0) { sig_fnames.push_back(fields[1]); dict_fnames.push_back(fields[1]); }
				if (fields[0].compare("CODES") == 0) code_to_signal_fname = fields[1];
				if (fields[0].compare("FNAMES") == 0) prefixes_fname = fields[1];
				if (fields[0].compare("SFILES") == 0) signal_to_files_fname = fields[1];
				if (fields[0].compare("REGISTRY") == 0) registry_fname = fields[1];
				if (fields[0].compare("DATA") == 0) in_data_fnames.push_back(fields[1]);
				if (fields[0].compare("DATA_S") == 0) in_strings_data_fnames.push_back(fields[1]);
				if (fields[0].compare("MODE") == 0) mode = stoi(fields[1]);
				if (fields[0].compare("SAFE_MODE") == 0) safe_mode = stoi(fields[1]);
				if (fields[0].compare("PREFIX") == 0) rep_files_prefix = fields[1];
				if (fields[0].compare("RELATIVE") == 0) relative = 1;
				if (fields[0].compare("DESCRIPTION") == 0) description = fields[1];
				if (fields[0].compare("FORCE_SIGNAL") == 0) {
					vector<string> fsigs;
					split(fsigs, fields[1], boost::is_any_of(","));
					for (int i = 0; i < fsigs.size(); i++) {
						MLOG_D("MedConvert: Will force signal %s\n", fsigs[i].c_str());
						forced.push_back(fsigs[i]);
					}
				}
				if (fields[0].compare("LOAD_ONLY") == 0) {
					vector<string> fsigs;
					split(fsigs, fields[1], boost::is_any_of(","));
					for (int i = 0; i < fsigs.size(); i++) {
						MLOG_D("MedConvert: Will load only signal %s\n", fsigs[i].c_str());
						load_only.push_back(fsigs[i]);
					}
				}
			}
		}
	}

	inf.close();
	return 0;
}

//------------------------------------------------
int MedConvert::read_code_to_signal(const string &fname)
{
	ifstream inf(fname);

	if (!inf) {
		MERR("MedConvert: read_code_to_signal: Can't open file %s\n", fname.c_str());
		return -1;
	}

	string curr_line;
	while (getline(inf, curr_line)) {
		if ((curr_line.size() > 1) && (curr_line[0] != '#')) {

			vector<string> fields;
			split(fields, curr_line, boost::is_any_of(" \t"));

			if (fields.size() >= 2) {
				//MLOG("Code[%s] = %s\n",fields[0].c_str(),fields[1].c_str());
				codes2names[fields[0]] = fields[1];
			}
		}
	}

	inf.close();
	return 0;

}

//------------------------------------------------
int MedConvert::read_prefix_names(const string &fname)
{
	ifstream inf(fname);

	if (!inf) {
		MERR("MedConvert: read_prefix_names: Can't open file %s\n", fname.c_str());
		return -1;
	}

	prefix_names.clear();
	string curr_line;
	while (getline(inf, curr_line)) {
		MLOG("read_prefix_name(): file %s line %s\n", fname.c_str(), curr_line.c_str());
		if ((curr_line.size() > 1) && (curr_line[0] != '#')) {

			vector<string> fields;
			split(fields, curr_line, boost::is_any_of(" \t"));

			if (fields.size() >= 2) {
				int fno = stoi(fields[0]);
				if (prefix_names.size() < fno + 1)
					prefix_names.resize(fno + 1);
				prefix_names[fno] = fields[1];
			}
		}
	}

	inf.close();
	MLOG("Finished reading prefix names file %s , got %d prefixes\n", fname.c_str(), prefix_names.size());
	return 0;
}

//------------------------------------------------
// assumes dictionary is already loaded !
int MedConvert::read_signal_to_files(const string &fname)
{
	ifstream inf(fname);

	if (!inf) {
		MERR("MedConvert: read_signal_to_files: Can't open file %s\n", fname.c_str());
		return -1;
	}

	string curr_line;
	int sid;
	int serial = 0;
	serial2siginfo.clear();
	while (getline(inf, curr_line)) {
		if ((curr_line.size() > 1) && (curr_line[0] != '#')) {

			vector<string> fields;
			split(fields, curr_line, boost::is_any_of(" \t"));

			if (fields.size() >= 2) {
				int fno = stoi(fields[0]);
				sid = dict.id(fields[1]);
				if (sid >= 0) {
					sid2fno[sid] = fno;
					serial2sid.push_back(sid);
					sig_info si;
					si.fno = fno;
					si.serial = (int)serial2siginfo.size();
					si.type = sigs.type(sid);
					si.sid = sid;
					serial2siginfo.push_back(si);
					sid2serial[sid] = serial++;
					MLOG("fno %d sig %s sid %d sid2serial %d sid3fno %d\n", fno, fields[1].c_str(), sid, sid2serial[sid], sid2fno[sid]);
				}
			}
		}
	}

	inf.close();

	return 0;
}
//------------------------------------------------
int MedConvert::prep_sids_to_load()
{
	if (load_only.size() > 0)
		load_only.insert(load_only.end(), forced.begin(), forced.end());
	sids_to_load.resize(MAX_SID_NUMBER, 0);
	if (mode >= 2 && load_only.size() > 0) {
		for (int i = 0; i < load_only.size(); i++) {
			int sid = sigs.sid(load_only[i]);
			if (sid <= 0) {
				MERR("ERROR: asked to load a non defined signal: %s\n", load_only[i].c_str());
				return -1;
			}
			sids_to_load[sid] = 1;
		}

	}
	else {
		for (int i = 0; i < sigs.signals_ids.size(); i++)
			sids_to_load[sigs.signals_ids[i]] = 1;
	}

	return 0;
}

//------------------------------------------------
int MedConvert::read_all(const string &config_fname)
{
#if defined (_MSC_VER) || defined (_WIN32)
	MLOG("Max open files %d\n", _getmaxstdio());
	_setmaxstdio(4096);
	MLOG("Max open files raised to %d\n", _getmaxstdio());
#endif

	if (read_config(config_fname) < 0) {
		MERR("MedConvert: read_all: read_config %s failed\n", config_fname.c_str());
		return -1;
	}

	MLOG("MedConvert: read_all: read config file\n");


	if (path.length() == 0)
		path = ".";

	if (out_path.length() == 0)
		out_path = path;

	// add path to all input fnames + fix names
	add_path_to_name_IM(path, code_to_signal_fname);
	add_path_to_name_IM(path, signal_to_files_fname);
	add_path_to_name_IM(path, in_data_fnames);
	add_path_to_name_IM(path, in_strings_data_fnames);
	add_path_to_name_IM(path, prefixes_fname);

	if (registry_fname != "")
		add_path_to_name_IM(path, registry_fname);


	// read dictionary
	if (dict.read(path, dict_fnames) < 0) {
		return -1;
	}

	MLOG("MedConvert: read_all: read dictionary files\n");

	// read signals
	if (sigs.read(path, sig_fnames) < 0) {
		return -1;
	}

	MLOG("MedConvert: read_all: read signal files\n");

	// now add as default all sigs name to their own
	for (auto& sig : sigs.signals_names)
		codes2names[sig] = sig;

	// read signal to file
	if (code_to_signal_fname != "" && read_code_to_signal(code_to_signal_fname) < 0) {
		return -1;
	}

	MLOG("MedConvert: read_all: read code_to_signal file [%s]\n", code_to_signal_fname.c_str());


	// mode 2 and up supports loading a subset of signals ! , older modes will always try to load all
	if (prep_sids_to_load() < 0)
		return -1;

	if (mode < 2) {
		// read prefix names
		if (prefixes_fname != "" && read_prefix_names(prefixes_fname) < 0) {
			return -1;
		}

		MLOG("MedConvert: read_all: read prefix_names file\n");

		// read maping of signals to output files (and build serial numbers for signals)
		if (signal_to_files_fname != "" && read_signal_to_files(signal_to_files_fname) < 0) {
			return -1;
		}

		MLOG("MedConvert: read_all: read signal_to_files file\n");

	}
	else {

		// in mode 2 we generate the prefix names on our own, one for each signal
		// and also generate the mapping from each signal to its file.
		generate_prefix_names();
	}

	// create and add path to all output file names
	index_fnames.resize(prefix_names.size());
	data_fnames.resize(prefix_names.size());
	for (int i = 0; i < prefix_names.size(); i++) {
		if (prefix_names[i] != "") {
			index_fnames[i] = prefix_names[i] + ".idx";
			data_fnames[i] = prefix_names[i] + ".data";

			MLOG("i=%d index %s data %s\n", i, index_fnames[i].c_str(), data_fnames[i].c_str());
		}
	}
	if (add_path_to_name_IM(out_path, repository_config_fname) == -1)
		return -1;

	// Create repository config file
	if (create_repository_config() < 0)
		MTHROW_AND_ERR("MedConvert: read_all(): failed generating repository config file\n");

	// Copy sig and dict file to output directory
	if (copy_files_IM(path, out_path, sig_fnames) < 0 || copy_files_IM(path, out_path, dict_fnames) < 0)
		MTHROW_AND_ERR("MedConvert : read_all() : failed copying files from in to out directory\n");

	// add path to more files + fix paths
	if (add_path_to_name_IM(path, sig_fnames) == -1 ||
		add_path_to_name_IM(path, dict_fnames) == -1 ||
		add_path_to_name_IM(out_path, index_fnames) == -1 ||
		add_path_to_name_IM(out_path, data_fnames) == -1)
		return -1;

	MLOG("MedConvert: read_all: prepared names\n");


	// actually do the work
	if (create_indexes() < 0)
		MTHROW_AND_ERR("MedConvert: read_all(): failed generating new data and indexes\n");

	return 0;
}

//------------------------------------------------
int MedConvert::get_next_signal(ifstream &inf, int file_type, pid_data &curr, int &fpid, file_stat &curr_fstat, map<pair<string, string>, int>& missing_dict_vals)
{
	if (!inf.is_open())
		return 0; // file is closed nothing to do

	string curr_line;
	streampos pos;

	bool get_next = true;
	collected_data cd;
	string vfield, vfield2;
	int i;
	int sid, vid;
	while (get_next) {
		pos = inf.tellg();
		if (getline(inf, curr_line)) {
			if (curr_line[curr_line.size() - 1] == '\r')
				curr_line.erase(curr_line.size() - 1);
			boost::trim(curr_line);
			curr_fstat.n_lines++;
			if ((curr_line.size() > 1) && (curr_line[0] != '#')) {
				//if (fpid == 5025392)
				//	MLOG("fpid %d : %s\n", fpid, curr_line.c_str());
				vector<string> fields;
				if (file_type == 1 || file_type == 3)
					split(fields, curr_line, boost::is_any_of("\t"));
				else
					split(fields, curr_line, boost::is_any_of(" \t"));
				curr_fstat.n_relevant_lines++;
				//if (fpid == 5025392)
	//				MLOG("working on: (fpid %d) (curr.pid %d) (file_type %d) (f[0] %s) (nfields %d) ##>%s<##\n",fpid,curr.pid,file_type,fields[0].c_str(),fields.size(),curr_line.c_str());
				//if (fields.size() > 6) MLOG("WEIRD f[6]= ##>%s<##\n",fields[5].c_str());
				if (((file_type == 1) && (fields.size() == 4)) ||
					((file_type == 2) && (fields.size() >= 3)) ||
					((file_type == 3) && (fields.size() >= 3))) {

					int line_pid;
					try {
						line_pid = stoi(fields[0]);
					}
					catch (...) {
						MERR("ERROR: bad format in file %s with first token of pid, in line %d:\n%s\n",
							curr_fstat.fname.c_str(), curr_fstat.n_parsed_lines, curr_line.c_str());
						throw;
					}
					//if (fpid == 5025392)
					//	MLOG("working on: (fpid %d) (curr.pid %d) (file_type %d) (line_pid %d) %s\n",fpid,curr.pid,file_type,line_pid,curr_line.c_str());
					if (line_pid == curr.pid) {
						cd.zero();
						if (file_type == 1) {

							// Registry file //format: pid , stage(string) , date, location(number)	// tab delimited 

							try {
								// Cancer_Location
								i = sid2serial[dict.id(string("Cancer_Location"))];
								cd.date = stoi(fields[2]);
								cd.val = (float)(dict.id(fields[3]));
								curr.raw_data[i].push_back(cd);

								// Cancer_Stage
								i = sid2serial[dict.id(string("Cancer_Stage"))];
								cd.date = stoi(fields[2]);
								cd.val = (float)(stoi(fields[1]));
								curr.raw_data[i].push_back(cd);

								curr_fstat.n_parsed_lines++;
							}
							catch (...) {
								MERR("ERROR: bad format in parsing registry file %s in line %d:\n%s\n",
									curr_fstat.fname.c_str(), curr_fstat.n_parsed_lines, curr_line.c_str());
								throw;
							}

						}
						else if (file_type == 2) {
							// regular data file // format: pid , signal_code, date/time/range , value
							if (codes2names.find(fields[1]) != codes2names.end()) {
								sid = dict.id(codes2names[fields[1]]);
								if (sid >= 0 && sids_to_load[sid]) {
									try {
										i = sid2serial[sid];
										switch (sigs.type(sid)) {

										case T_Value:
											if (fields.size() == 3)
												cd.val = med_stof(fields[2]);
											else
												cd.val = med_stof(fields[3]); // backward compatible with date 0 trick to load value only data
											break;

										case T_DateVal:
											if (fields.size() < 4) {
												MERR("Convert ERROR : line with too few fields :: %s\n", curr_line.c_str());
												exit(-1);
											}
											cd.date = med_stoi(fields[2]);
											cd.val = med_stof(fields[3]);
											break;

										case T_DateRangeVal:
											cd.date = med_stoi(fields[2]);
											cd.date2 = med_stoi(fields[3]);
											cd.val = med_stof(fields[4]);
											break;

										case T_TimeVal:
											cd.time = stoll(fields[2]);
											cd.val = med_stof(fields[3]);
											break;

										case T_TimeRangeVal:
											cd.time = stoll(fields[2]);
											cd.time2 = stoll(fields[3]);
											cd.val = stof(fields[4]);
											break;

										case T_TimeStamp:
											cd.time = stoll(fields[2]);
											break;

										case T_DateVal2:
											cd.date = med_stoi(fields[2]);
											cd.val = med_stof(fields[3]);
											cd.val2 = (unsigned short)med_stoi(fields[4]);
											break;

										case T_TimeLongVal:
											cd.time = stoll(fields[2]);
											cd.longVal = stoll(fields[3]);
											break;

										case T_DateShort2:
											cd.date = med_stoi(fields[2]);
											cd.val1 = (short)med_stoi(fields[3]);
											cd.val2 = (short)med_stoi(fields[4]);
											break;

										case T_ValShort2:
											cd.val1 = (short)med_stoi(fields[2]);
											cd.val2 = (short)med_stoi(fields[3]);
											break;

										case T_ValShort4:
											cd.val1 = (short)med_stoi(fields[2]);
											cd.val2 = (short)med_stoi(fields[3]);
											cd.val3 = (short)med_stoi(fields[4]);
											cd.val4 = (short)med_stoi(fields[5]);
											break;

										case T_CompactDateVal:
											cd.date = (int)med_stoi(fields[2]);
											cd.val1 = (unsigned short)med_stoi(fields[3]);
											break;

										default:
											MERR("MedConvert: get_next_signal: unknown signal type for sid %d\n", sid);
											return -1;
										}

										curr_fstat.n_parsed_lines++;
										curr.raw_data[i].push_back(cd);
									}
									catch (...) {
										curr_fstat.n_bad_format_lines++;
										if (curr_fstat.n_bad_format_lines < 10)
											MWARN("MedConvert: WARNING: bad format in parsing DATA file %s with type=%d in line %d:\n%s\n",
												curr_fstat.fname.c_str(), sigs.type(sid)
												, curr_fstat.n_parsed_lines, curr_line.c_str());
									}
								}
							}
							else {
								if (safe_mode) {
									MERR("MedConvert: ERROR: unrecognized signal name %s (need to add to codes_to_signals file) in file %s :: curr_line is %s\n",
										fields[1].c_str(), curr_fstat.fname.c_str(), curr_line.c_str());
									return -1;
								}
							}
						}
						else if (file_type == 3) {
							// String-valued data file // format : pid , signal_code, date/time , string in dictionary
							if (codes2names.find(fields[1]) != codes2names.end()) {
								sid = dict.id(codes2names[fields[1]]);
								if (sid >= 0 && sids_to_load[sid]) {
									try {
										i = sid2serial[sid];
										switch (sigs.type(sid)) {

										case T_Value:
											cd.date = 0;
											vfield = fields[2];
											break;

										case T_DateVal:
											cd.date = med_stoi(fields[2]);
											vfield = fields[3];
											break;

										case T_DateRangeVal:
											cd.date = med_stoi(fields[2]);
											cd.date2 = med_stoi(fields[3]);
											vfield = fields[4];
											break;

										case T_TimeVal:
											cd.time = stoll(fields[2]);
											vfield = fields[3];
											break;

										case T_TimeRangeVal:
											cd.time = stoll(fields[2]);
											cd.time2 = stoll(fields[3]);
											vfield = fields[4];
											break;

										case T_TimeLongVal:
											cd.time = stoll(fields[2]);
											vfield = fields[3];
											break;

										case T_TimeStamp:
											cd.time = stoll(fields[2]);
											break;
										case T_DateShort2:
											cd.date = med_stoi(fields[2]);
											vfield = fields[3];
											vfield2 = fields[4];
											break;

										case T_DateVal2:
											// NOT SUPPORTED !!!!!
											MERR("This type is NOT supported in string input mode yet !!!!!\n");
											break;

										default:
											MERR("MedConvert: get_next_signal: unknown signal type for sid %d\n", sid);
											return -1;

										}

										int section = dict.section_id(sigs.name(sid));
										vid = dict.id(section, vfield);
										int vid2 = dict.id(section, vfield2);
										if (vid >= 0) {
											if (sigs.type(sid) == T_TimeLongVal)
												cd.longVal = (long long)vid;
											else if (sigs.type(sid) == T_DateShort2) {
												cd.val1 = (short)vid;
												cd.val2 = (short)vid2;
											}
											else
												cd.val = (float)vid;
											curr.raw_data[i].push_back(cd);
											curr_fstat.n_parsed_lines++;

										}
										else {
											pair<string, string> my_key = make_pair(sigs.name(sid), vfield);
											if (missing_dict_vals.find(my_key) == missing_dict_vals.end()) {
												MWARN("MedConvert::get_next_signal: signal string [%s] is missing from dictionary (sig [%s]) : file [%s] : line [%s] \n",
													vfield.c_str(), sigs.name(sid).c_str(), curr_fstat.fname.c_str(), curr_line.c_str());
												missing_dict_vals[my_key] = 1;
											}
											else
												missing_dict_vals[my_key]++;
										}
									}
									catch (...) {
										MERR("ERROR: bad format in parsing DATA_S file %s (file_type=%d) in line %d:\n%s\n",
											curr_fstat.fname.c_str(), file_type, curr_fstat.n_parsed_lines, curr_line.c_str());
										throw;
									}
								}
							}
							else {
								if (safe_mode) {
									MTHROW_AND_ERR("MedConvert: ERROR: unrecognized signal name %s (need to add to codes_to_signals file) in file %s :: curr_line is %s\n",
										fields[1].c_str(), curr_fstat.fname.c_str(), curr_line.c_str());
								}

							}
						}
					}
					else if (line_pid < fpid) {
						MWARN("MedConvert: get_next_signal: fpid is %d , but got line: %s\n", fpid, curr_line.c_str());
						if (safe_mode) {
							MERR("MedConvert: ERROR: file %s seems to be not sorted by pid\n", curr_fstat.fname.c_str());
							return -1;
						}
					}
					else {
						fpid = line_pid;
						inf.seekg(pos, ios::beg); // roll file back to the start of curr line
						curr_fstat.n_lines--;
						curr_fstat.n_relevant_lines--;
						get_next = false;
					}
				}
			}
		}
		else
			get_next = false;
	}
	if (inf.eof() || (fpid > MAX_PID_TO_TAKE)) {
		fpid = -1;
		inf.close();
		n_open_in_files--;
	}

	return 0;
}

//------------------------------------------------
int MedConvert::create_repository_config()
{
	if (repository_config_fname == "") {
		MWARN("No repository_config (CONFIG) file specified, not creating it\n");
		return 1;
	}
	// Open output file
	repository_config_f.open(repository_config_fname, ios::out);
	if (!repository_config_f) {
		MERR("MedConvert:: create_repository_config:: can't open output file [%s]\n", repository_config_fname.c_str());
		return -1;
	}

	if (description != "")
		repository_config_f << "DESCRIPTION\t" << description << endl;

	if (relative)
		repository_config_f << "DIR\t." << endl;
	else
		repository_config_f << "DIR\t" << out_path << endl;

	for (unsigned int i = 0; i < dict_fnames.size(); i++)
		repository_config_f << "DICTIONARY\t" << dict_fnames[i].c_str() << endl;

	for (unsigned int i = 0; i < sig_fnames.size(); i++)
		repository_config_f << "SIGNAL\t" << sig_fnames[i].c_str() << endl;

	repository_config_f << "MODE\t" << mode << endl;
	if (mode < 3) {
		for (unsigned int i = 0; i < data_fnames.size(); i++)
			repository_config_f << "DATA\t" << i << "\t" << data_fnames[i].c_str() << endl;

		for (unsigned int i = 0; i < index_fnames.size(); i++)
			repository_config_f << "INDEX\t" << index_fnames[i].c_str() << endl;
	}
	else {
		repository_config_f << "PREFIX\t" << rep_files_prefix.c_str() << endl;
	}

	repository_config_f.close();
	return 0;
}
//------------------------------------------------
int MedConvert::create_indexes()
{
	int i;
	pid_data curr;

	int n_files = (int)in_data_fnames.size() + (int)in_strings_data_fnames.size() + 1; // all input data files  + registry
	vector<ifstream> infs(n_files);
	fstats.resize(n_files);
	vector<int> file_type;

	pid_in_file.resize(n_files);
	file_type.resize(n_files);

	fill(pid_in_file.begin(), pid_in_file.end(), -1);

	// open all files

	n_open_in_files = 0;

	// registry
	if (registry_fname != "") {
		infs[n_open_in_files].open(registry_fname, ios::in | ios::binary);
		if (!infs[n_open_in_files]) {
			MERR("MedConvert: create_indexes: can't open registry file %s\n", registry_fname.c_str());
			return -1;
		}
		file_type[n_open_in_files] = 1;
		fstats[n_open_in_files].fname = registry_fname;
		fstats[n_open_in_files].id = n_open_in_files;
		n_open_in_files++;
	}

	// all data files
	for (i = 0; i < in_data_fnames.size(); i++) {
		if (in_data_fnames[i] != "") {
			infs[n_open_in_files].open(in_data_fnames[i], ios::in | ios::binary);
			if (!infs[n_open_in_files]) {
				MERR("MedConvert: create_indexes: can't open input data file %s\n", in_data_fnames[i].c_str());
				return -1;
			}
		}
		file_type[n_open_in_files] = 2;
		fstats[n_open_in_files].fname = in_data_fnames[i];
		fstats[n_open_in_files].id = n_open_in_files;
		MLOG("MedConvert: opened file %s for input file (%d) , of type %d\n", in_data_fnames[i].c_str(), n_open_in_files, file_type[n_open_in_files]);
		n_open_in_files++;
	}

	for (i = 0; i < in_strings_data_fnames.size(); i++) {
		if (in_strings_data_fnames[i] != "") {
			infs[n_open_in_files].open(in_strings_data_fnames[i], ios::in | ios::binary);
			if (!infs[n_open_in_files]) {
				MERR("MedConvert: create_indexes: can't open input data file %s\n", in_strings_data_fnames[i].c_str());
				return -1;
			}
		}
		file_type[n_open_in_files] = 3;
		fstats[n_open_in_files].fname = in_strings_data_fnames[i];
		fstats[n_open_in_files].id = n_open_in_files;
		MLOG("MedConvert: opened file %s for input file (%d) , of type %d\n", in_strings_data_fnames[i].c_str(), n_open_in_files, file_type[n_open_in_files]);
		n_open_in_files++;
	}

	int c_pid = -1;
	int n_files_opened = n_open_in_files;

	MLOG("MedConvert: create_indexes: n_open_in_files %d\n", n_open_in_files);

	if (open_indexes() < 0) {
		MERR("MedConvert: create_indexes: couldn't open index and data files\n");
		return -1;
	}

	int n_pids_extracted = 0;
	map<pair<string, string>, int> missing_dict_vals;
	vector<int> all_pids;  // a list of all pids in the repository to be written to file.
	all_pids.push_back(0); // reserved place for later placing of total number of pids
	time_t start = time(NULL);
	time_t last_time_print = start;
	while (n_open_in_files > 0) {

		// find current pid to extract
		c_pid = -1;
		for (i = 0; i < n_files_opened; i++) {
			if (pid_in_file[i] > 0 && c_pid < 0) c_pid = pid_in_file[i];
			if (c_pid >= 0 && pid_in_file[i] > 0 && pid_in_file[i] < c_pid) c_pid = pid_in_file[i];
		}

		if (c_pid % 100000 == 0)
			MLOG("Current pid to extract is %d <<<<< >>>>> n_extracted %d n_open_in_files %d\n", c_pid, n_pids_extracted, n_open_in_files);

		// read data from files
		curr.raw_data.clear();
		curr.raw_data.resize(serial2sid.size());
		curr.pid = c_pid;



		for (i = 0; i < n_files_opened; i++) {
			int fpid = c_pid;
			if (infs[i].is_open() && pid_in_file[i] <= c_pid) {
				//MLOG("file %d :: pid_int_file %d fpid %d\n", i, pid_in_file[i], fpid);
				if (get_next_signal(infs[i], file_type[i], curr, fpid, fstats[i], missing_dict_vals) == -1) {
					MERR("create_indexes : get_next_signal failed for file %d/%d\n", i, n_files_opened);
					return -1;
				}
				pid_in_file[i] = fpid; // current pid after the one we wanted
			}
			//MLOG("i=%d c_pid=%d fpid=%d curr %d %d %d\n",i,c_pid,fpid,curr.pid,n_files_opened,n_open_in_files);
		}

		// write data to output files
		if (curr.pid >= 0) {
			if (write_indexes(curr) < 0) {
				//MERR("MedConvert: create_indexes: curr packet for pid %d was not written...\n", curr.pid);

			}
			else
				all_pids.push_back(curr.pid);
		}

		n_pids_extracted++;
		//if (n_pids_extracted % 100000 == 0)
		//	MLOG("MedConvert: create_indexes: extracted %d pids (%d open files)\n", n_pids_extracted, n_open_in_files);
		if (n_pids_extracted % 10000 == 0 && (int)difftime(time(NULL), last_time_print) >= 30) {
			last_time_print = time(NULL);
			float time_elapsed = (float)difftime(time(NULL), start);
			float estimate_time = float(n_open_in_files) / (n_files_opened - n_open_in_files) * time_elapsed;
			MLOG("Processed %d out of %d(%2.2f%) time elapsed: %2.1f Minutes, estimate time to finish %2.1f Minutes."
				" extracted %d pids\n",
				n_files_opened - n_open_in_files, n_files_opened,
				100.0*((n_files_opened - n_open_in_files) / float(n_files_opened)),
				time_elapsed / 60, estimate_time / 60.0, n_pids_extracted);
		}
	}
	bool too_many_missing = false;
	for (auto& entry : missing_dict_vals) {
		MWARN("MedConvert: saw missing entry [%s]:[%s] %d times\n", entry.first.first.c_str(), entry.first.second.c_str(), entry.second);
		if (safe_mode && entry.second > 50) {
			MERR("%d > 50 missing entries is too much... refusing to create repo!\n", entry.second);
			too_many_missing = true;
		}
	}
	if (too_many_missing)
		MTHROW_AND_ERR("too many missing values...\n");
	for (auto& entry : missing_forced_signals) {
		MWARN("MedConvert: saw missing_forced_signal [%s] %d times\n", entry.first.c_str(), entry.second);
		if (safe_mode && 1.0*entry.second / n_pids_extracted > 0.05) {
			MERR("%d / %d missing_forced_signal is too much... refusing to create repo!\n", entry.second, n_pids_extracted);
			too_many_missing = true;
		}
	}
	// all files are closed, all are written correctly

	// print statistics for data files
	MLOG("Statistics for %d data files\n", fstats.size());
	for (auto& stat : fstats) {
		float ratio = (float)(stat.n_parsed_lines + 1) / (float)(stat.n_relevant_lines + 1);
		float bad_ratio = (float)(stat.n_bad_format_lines + 1) / (float)(stat.n_relevant_lines + 1);
		MLOG("file [%d] : %s : n_lines %d , n_relevant_lines %d , n_bad_format_lines %d n_loaded_lines %d : loaded %g\n",
			stat.id, stat.fname.c_str(), stat.n_lines, stat.n_relevant_lines, stat.n_bad_format_lines, stat.n_parsed_lines,
			ratio);
		if (ratio < 0.01 || bad_ratio > 0.05) {
			if (stat.n_relevant_lines > 1000) {
				MTHROW_AND_ERR("%d/%d lines loaded for file [%s]\n", stat.n_parsed_lines, stat.n_relevant_lines, stat.fname.c_str());
			}
			else MWARN("%d/%d lines loaded for file [%s]\n", stat.n_parsed_lines, stat.n_relevant_lines, stat.fname.c_str());
		}

	}

	MLOG("Finished reading all pids (%d pids extracted) - closing index and data files\n", n_pids_extracted);
	if (mode < 3)
		close_indexes();
	else {
		write_all_indexes(all_pids);
	}

	return 0;
}

//------------------------------------------------
int MedConvert::open_indexes()
{
	int i;

	if (mode < 3) {
		index_f.resize(index_fnames.size());
		fill(index_f.begin(), index_f.end(), (ofstream *)NULL);
		unsigned long long magic = MED_MAGIC_NUM;
		int index_mode = 0;
		for (i = 0; i < index_fnames.size(); i++)
			if (sids_to_load[serial2sid[i]]) {
				index_f[i] = (ofstream *)new ofstream;
				index_f[i]->open(index_fnames[i], ios::out | ios::binary);
				if (!index_f[i]->is_open()) {
					MERR("MedConvert:: open_indexes:: can't open output file %s\n", index_fnames[i].c_str());
					return -1;
				}
				// writing index file header (current mode is 0)
				index_f[i]->write((char *)&magic, sizeof(unsigned long long));
				index_f[i]->write((char *)&index_mode, sizeof(int));
				index_f[i]->flush();
			}
	}
	else {
		indexes.resize(index_fnames.size());
		for (i = 0; i < indexes.size(); i++) {
			indexes[i].base = 4; // 4 bytes at start are for format version of data file
			indexes[i].sid = serial2sid[i];
			indexes[i].factor = sigs.Sid2Info[serial2sid[i]].bytes_len;
			indexes[i].last_len = 0;
			indexes[i].work_area = NULL;
		}
	}

	data_f.resize(data_fnames.size());
	fill(data_f.begin(), data_f.end(), (ofstream *)NULL);
	data_f_pos.resize(data_fnames.size());
	for (i = 0; i < data_fnames.size(); i++)
		if (sids_to_load[serial2sid[i]]) {
			data_f[i] = (ofstream *)new ofstream;
			data_f[i]->open(data_fnames[i], ios::out | ios::binary);
			if (!data_f[i]->is_open()) {
				MERR("MedConvert:: open_indexes:: can't open output file %s\n", data_fnames[i].c_str());
				return -1;
			}
			MLOG("data_f file %d %s opened\n", i, data_fnames[i].c_str());
			// writing repository stripped format bits to data fo;es
			int data_format = REPOSITORY_STRIPPED_FORMAT;
			data_f[i]->write((char *)&data_format, sizeof(int));
			data_f_pos[i] = sizeof(int);
			data_f[i]->flush();
		}

	MLOG("opened %d index files and %d data files\n", index_fnames.size(), data_fnames.size());
	return 0;
}

//------------------------------------------------
int MedConvert::close_indexes()
{
	int i;

	for (i = 0; i < index_f.size(); i++) {
		if (index_f[i] != NULL) {
			index_f[i]->close();
			delete index_f[i];
			index_f[i] = NULL;
		}
	}
	for (i = 0; i < data_f.size(); i++) {
		if (data_f[i] != NULL) {
			data_f[i]->close();
			delete data_f[i];
			data_f[i] = NULL;
		}
	}
	return 0;

}

//------------------------------------------------
int MedConvert::write_all_indexes(vector<int> &all_pids)
{
	for (int i = 0; i < indexes.size(); i++) {
		if (sids_to_load[serial2sid[i]]) {
			if (indexes[i].write_to_file(index_fnames[i]) < 0)
				return -1;
		}
	}
	for (int i = 0; i < data_f.size(); i++) {
		if (data_f[i] != NULL) {
			data_f[i]->close();
			delete data_f[i];
			data_f[i] = NULL;
		}
	}

	if (load_only.size() == 0) {
		// writing all_pids to a file - a list of all available pids significantly speeds up repository usage in some cases
		// format - number of all, followed by the pids
		// in case of secondary loads - we do not update the all_pids file.
		all_pids[0] = (int)all_pids.size() - 1;
		string fname_pids = out_path + "/" + rep_files_prefix + "_all_pids.list";

		if (write_bin_file_IM(fname_pids, (unsigned char *)&all_pids[0], sizeof(int)*all_pids.size()) < 0) {
			MERR("ERROR: Could not write file %s for %d pids\n", fname_pids.c_str(), all_pids[0]);
			return -1;
		}
	}

	return 0;
}

//------------------------------------------------
int MedConvert::write_indexes(pid_data &curr)
{
	if (curr.pid < 0)
		MTHROW_AND_ERR("MedConvert::write_indexes negative pid %d", curr.pid);
	// first we sort all elements by time
	int i;
	for (i = 0; i < curr.raw_data.size(); i++) {
		sort(curr.raw_data[i].begin(), curr.raw_data[i].end(),
			[](const collected_data &v1, const collected_data &v2)
		{
			if (v1.date == v2.date) {
				if (v1.date2 == v2.date2) {
					if (v1.time == v2.time) {
						if (v1.time2 == v2.time2)
							return (v1.val < v2.val) || (v1.longVal < v2.longVal);
						else
							return v1.time2 < v2.time2;
					}
					else
						return v1.time < v2.time;
				}
				else
					return v1.date2 < v2.date2;
			}
			else
				return v1.date < v2.date;
		});
	}

	// getting rid of duplicates
	vector<collected_data>::iterator it;
	for (i = 0; i < curr.raw_data.size(); i++) {
		it = unique(curr.raw_data[i].begin(), curr.raw_data[i].end(), [](const collected_data &v1, const collected_data &v2) {return ((v1.longVal == v2.longVal) && (v1.val == v2.val) && (v1.date == v2.date) && (v1.date2 == v2.date2) && (v1.time == v2.time) && (v1.time2 == v2.time2) && (v1.val1 == v2.val1) && (v1.val2 == v2.val2)); });
		curr.raw_data[i].resize(distance(curr.raw_data[i].begin(), it));
	}


	// sanity checks - things we force to have, and things we force to have as single.
//	if (curr.raw_data[sid2serial[dict.id(string("GENDER"))]].size() != 1) return -1;
//	if (curr.raw_data[sid2serial[dict.id(string("BYEAR"))]].size() != 1) return -1;

	// forced signals
	for (i = 0; i < forced.size(); i++) {
		if (curr.raw_data[sid2serial[dict.id(forced[i])]].size() != 1) {
			if (missing_forced_signals.find(forced[i]) == missing_forced_signals.end())
				missing_forced_signals[forced[i]] = 1;
			else
				missing_forced_signals[forced[i]] += 1;
			if (missing_forced_signals[forced[i]] < 10)
				MLOG("MedConvert: pid %d is missing forced signal %s (%d,%d,%d)\n", curr.pid, forced[i].c_str(),
					dict.id(forced[i]), sid2serial[dict.id(forced[i])], curr.raw_data[sid2serial[dict.id(forced[i])]].size());
			return -1;
		}
	}

	// writing indexes
	int fno;
	int n_pid_sigs;
	for (fno = 0; fno < index_fnames.size(); fno++)
		if (data_f[fno] != NULL)
		{
			n_pid_sigs = 0;
			if (mode < 3) {
				for (i = 0; i < curr.raw_data.size(); i++)
					if (curr.raw_data[i].size() > 0 && serial2siginfo[i].fno == fno &&
						(serial2siginfo[i].type >= 0 && serial2siginfo[i].type < T_Last))
						//if (curr.raw_data[i].size() > 0 && sid2fno[serial2sid[i]] == fno &&
						//(sigs.type(serial2sid[i]) >= 0 &&   sigs.type(serial2sid[i])<T_Last))
						n_pid_sigs++;
			}
			else {
				// in this mode fno is i... and there's one option for it
				if (curr.raw_data[fno].size() > 0 && serial2siginfo[fno].type >= 0 && serial2siginfo[fno].type < T_Last)
					n_pid_sigs++;
				//MLOG("i=%d/%d n_pid_sigs %d\n", fno, index_fnames.size(), n_pid_sigs);
			}


			if (n_pid_sigs > 0) {

				// write packet header: magic number + pid + number of signals
				unsigned long long magic = MED_MAGIC_NUM;
				int pid = curr.pid;
				if (mode < 3) {
					index_f[fno]->write((char *)&magic, sizeof(unsigned long long));
					index_f[fno]->write((char *)&pid, sizeof(int));
					index_f[fno]->write((char *)&n_pid_sigs, sizeof(int));
				}

				// write data and index pointer for each signal
				for (i = 0; i < curr.raw_data.size(); i++) {

					int ilen = (int)curr.raw_data[i].size();
					if (ilen > 0) {
						int sid = serial2sid[i];
						int sid_type = serial2siginfo[i].type; //sigs.type(sid);
						int sid_fno = serial2siginfo[i].fno; //sid2fno[sid];

						if ((ilen>0) && (sid_fno == fno) && (sid_type >= 0 && sid_type < T_Last)) {

							//int sid = serial2sid[i];
							unsigned short file_n = fno;
							unsigned long long pos = data_f_pos[fno];
							int len = 0;

							if (sid_type == T_Value) {
								len = (int)sizeof(SVal)*ilen;
								SVal sv;
								for (int j = 0; j < ilen; j++) {
									sv.val = curr.raw_data[i][j].val;
									data_f[fno]->write((char *)&sv, sizeof(SVal));
								}
							}

							if (sid_type == T_DateVal) {
								len = (int)sizeof(SDateVal)*ilen;
								SDateVal sdv;
								for (int j = 0; j < ilen; j++) {
									sdv.date = curr.raw_data[i][j].date;
									sdv.val = curr.raw_data[i][j].val;
									data_f[fno]->write((char *)&sdv, sizeof(SDateVal));
								}
							}

							if (sid_type == T_DateRangeVal) {
								len = (int)sizeof(SDateRangeVal)*ilen;
								SDateRangeVal sdrv;
								for (int j = 0; j < ilen; j++) {
									sdrv.date_start = curr.raw_data[i][j].date;
									sdrv.date_end = curr.raw_data[i][j].date2;
									sdrv.val = curr.raw_data[i][j].val;
									data_f[fno]->write((char *)&sdrv, sizeof(SDateRangeVal));
								}
							}

							if (sid_type == T_TimeVal) {
								len = (int)sizeof(STimeVal)*ilen;
								STimeVal stv;
								for (int j = 0; j < ilen; j++) {
									stv.time = curr.raw_data[i][j].time;
									stv.val = curr.raw_data[i][j].val;
									data_f[fno]->write((char *)&stv, sizeof(STimeVal));
								}
							}

							if (sid_type == T_TimeRangeVal) {
								len = (int)sizeof(STimeRangeVal)*ilen;
								STimeRangeVal strv;
								for (int j = 0; j < ilen; j++) {
									strv.time_start = curr.raw_data[i][j].time;
									strv.time_end = curr.raw_data[i][j].time2;
									strv.val = curr.raw_data[i][j].val;
									data_f[fno]->write((char *)&strv, sizeof(STimeRangeVal));
								}
							}

							if (sid_type == T_TimeStamp) {
								len = (int)sizeof(STimeStamp)*ilen;
								STimeStamp sts;
								for (int j = 0; j < ilen; j++) {
									sts.time = curr.raw_data[i][j].time;
									data_f[fno]->write((char *)&sts, sizeof(STimeStamp));
								}
							}

							if (sid_type == T_DateVal2) {
								len = (int)sizeof(SDateVal2)*ilen;
								SDateVal2 sdv2;
								for (int j = 0; j < ilen; j++) {
									sdv2.date = curr.raw_data[i][j].date;
									sdv2.val = curr.raw_data[i][j].val;
									sdv2.val2 = curr.raw_data[i][j].val2;
									data_f[fno]->write((char *)&sdv2, sizeof(SDateVal2));
								}
							}

							if (sid_type == T_TimeLongVal) {
								len = (int)sizeof(STimeLongVal)*ilen;
								STimeLongVal stv;
								for (int j = 0; j < ilen; j++) {
									stv.time = curr.raw_data[i][j].time;
									stv.val = curr.raw_data[i][j].longVal;
									data_f[fno]->write((char *)&stv, sizeof(STimeLongVal));
								}
							}

							if (sid_type == T_DateShort2) {
								len = (int)sizeof(SDateShort2)*ilen;
								SDateShort2 sds2;
								for (int j = 0; j < ilen; j++) {
									sds2.date = curr.raw_data[i][j].date;
									sds2.val1 = curr.raw_data[i][j].val1;
									sds2.val2 = curr.raw_data[i][j].val2;
									data_f[fno]->write((char *)&sds2, sizeof(SDateShort2));
								}
							}

							if (sid_type == T_ValShort2) {
								len = (int)sizeof(SValShort2)*ilen;
								SValShort2 svs2;
								for (int j = 0; j < ilen; j++) {
									svs2.val1 = curr.raw_data[i][j].val1;
									svs2.val2 = curr.raw_data[i][j].val2;
									data_f[fno]->write((char *)&svs2, sizeof(SValShort2));
								}
							}

							if (sid_type == T_ValShort4) {
								len = (int)sizeof(SValShort4)*ilen;
								SValShort4 svs4;
								for (int j = 0; j < ilen; j++) {
									svs4.val1 = curr.raw_data[i][j].val1;
									svs4.val2 = curr.raw_data[i][j].val2;
									svs4.val3 = curr.raw_data[i][j].val3;
									svs4.val4 = curr.raw_data[i][j].val4;
									data_f[fno]->write((char *)&svs4, sizeof(SValShort4));
								}
							}

							if (sid_type == T_CompactDateVal) {
								len = (int)sizeof(SCompactDateVal)*ilen;
								SCompactDateVal scdv;
								for (int j = 0; j < ilen; j++) {
									scdv.compact_date = date_to_compact_date(curr.raw_data[i][j].date);
									scdv.val = curr.raw_data[i][j].val1;
									data_f[fno]->write((char *)&scdv, sizeof(SCompactDateVal));
								}
							}

							//MLOG("writing to fno %d : sid %d file_n %d pos %ld len %d\n", fno, sid, file_n, pos, len);
							if (mode < 3) {
								index_f[fno]->write((char *)&sid, sizeof(int));
								index_f[fno]->write((char *)&file_n, sizeof(short));
								index_f[fno]->write((char *)&pos, sizeof(unsigned long long));
								index_f[fno]->write((char *)&len, sizeof(int));
							}
							else {
								indexes[fno].insert(pid, ilen);
							}
							data_f_pos[fno] += len;

						}
					}
				}
			}

		}
	return 0;
}

/////////////////////////////////////////////////////////////////
// mode >=2 related
//-----------------------------------------------------------------------------------------------------------------
int MedConvert::generate_prefix_names()
{
	prefix_names.clear();
	sid2fno.clear();
	serial2sid.clear();
	sid2serial.clear();
	serial2siginfo.clear();
	for (int i = 0; i < sigs.signals_names.size(); i++) {
		string fixed_sig_name = sigs.signals_names[i];
		replace_all(fixed_sig_name, "/", "_div_");
		replace_all(fixed_sig_name, ":", "_over_");
		replace_all(fixed_sig_name, "%", "_percent_");
		string name = rep_files_prefix + "_" + fixed_sig_name; //sigs.signals_names[i];
		int sid = sigs.sid(sigs.signals_names[i]); // signals_ids[i];
		prefix_names.push_back(name);
		sid2fno[sid] = i;
		serial2sid.push_back(sid);
		sig_info si;
		si.fno = i;
		si.serial = i;
		si.type = sigs.type(sid);
		si.sid = sid;
		serial2siginfo.push_back(si);
		sid2serial[sid] = i;
	}

	return 0;
}
