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
#include <boost/filesystem.hpp>

using namespace boost;

#define LOCAL_SECTION LOG_CONVERT
#define LOCAL_LEVEL	LOG_DEF_LEVEL
extern MedLogger global_logger;

//#define MAX_PID_TO_TAKE	1000
#define MAX_PID_TO_TAKE	1000000000
//#define MAX_PID_TO_TAKE	5010000

vector<mutex> sid_mutex(MAX_SID_NUMBER);


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
	default_time_unit = MedTime::Date;
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
				if (fields[0].compare("MODE") == 0) mode = med_stoi(fields[1]);
				if (fields[0].compare("SAFE_MODE") == 0) safe_mode = med_stoi(fields[1]);
				if (fields[0].compare("PREFIX") == 0) rep_files_prefix = fields[1];
				if (fields[0].compare("RELATIVE") == 0) relative = 1;
				if (fields[0].compare("TIMEUNIT") == 0 || fields[0].compare("TIME_UNIT") == 0) {
					default_time_unit = med_stoi(fields[1]);
					MLOG("MedConvert: Will convert all dates field to MedTime::[%d] format\n", default_time_unit);
				}
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
				int fno = med_stoi(fields[0]);
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
	MLOG("MedConvert: read_signal_to_files: %s\n", fname.c_str());
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
				int fno = med_stoi(fields[0]);
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

//fetch file name
void get_rel_filename(vector<string> &paths) {
	for (size_t i = 0; i < paths.size(); i++)
	{
		boost::filesystem::path pt(paths[i]);
		paths[i] = pt.filename().string();
	}
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

		MLOG("MedConvert: read_all: read signal_to_files file [%s]\n", signal_to_files_fname.c_str());

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

			if (verbose_open_files)
				MLOG("i=%d index %s data %s\n", i, index_fnames[i].c_str(), data_fnames[i].c_str());
		}
	}
	if (add_path_to_name_IM(out_path, repository_config_fname) == -1 ||
		add_path_to_name_IM(out_path, sig_fnames) == -1)
		return -1;
	if (create_signals_config() < 0)
		MTHROW_AND_ERR("MedConvert: read_all(): failed generating signals config file\n");
	get_rel_filename(sig_fnames);

	// Create repository config file
	if (create_repository_config() < 0)
		MTHROW_AND_ERR("MedConvert: read_all(): failed generating repository config file\n");

	// Copy dict files as-is to output directory
	//validate dir exits:
	for (string d : dict_fnames)
	{
		add_path_to_name_IM(out_path, d);
		boost::filesystem::path pt_dict(d);
		boost::filesystem::create_directories(pt_dict.parent_path().string());
	}
	if (copy_files_IM(path, out_path, dict_fnames) < 0)
		MTHROW_AND_ERR("MedConvert : read_all() : failed copying files from in to out directory\n");

	// add path to more files + fix paths
	if (add_path_to_name_IM(path, dict_fnames) == -1 ||
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
bool read_file_to_buffer(ifstream &inf, vector<string> &buffered_lines, int read_lines_buffer) {
	if (!inf.is_open())
		return true; // file is closed (and no lines in buffer) nothing to do
	std::ios::sync_with_stdio(false);
	std::ios_base::sync_with_stdio(false);


	string line;
	int curr_ln = 0;
	//no more than 1 thread reading files
#pragma omp critical
	while ((read_lines_buffer <= 0 || curr_ln < read_lines_buffer) && getline(inf, line)) {
		if (line[line.size() - 1] == '\r')
			line = line.substr(0, line.size() - 1);
		boost::trim(line);
		if (line.empty() || (line[0] == '#'))
			continue;
		buffered_lines.push_back(line);
		++curr_ln;
	}
	if (inf.eof())
		inf.close();
	if (buffered_lines.size() == 0)
		return true;
	return false;
}

//-------------------------------------------------------------------------------------------------------------------------------------------
void MedConvert::collect_lines(vector<string> &lines, vector<int> &f_i, int file_i, vector<string> &buffered_lines, int &buffer_pos, ifstream &inf, int file_type, pid_data &curr, int &fpid, file_stat& curr_fstat, map<pair<string, string>, int>&)
{
	bool get_next = true;
	collected_data cd;
	GenericSigVec cd_sv;
	cd_sv.set_data(&cd.buf[0], 1);
	string vfield, vfield2;
	while (get_next) {
		if (buffer_pos >= buffered_lines.size() || buffered_lines.empty()) {
			//need to refresh and read more/again from the buffer:
			buffer_pos = 0;
			buffered_lines.clear();
			if (read_lines_buffer > 0)
				buffered_lines.reserve(read_lines_buffer);
			bool finished = read_file_to_buffer(inf, buffered_lines, read_lines_buffer);
			if (finished) {
				get_next = false;
				fpid = -2; // signing the true end of this file
				MLOG("MedConvert: Closing file %s - reached end of file (n_lines %d) n_open_in_files is %d\n", curr_fstat.fname.c_str(), curr_fstat.n_lines, n_open_in_files);
				--n_open_in_files;
				return;
			}
		}

		string &curr_line = buffered_lines[buffer_pos++];
		curr_fstat.n_lines++;
		curr_fstat.n_relevant_lines++;
		int line_pid = -1;
		try {
			line_pid = med_stoi(curr_line);
		}
		catch (...) {
			MERR("ERROR: bad format in file %s with first token of pid, in line %d:\n%s\n",
				curr_fstat.fname.c_str(), curr_fstat.n_parsed_lines, curr_line.c_str());
			if (!full_error_file.empty() && (err_log_file.good()))
				err_log_file << "ERROR: bad format in file " << curr_fstat.fname << " with first token of pid, in line " << curr_fstat.n_parsed_lines << " line: " << curr_line << "\n";
			continue;
		}

		if (line_pid == curr.pid) {
#pragma omp critical
			{
				//MLOG("pid is %d : file %d : pushing line : %s\n", line_pid, file_i, curr_line.c_str());
				lines.push_back(curr_line);
				f_i.push_back(file_i);
			}
		}
		else if (line_pid < fpid) {
			MWARN("MedConvert: get_next_signal: fpid is %d , but got line: %s\n", fpid, curr_line.c_str());
			if (!full_error_file.empty() && (err_log_file.good()))
				err_log_file << "MedConvert: get_next_signal: fpid is " << fpid << " , but got line: " << curr_line << "\n";

			if (safe_mode)
				MERR("MedConvert: ERROR: file %s seems to be not sorted by pid\n", curr_fstat.fname.c_str());
		}
		else {
			fpid = line_pid;
			--buffer_pos; // roll file back one line
			--curr_fstat.n_lines;
			--curr_fstat.n_relevant_lines;
			get_next = false;
		}
	}

	if (fpid > MAX_PID_TO_TAKE) {
		fpid = -1;
		--n_open_in_files;

	}
}


//-------------------------------------------------------------------------------------------------------------------------------------------------
void MedConvert::parse_fields_into_gsv(string &curr_line, vector<string> &fields, int sid, GenericSigVec &cd_sv)
{
	char convert_mode = safe_mode > 0 ? 2 : 1;
	int section = dict.section_id(sigs.name(sid));
	SignalInfo& info = sigs.Sid2Info[sid];
	cd_sv.init(info);
	if (cd_sv.size() > MAX_COLLECTED_DATA_SIZE) {
		if (!full_error_file.empty() && (err_log_file.good()))
			err_log_file << "ERROR: cd_sv.size() (" << (int)cd_sv.size() << ") > MAX_COLLECTED_DATA_SIZE (" << (int)MAX_COLLECTED_DATA_SIZE << "), Please Increase MAX_COLLECTED_DATA_SIZE\n";
		MTHROW_AND_ERR("ERROR: cd_sv.size() (%d) > MAX_COLLECTED_DATA_SIZE (%d), Please Increase MAX_COLLECTED_DATA_SIZE\n", (int)cd_sv.size(), (int)MAX_COLLECTED_DATA_SIZE);
	}

	if ((cd_sv.n_time + cd_sv.n_val) > fields.size() - 2) {
		if (!full_error_file.empty() && (err_log_file.good()))
			err_log_file << "ERROR: in signal " << fields[1] << " expecting " << cd_sv.n_time << " time channels, and " << cd_sv.n_val << " value channels, but only " << fields.size() << " fields in line " << curr_line << "\n";
		MTHROW_AND_ERR("ERROR: in signal %s expecting %d time channels, and %d value channels, but only %d fields in line %s\n", fields[1].c_str(), cd_sv.n_time, cd_sv.n_val, (int)fields.size(), curr_line.c_str());

	}

	int time_unit = info.time_unit == MedTime::Undefined ? default_time_unit : info.time_unit;

	int field_i = 2;

	for (int tchan = 0; tchan < cd_sv.n_time; tchan++) {
		switch (cd_sv.time_channel_types[tchan]) {

		case GenericSigVec::type_enc::UINT8:   //unsigned char
		case GenericSigVec::type_enc::UINT32:  //unsigned int
		case GenericSigVec::type_enc::UINT64:  //unsigned long long
		case GenericSigVec::type_enc::INT8:    //char
		case GenericSigVec::type_enc::INT16:   //short
		case GenericSigVec::type_enc::FLOAT32: //float
		case GenericSigVec::type_enc::FLOAT64: //double
		case GenericSigVec::type_enc::FLOAT80: //long double
		case GenericSigVec::type_enc::INT32:   //int
			cd_sv.setTime(0, tchan, stoi(fields[field_i]));
			//cd_sv.setTime(0, tchan, med_time_converter.convert_datetime_safe(time_unit, fields[field_i], convert_mode)); break;

		case GenericSigVec::type_enc::INT64:   //long long
			cd_sv.setTime<long long>(0, tchan, med_time_converter.convert_datetime_safe(time_unit, fields[field_i], convert_mode));	break;

		case GenericSigVec::type_enc::UINT16:  //unsigned short
			int value = (int)med_time_converter.convert_datetime_safe(time_unit, fields[field_i], convert_mode);
			if (value < 0) {
				if (!full_error_file.empty() && (err_log_file.good()))
					err_log_file << "MedConvert: get_next_signal: Detected attempt to assign negative number (" << value << ") into unsigned time channel " << tchan << " :: curr_line is :" << curr_line << "\n";
				MTHROW_AND_ERR("MedConvert: get_next_signal: Detected attempt to assign negative number (%d) into unsigned time channel %d :: curr_line is '%s'\n", value, tchan, curr_line.c_str());
			}
			cd_sv.setTime<unsigned short>(0, tchan, value);
			break;
		}
		field_i++;
	}


	for (int vchan = 0; vchan < cd_sv.n_val; vchan++) {
		switch (cd_sv.val_channel_types[vchan]) {

		case GenericSigVec::type_enc::UINT8:   //unsigned char
		case GenericSigVec::type_enc::UINT32:  //unsigned int
		case GenericSigVec::type_enc::INT8:    //char
		case GenericSigVec::type_enc::INT32:   //int
		case GenericSigVec::type_enc::INT64:   //long long
		case GenericSigVec::type_enc::FLOAT64: //double
		case GenericSigVec::type_enc::FLOAT80: //long double
		case GenericSigVec::type_enc::FLOAT32: //float
			if (sigs.is_categorical_channel(sid, vchan))
				cd_sv.setVal(0, vchan, dict.get_id_or_throw(section, fields[field_i]));
			else
				cd_sv.setVal(0, vchan, med_stof(fields[field_i]));
			break;

		case GenericSigVec::type_enc::UINT16:  //unsigned short
			if (sigs.is_categorical_channel(sid, vchan))
				cd_sv.setVal<unsigned short>(0, vchan, dict.get_id_or_throw(section, fields[field_i]));
			else {
				auto value = med_stoi(fields[field_i]);
				if (value < 0) {
					if (!full_error_file.empty() && (err_log_file.good()))
						err_log_file << "MedConvert: get_next_signal: Detected attempt to assign negative number (" << value << ") into unsigned value channel " << vchan << "  :: curr_line is " << curr_line << "\n";
					MTHROW_AND_ERR("MedConvert: get_next_signal: Detected attempt to assign negative number (%d) into unsigned value channel %d :: curr_line is '%s'\n", value, vchan, curr_line.c_str());
				}
				cd_sv.setVal<unsigned short>(0, vchan, (unsigned short)value);
			}
			break;

		case GenericSigVec::type_enc::UINT64:  //unsigned long long
			if (sigs.is_categorical_channel(sid, vchan))
				cd_sv.setVal<long long>(0, vchan, (long long)dict.get_id_or_throw(section, fields[field_i]));
			else cd_sv.setVal<long long>(0, vchan, med_stof(fields[field_i]));
			break;

		case GenericSigVec::type_enc::INT16:   //short
			if (sigs.is_categorical_channel(sid, vchan))
				cd_sv.setVal<short>(0, vchan, dict.get_id_or_throw(section, fields[field_i]));
			else cd_sv.setVal<short>(0, vchan, med_stof(fields[field_i]));
			break;
		}
		field_i++;
	}
}


//-------------------------------------------------------------------------------------------------------------------------------------------------
void MedConvert::get_next_signal_all_lines(vector<string> &lines, vector<int> &f_i, pid_data &curr, vector<file_stat> &fstat, map<pair<string, string>, int>& missing_dict_vals)
{
	//MLOG("===> lines %d\n", lines.size());
//#pragma omp parallel for schedule(dynamic) if (run_parallel)

	vector<collected_data> cds(lines.size());
#pragma omp parallel for if (run_parallel)
	for (int k = 0; k < lines.size(); k++) {
		//MLOG("k=%d line %s\n", k, lines[k].c_str());
		collected_data &cd = cds[k];
		GenericSigVec cd_sv;
		cd_sv.set_data(&cd.buf[0], 1);
		string vfield, vfield2;
		int sid;

		string &curr_line = lines[k];
		file_stat &curr_fstat = fstat[f_i[k]];
		vector<string> fields;
		split(fields, curr_line, boost::is_any_of("\t"));
		//split(fields, curr_line, boost::is_from_range('\t','\t'));

		cd.zero();

		if (fields.size() < 3) {
			MERR("MedConvert: ERROR: Too few fields in file %s, line : %s\n", curr_fstat.fname.c_str(), curr_line.c_str());
			if (!full_error_file.empty() && (err_log_file.good()))
				err_log_file << "TOO_FEW_FIELDS_IN_LINE" << "\tfile: " << curr_fstat.fname << "\tline: " << curr_line << "\n";
			continue;
		}

		if (codes2names.find(fields[1]) == codes2names.end()) {
			MERR("MedConvert: ERROR: unrecognized signal name %s (need to add to codes_to_signals file) in file %s :: curr_line is %s\n",
				fields[1].c_str(), curr_fstat.fname.c_str(), curr_line.c_str());
			if (!full_error_file.empty() && (err_log_file.good()))
				err_log_file << "UNRECOGNIZED_SIGNAL" << "\t(" << fields[1] << ")\tfile: " << curr_fstat.fname << "\tline: " << curr_line << "\n";
			continue;
		}

		if ((sid = sigs.sid(codes2names[fields[1]])) < 0) {
			MERR("MedConvert: ERROR: signal name %s converted to %s is not in dict in file %s :: curr_line is %s\n",
				fields[1].c_str(), codes2names[fields[1]].c_str(), curr_fstat.fname.c_str(), curr_line.c_str());
			if (!full_error_file.empty() && (err_log_file.good()))
				err_log_file << "NO_SID_FOR_SIGNAL" << "\t(" << fields[1] << ")\tfile: " << curr_fstat.fname << "\tline: " << curr_line << "\n";
			continue;
		}
		//MLOG("sig %s %s sid %d to_load %d sigs.name %d %d\n", fields[1].c_str(), codes2names[fields[1]].c_str(), sid, sids_to_load[sid], (int)sigs.Name2Sid.size(), sigs.Name2Sid[codes2names[fields[1]]]);
		if (!sids_to_load[sid])
			continue;

		try {
			int i = sid2serial[sid];
			parse_fields_into_gsv(curr_line, fields, sid, cd_sv);
			cd.serial = i;
			/*
#pragma omp critical
			{
				//lock_guard<mutex> guard(sid_mutex[i]);
				//MLOG("inserted\n", curr_line.c_str(), i, sid, sigs.name(sid).c_str(), f_i[k]);
				curr.raw_data[i].push_back(cd);
				curr_fstat.n_parsed_lines++;
			}
			*/
		} catch (invalid_argument e) {

			pair<string, string> my_key = make_pair(sigs.name(sid), string(e.what()));
			if (missing_dict_vals.find(my_key) == missing_dict_vals.end()) {
#pragma omp critical
				missing_dict_vals[my_key] = 1;
			} else {
#pragma omp atomic
				++missing_dict_vals[my_key];
			}
			if (missing_dict_vals.size() < 10)
				MWARN("MedConvert::get_next_signal: missing from dictionary (sig [%s], type %d) : file [%s] : line [%s] \n",
					sigs.name(sid).c_str(), sigs.type(sid), curr_fstat.fname.c_str(), curr_line.c_str());
			if (!full_error_file.empty()) {
				if (err_log_file.good()) {
					err_log_file << "MISSING_FROM_DICTIONARY"
						<< "\t" << sigs.name(sid) << " (" << sigs.type(sid) << ")"
						<< "\t" << curr_fstat.fname
						<< "\t" << "\"" << curr_line << "\""
						<< "\n";
				}
			}
		} catch (...) {

			curr_fstat.n_bad_format_lines++;
			if (curr_fstat.n_bad_format_lines < 10) {
				MWARN("MedConvert::get_next_signal: bad format in parsing file %s in line %d:\n%s\n",
					curr_fstat.fname.c_str(), curr_fstat.n_parsed_lines, curr_line.c_str());
			}
			if (!full_error_file.empty()) {
				if (err_log_file.good()) {
					err_log_file << "BAD_FILE_FORMAT"
						<< "\t" << curr_fstat.fname //<< " (" << file_type << ")"
						<< "\t" << curr_fstat.n_parsed_lines
						<< "\t" << "\"" << curr_line << "\""
						<< "\n";
				}
			}
		}
	}

	for (int k = 0; k < cds.size(); k++) {
		file_stat &curr_fstat = fstat[f_i[k]];
		if (cds[k].serial >= 0) {
			{
				curr.raw_data[cds[k].serial].push_back(cds[k]);
				curr_fstat.n_parsed_lines++;
			}
		}
	}
}

//------------------------------------------------
int MedConvert::create_signals_config()
{
	assert(sig_fnames.size() == 1);
	MLOG("MedConvert::create_signals_config [%s]\n", sig_fnames[0].c_str());
	signals_config_f.open(sig_fnames[0].c_str(), ios::out);
	for (unsigned int i = 0; i < sig_fnames.size(); i++) {
		for (SignalInfo& info : this->sigs.Sid2Info) {
			if (info.sid < 0)
				continue;
			signals_config_f << "SIGNAL\t" << info.name << "\t" << info.sid << "\t" << info.type << "\t" << info.description << "\t";
			for (int j = 0; j < info.n_val_channels; j++)
				signals_config_f << info.is_categorical_per_val_channel[j];
			signals_config_f << "\t";
			for (int j = 0; j < info.n_val_channels; j++) {
				signals_config_f << info.unit_of_measurement_per_val_channel[j];
				if (j < info.n_val_channels - 1)
					signals_config_f << '|';
			}
			if (info.time_unit != MedTime::Undefined)
				signals_config_f << "\t" << info.time_unit;
			signals_config_f << endl;
		}
	}
	signals_config_f.close();
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

	// support only a single signals file
	assert(sig_fnames.size() == 1);
	repository_config_f << "SIGNAL\t" << sig_fnames[0].c_str() << endl;

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
	repository_config_f << "TIMEUNIT\t" << default_time_unit << endl;
	repository_config_f.close();
	return 0;
}

//------------------------------------------------
int MedConvert::create_indexes()
{
	std::ios::sync_with_stdio(false);
	std::ios_base::sync_with_stdio(false);

	pid_data curr;
	vector<pid_data> curr_f; // for each file

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
	for (int i = 0; i < in_data_fnames.size(); i++) {
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
		if (verbose_open_files)
			MLOG("MedConvert: opened file %s for input file (%d) , of type %d\n", in_data_fnames[i].c_str(), n_open_in_files, file_type[n_open_in_files]);
		n_open_in_files++;
	}

	for (int i = 0; i < in_strings_data_fnames.size(); i++) {
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
		if (verbose_open_files)
			MLOG("MedConvert: opened file %s for input file (%d) , of type %d\n", in_strings_data_fnames[i].c_str(), n_open_in_files, file_type[n_open_in_files]);
		n_open_in_files++;
	}

	int c_pid = -1;
	int n_files_opened = n_open_in_files;
	curr_f.resize(n_files_opened);

	MLOG("MedConvert: create_indexes: n_open_in_files %d\n", n_open_in_files);

	//check not dry run
	if (test_run_max_pids == 0) {
		if (open_indexes() < 0) {
			MERR("MedConvert: create_indexes: couldn't open index and data files\n");
			return -1;
		}
	}
	else
		MLOG("###!!! DRY RUN !!!###\n");

	int n_pids_extracted = 0;
	map<pair<string, string>, int> missing_dict_vals;
	vector<int> all_pids;  // a list of all pids in the repository to be written to file.
	all_pids.push_back(0); // reserved place for later placing of total number of pids
	MedTimer timer_action, inside_timer;
	vector<double> tot_time(5, 0);
	int curr_errors = 0;
	map<string, int> prev_forced_errs = missing_forced_signals;
	//stores in memory next lines. first index is file id. second is line
	vector<vector<string>> file_to_lines(n_files_opened);
	vector<int> file_buffer_pos(n_files_opened);
	if (read_lines_buffer > 0)
		for (size_t i = 0; i < n_files_opened; ++i)
			file_to_lines[i].reserve(read_lines_buffer);
	//read first buffers - without parallel:
	MLOG("Reading first line buffer for all files\n");

	MedProgress prog_read_buf("Buffering files", n_files_opened, 15, 1);
	prog_read_buf.alway_print_total = true;
	for (size_t i = 0; i < n_files_opened; ++i) {
		read_file_to_buffer(infs[i], file_to_lines[i], read_lines_buffer);
		prog_read_buf.update();
	}
	if (!full_error_file.empty()) {
		err_log_file.open(full_error_file);
		if (!err_log_file.good())
			MTHROW_AND_ERR("Error - unable to open error log file %s in write mode\n", full_error_file.c_str());
	}

	MedProgress load_progress("MedConvert::create_indexes", 0, 30);


	while (n_open_in_files > 0) {

		// find current pid to extract
		c_pid = -1;
		for (int i = 0; i < n_files_opened; i++) {
			if (pid_in_file[i] > 0) {
				if (c_pid < 0)
					c_pid = pid_in_file[i];
				else if (pid_in_file[i] < c_pid)
					c_pid = pid_in_file[i];
			}
		}

		if (c_pid % 10000 == 0) {
			MLOG("Current pid to extract is %d <<<<< >>>>> n_extracted %d n_open_in_files %d. Times [%2.1f (%2.1f, %2.1f), %2.1f, %2.1f]\n",
				c_pid, n_pids_extracted, n_open_in_files, tot_time[0], tot_time[3], tot_time[4], tot_time[1], tot_time[2]);
			if (err_log_file.is_open()) err_log_file.flush();
		}

		// read data from files
		curr.raw_data.clear();
		curr.raw_data.resize(serial2sid.size());
		curr.pid = c_pid;

		timer_action.start();

		inside_timer.start();
		vector<string> lines;
		vector<int> f_i;

#pragma omp parallel for schedule(dynamic) if (run_parallel_files)
		for (int i = 0; i < n_files_opened; i++) {
			int fpid = c_pid;
			if (pid_in_file[i] >= -1 && pid_in_file[i] <= c_pid) {
				collect_lines(lines, f_i, i, file_to_lines[i], file_buffer_pos[i], infs[i], file_type[i], curr, fpid, fstats[i], missing_dict_vals);
				pid_in_file[i] = fpid; // current pid after the one we wanted
			}
		}

		inside_timer.take_curr_time();
		tot_time[3] += inside_timer.diff_sec();

		inside_timer.start();
		get_next_signal_all_lines(lines, f_i, curr, fstats, missing_dict_vals);
		inside_timer.take_curr_time();
		tot_time[4] += inside_timer.diff_sec();


		timer_action.take_curr_time();
		tot_time[0] += timer_action.diff_sec();

		// write data to output files
		timer_action.start();
		if (curr.pid >= 0) {
			if (write_indexes_new_modes(curr) < 0) {
				//MERR("MedConvert: create_indexes: curr packet for pid %d was not written...\n", curr.pid);

			}
			else
				all_pids.push_back(curr.pid);
		}
		timer_action.take_curr_time();
		tot_time[1] += timer_action.diff_sec();

		++n_pids_extracted;
		load_progress.update();
		if (check_for_error_pid_cnt > 0 && n_pids_extracted % check_for_error_pid_cnt == 0) {
			timer_action.start();
			test_for_load_error(missing_dict_vals, n_pids_extracted, false, curr_errors, curr_errors,
				prev_forced_errs);
			prev_forced_errs = missing_forced_signals;
			timer_action.take_curr_time();
			tot_time[2] += timer_action.diff_sec();
		}

		if (test_run_max_pids > 0 && n_pids_extracted >= test_run_max_pids)
			break;
	}

	MLOG("Current pid to extract is %d <<<<< >>>>> n_extracted %d n_open_in_files %d. Times [%2.1f (%2.1f, %2.1f), %2.1f, %2.1f]\n",
		c_pid, n_pids_extracted, n_open_in_files, tot_time[0], tot_time[3], tot_time[4], tot_time[1], tot_time[2]);
	if (err_log_file.is_open()) err_log_file.flush();

	if (test_run_max_pids > 0) {
		//close all input files:
		for (int i = 0; i < n_files_opened; i++)
			if (infs[i].is_open())
				infs[i].close();
	}
	if (!full_error_file.empty())
		err_log_file.close();
	map<string, int> empty_cnts;
	test_for_load_error(missing_dict_vals, n_pids_extracted, true, 0, curr_errors, empty_cnts);

	MLOG("Finished reading all pids (%d pids extracted) - closing index and data files\n", n_pids_extracted);
	if (mode < 3)
		close_indexes();
	else {
		if (test_run_max_pids == 0)
			write_all_indexes(all_pids);
	}

	return 0;
}



//------------------------------------------------
void MedConvert::test_for_load_error(const map<pair<string, string>, int> &missing_dict_vals,
	int n_pids_extracted, bool final_test,
	int prev_total_missings, int &total_missing, const map<string, int> &prev_missing_forced_signals) const {
	total_missing = 0;
	for (auto& entry : missing_dict_vals) {
		total_missing += entry.second;
		if (prev_total_missings < total_missing)
			MWARN("MedConvert: saw missing entry [%s]:[%s] %d times, total %d missing\n", entry.first.first.c_str(),
				entry.first.second.c_str(), entry.second, total_missing);
		if (safe_mode && total_missing > allowed_unknown_catgory_cnt) {
			MTHROW_AND_ERR("%d > %d missing entries is too much... refusing to create repo!\n", 
				total_missing, allowed_unknown_catgory_cnt);
		}
	}
	for (auto& entry : missing_forced_signals) {
		if (prev_missing_forced_signals.find(entry.first) == prev_missing_forced_signals.end() ||
			prev_missing_forced_signals.at(entry.first) < entry.second)
			MWARN("MedConvert: saw missing_forced_signal [%s] %d times\n", entry.first.c_str(), entry.second);
		if (n_pids_extracted > 0 && safe_mode &&
			(double(entry.second) / n_pids_extracted > allowed_missing_pids_from_forced_ratio ||
			(allowed_missing_pids_from_forced_cnt > 0 && entry.second > allowed_missing_pids_from_forced_cnt)))
			MTHROW_AND_ERR("%d / %d missing_forced_signal is too much... refusing to create repo!\n", entry.second, n_pids_extracted);
	}
	// all files are closed, all are written correctly

	// print statistics for data files
	if (final_test)
		MLOG("Statistics for %d data files\n", fstats.size());
	for (auto& stat : fstats) {
		float ratio = (float)(stat.n_parsed_lines + 1) / (float)(stat.n_relevant_lines + 1);
		float bad_ratio = (float)(stat.n_bad_format_lines + 1) / (float)(stat.n_relevant_lines + 1);
		if (final_test)
			MLOG("file [%d] : %s : n_lines %d , n_relevant_lines %d , n_bad_format_lines %d n_parsed_lines %d : parsed %g\n",
				stat.id, stat.fname.c_str(), stat.n_lines, stat.n_relevant_lines, stat.n_bad_format_lines, stat.n_parsed_lines,
				ratio);
		if (ratio < min_parsed_line_ratio || bad_ratio > max_bad_line_ratio) {
			if (stat.n_relevant_lines > 1000) {
				MTHROW_AND_ERR("%d/%d lines loaded for file [%s]\n", stat.n_parsed_lines, stat.n_relevant_lines, stat.fname.c_str());
			}
			else {
				if (final_test)
					MWARN("%d/%d lines loaded for file [%s]\n", stat.n_parsed_lines, stat.n_relevant_lines, stat.fname.c_str());
			}
		}

	}
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
			if (verbose_open_files)
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
/*
auto cd_to_tuple(const collected_data & v1) -> decltype(std::tie(v1.date, v1.date2, v1.time, v1.time, v1.time2, v1.val, v1.longVal, v1.val1, v1.val2, v1.val3, v1.val4, v1.f_val2)) {
	return std::tie(v1.date, v1.date2, v1.time, v1.time, v1.time2, v1.val, v1.longVal, v1.val1, v1.val2, v1.val3, v1.val4, v1.f_val2);
}
*/

//------------------------------------------------
int MedConvert::write_indexes_new_modes(pid_data &curr)
{
	if (curr.pid < 0)
		MTHROW_AND_ERR("MedConvert::write_indexes negative pid %d", curr.pid);
	// sort unique per signal
	int i;
	int n_sids = (int)curr.raw_data.size();
#pragma omp parallel for if (n_sids > 4) schedule(dynamic)
	for (i = 0; i < curr.raw_data.size(); i++) {
		GenericSigVec gsv1;
		auto& info = sigs.Sid2Info[serial2siginfo[i].sid];
		gsv1.init(info);
		// sort
		sort(curr.raw_data[i].begin(), curr.raw_data[i].end(),
			[&](const collected_data &v1, const collected_data &v2)
		{
			return gsv1.compareTimeLt(&v1.buf[0], 0, &v2.buf[0], 0);
		});

		// gettomg rid of duplicates
		vector<collected_data>::iterator it;
		int struct_size = sigs.Sid2Info[serial2siginfo[i].sid].bytes_len;
		it = unique(curr.raw_data[i].begin(), curr.raw_data[i].end(), [=](const collected_data &v1, const collected_data &v2) {
			return memcmp(v1.buf, v2.buf, struct_size) == 0;
			//return cd_to_tuple(v1) == cd_to_tuple(v2); 
		});
		curr.raw_data[i].resize(distance(curr.raw_data[i].begin(), it));

	}

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
#pragma omp parallel for schedule(dynamic)
	for (i = 0; i < curr.raw_data.size(); i++) {

		int ilen = (int)curr.raw_data[i].size();
		if (ilen > 0) {
			int sid = serial2sid[i];
			int sid_type = serial2siginfo[i].type; //sigs.type(sid);
			int fno = serial2siginfo[i].fno; //sid2fno[sid];

			if (sid_type >= 0 && sid_type < T_Last) {

				int struct_len = (int)sigs.Sid2Info[sid].bytes_len;
				for (int j = 0; j < ilen; j++) {
					data_f[fno]->write((char *)&(curr.raw_data[i][j].buf[0]), struct_len);
				}

				indexes[fno].insert(curr.pid, ilen);
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

void MedConvert::init_load_params(const string &init_str) {
	map<string, string> mapper;
	if (MedSerialize::init_map_from_string(init_str, mapper) < 0)
		MTHROW_AND_ERR("Error Init from String %s\n", init_str.c_str());

	for (const auto &it : mapper)
	{
		if (it.first == "check_for_error_pid_cnt")
			check_for_error_pid_cnt = med_stoi(it.second);
		else if (it.first == "test_run_max_pids")
			test_run_max_pids = med_stoi(it.second);
		else if (it.first == "allowed_missing_pids_from_forced_ratio")
			allowed_missing_pids_from_forced_ratio = med_stof(it.second);
		else if (it.first == "max_bad_line_ratio")
			max_bad_line_ratio = med_stof(it.second);
		else if (it.first == "min_parsed_line_ratio")
			min_parsed_line_ratio = med_stof(it.second);
		else if (it.first == "allowed_unknown_catgory_cnt")
			allowed_unknown_catgory_cnt = med_stoi(it.second);
		else if (it.first == "allowed_missing_pids_from_forced_cnt")
			allowed_missing_pids_from_forced_cnt = med_stoi(it.second);
		else if (it.first == "read_lines_buffer")
			read_lines_buffer = med_stoi(it.second);
		else if (it.first == "verbose_open_files")
			verbose_open_files = med_stoi(it.second) > 0;
		else if (it.first == "run_parallel")
			run_parallel = (med_stoi(it.second) > 0);
		else if (it.first == "run_parallel_files")
			run_parallel_files = (med_stoi(it.second) > 0);
		else if (it.first == "full_error_file")
			full_error_file = it.second;
		else
			MTHROW_AND_ERR("Error in MedConvert::init_load_params - unknown parameter %s\n",
				it.first.c_str());
	}

}