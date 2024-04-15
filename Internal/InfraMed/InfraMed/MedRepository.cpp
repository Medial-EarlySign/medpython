//
// MedRepository.c
//

#define __INFRAMED_DLL

#ifndef _SCL_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS
#endif

#include "InfraMed.h"
#include "Utils.h"
#include "MedPidRepository.h"
#include <fstream>
#include <Logger/Logger/Logger.h>
#include <MedUtils/MedUtils/MedUtils.h>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <omp.h>

#define LOCAL_SECTION LOG_REP
#define LOCAL_LEVEL	LOG_DEF_LEVEL
extern MedLogger global_logger;


mutex index_table_locks[MAX_SID_NUMBER];
mutex index_table_read_locks[MAX_SID_NUMBER];

//===========================================================
// MedRepository
//===========================================================

//-----------------------------------------------------------
void MedRepository::clear()
{
	desc = "";
	config_fname = "";
	dictionary_fnames.clear();
	signal_fnames.clear();
	data_fnames.clear();
	index_fnames.clear();
	format = 0;
	time_unit = MedTime::Date;
	index.clear();
	if (dict.read_state != 1) dict.clear();
	free_all_sigs();
	sigs.clear();
	pids.clear();
	if (work_area != NULL)
		delete[] work_area;
	work_area = NULL;
}
//-----------------------------------------------------------
int MedRepository::read_config(const string &fname)
{
	string fixed_name;
	//if (fix_path(fname,fixed_name) != 0) {
	//	fprintf(stderr,"MedRepository : read_config: Can't fix path name %s\n",fname.c_str()) ;
	//	fixed_name = fname;
	//	//return -1 ;
	//}

	fixed_name = fname;
	ifstream inf(fixed_name);

	if (!inf) {
		MERR("MedRepository: read_config: Can't open file [%s]\n", fixed_name.c_str());
		return -1;
	}

	clear();

	string curr_line;
	while (getline(inf, curr_line)) {
		if ((curr_line.size() > 1) && (curr_line[0] != '#')) {

			if (curr_line[curr_line.size() - 1] == '\r')
				curr_line.erase(curr_line.size() - 1);

			vector<string> fields;
			split(fields, curr_line, boost::is_any_of("\t"));

			if (fields[0].compare("DESCRIPTION") == 0) {
				desc = fields[1];
			}
			else if (fields[0].compare("DIR") == 0) {
				path = fields[1];
				if (path.compare(".") == 0) {
					// in this case we fix our path to be from where we were called
					size_t found = fixed_name.find_last_of("/\\");
					path = fixed_name.substr(0, found);
				}
			}
			else if (fields[0].compare("METADATA_DIR") == 0) {
				metadata_path = fields[1];
				if (metadata_path.compare(".") == 0) {
					// in this case we fix our path to be from where we were called
					size_t found = fixed_name.find_last_of("/\\");
					metadata_path = fixed_name.substr(0, found);
				}
			}
			else if (fields[0].compare("DICTIONARY") == 0) {
				dictionary_fnames.push_back(fields[1]);
			}
			else if (fields[0].compare("SIGNAL") == 0) {
				signal_fnames.push_back(fields[1]);
				dictionary_fnames.push_back(fields[1]);
			}
			else if (fields[0].compare("SFILES") == 0) {
				fsignals_to_files = fields[1];
			}
			else if (fields[0].compare("MODE") == 0) {
				rep_mode = stoi(fields[1]);
				index.rep_mode = rep_mode;
			}
			else if (fields[0].compare("PREFIX") == 0) {
				rep_files_prefix = fields[1];
			}
			else if (fields[0].compare("TIME_UNIT") == 0 || fields[0].compare("TIMEUNIT") == 0) {
				time_unit = med_stoi(fields[1]);
			}
			else if (fields[0].compare("DATA") == 0) {
				if (fields.size() < 3) {
					MERR("NedRepository: read_config: DATA line in wrong format: %s \n", curr_line.c_str());
					return -1;
				}
				int fno = stoi(fields[1]);
				if (data_fnames.size() < fno + 1)
					data_fnames.resize(fno + 1);
				data_fnames[fno] = fields[2];
			}
			else if (fields[0].compare("INDEX") == 0) {
				index_fnames.push_back(fields[1]);
			}
			else MLOG("Ignoring line: [%s]\n", curr_line.c_str());

		}
	}
	if (path != "" && metadata_path == "")
		metadata_path = path;
	// adding path to names + fixing paths.
	if (add_path_to_name_IM(metadata_path, dictionary_fnames) == -1 || add_path_to_name_IM(metadata_path, signal_fnames) == -1 || add_path_to_name_IM(path, data_fnames) == -1 || add_path_to_name_IM(path, index_fnames) == -1)
		return -1;

	config_fname = fname;
	inf.close();

	for (const string &dict_path : _addtional_dict_path)
		dictionary_fnames.push_back(dict_path);
	return 0;
}

//-----------------------------------------------------------
int MedRepository::read_dictionary()
{
	if (dictionary_fnames.size() > 0) {
		for (int i = 0; i < dictionary_fnames.size(); i++) {
			if (read_dictionary(dictionary_fnames[i]) < 0)
				return -1;
		}
	}

	return 0;
}

//-----------------------------------------------------------
int MedRepository::read_dictionary(const string &fname)
{
	return(dict.read(fname));
}

void MedRepository::load_additional_dict(const string &dict_path) {
	MLOG("Adding new dictionary in runtime %s\n", dict_path.c_str());
	_addtional_dict_path.push_back(dict_path);
}

void MedRepository::clear_additional_dict() {
	_addtional_dict_path.clear();
}

//-----------------------------------------------------------
int MedRepository::read_signals()
{
	if (signal_fnames.size() > 0) {
		for (int i = 0; i < signal_fnames.size(); i++) {
			if (read_signals(signal_fnames[i]) < 0)
				return -1;
		}
	}
	return 0;
}

//-----------------------------------------------------------
int MedRepository::read_signals(const string &fname)
{
	return(sigs.read(fname));
}

//-----------------------------------------------------------
int MedRepository::read_index()
{
	return(read_index(vector<int>(), vector<int>()));
}

//-----------------------------------------------------------
int MedRepository::read_index(const vector<int> &pids_to_take, const vector<int> &signals_to_take)
{
	if (index_fnames.size() > 0) {
		return (read_index(index_fnames, pids_to_take, signals_to_take));
	}
	return 0;
}

//-----------------------------------------------------------
int MedRepository::read_index(vector<string> &fnames)
{
	return(read_index(fnames, vector<int>(), vector<int>()));
}

//-----------------------------------------------------------
int MedRepository::read_index(vector<string> &fnames, const vector<int> &pids_to_take, const vector<int> &signals_to_take)
{
	return(index.read_sub_index(fnames, pids_to_take, signals_to_take));
}

//-----------------------------------------------------------
int MedRepository::get_data_mode(const string &fname)
{
	ifstream inf;

	inf.open(fname, ios::in | ios::binary);

	if (!inf) {
		MERR("MedRepository: get_data_mode: Can't open file %s\n", fname.c_str());
		return -1;
	}

	int mode;
	inf.read((char *)&mode, sizeof(int));

	inf.close();

	return mode;
}

//-----------------------------------------------------------
void MedRepository::free_data()
{
	if (work_area != NULL) {
		delete[] work_area;
		work_area = NULL;
	}
	if (work_size > 0)
		work_size = 0;
	index.set_mem_ptrs_off();
}

//-----------------------------------------------------------
int MedRepository::read_data()
{
	return(read_data(data_fnames));
}

//-----------------------------------------------------------
int MedRepository::read_data(const string &fname)
{
	vector<string> names;
	names.push_back(fname);
	return(read_data(names));
}

//-----------------------------------------------------------
int MedRepository::read_data(vector<string> &fnames)
{
	free_data();
	return(index.read_all_data(work_area, work_size, fnames));
}

//-----------------------------------------------------------
int MedRepository::read_all(const string &conf_name)
{
	return(read_all(conf_name, vector<int>(), vector<int>(), 2)); // default: read all in the fastest possible way
}

//-----------------------------------------------------------
int MedRepository::read_all(const string &conf_fname, const vector<int> &pids_to_take, const vector<int> &signals_to_take)
{
	return(read_all(conf_fname, pids_to_take, signals_to_take, 1));
}
//-----------------------------------------------------------
int MedRepository::read_all(const string &conf_fname, const vector<int> &pids_to_take, const vector<string> &signals_to_take)
{
	return(read_all(conf_fname, pids_to_take, signals_to_take, 1));
}
//-----------------------------------------------------------
int MedRepository::read_all(const string &conf_fname, const vector<int> &pids_to_take, const vector<string> &signals_to_take, int read_data_flag)
{
	free_data();
	clear();

	// read config
	if (read_config(conf_fname) < 0) {
		MERR("MedRepository: read_all: error: read_config [%s] failed\n", conf_fname.c_str());
		return -1;
	}

	MLOG("MedRepository: read config file %s\n", conf_fname.c_str());

	// read dictionaries
	if (dict.read(dictionary_fnames) < 0) {
		MERR("MedRepository: read_all: error: read dictionary failed\n");
		return -1;
	}
	dict.read_state = 1;
	MLOG_D("MedRepository: read %d dictionary files\n", dictionary_fnames.size());

	MLOG("MedRepository: reading signals: ");
	// Get signals ids
	vector<int> signal_ids_to_take;
	for (unsigned int i = 0; i < signals_to_take.size(); i++) {
		MLOG("%s,", signals_to_take[i].c_str());
		if (dict.id(signals_to_take[i]) < 0) {
			MERR("MedRepository: requested unknown signal %s\n", signals_to_take[i].c_str());
			return -1;
		}
		else {
			signal_ids_to_take.push_back(dict.id(signals_to_take[i]));
		}
	}
	MLOG("\n");

	return read_all(conf_fname, pids_to_take, signal_ids_to_take, read_data_flag);

}
//-----------------------------------------------------------
int MedRepository::read_pid_list()
{
	if (rep_mode < 3)
		return 0;

	string fname_pids = path + "/" + rep_files_prefix + "_all_pids.list";

	all_pids_list.clear();
	unsigned char *data = NULL;
	unsigned long long size = 0;
	if (read_bin_file_IM(fname_pids, data, size) < 0) {
		MERR("ERROR: Failed reading %s ... it is recommended to repeat conversion\n", fname_pids.c_str());
		return 0; // -1
	}

	int *list = (int *)data;
	all_pids_list.resize(list[0]);
	for (int i = 0; i < all_pids_list.size(); i++)
		all_pids_list[i] = list[i + 1];

	if (size > 0) delete[] data;

	return 0;
}

//-----------------------------------------------------------------------------------
int MedRepository::init(const string &conf_fname)
{
	// read config
	if (read_config(conf_fname) < 0) {
		MERR("MedRepository: init: error: read_config %s failed\n", conf_fname.c_str());
		return -1;
	}

	if (rep_mode < 3) {
		MERR("MedRepository: init: error: init() is possible only from mode 3 and up\n");
		return -1;
	}

	MLOG_D("MedRepository: read config file %s\n", conf_fname.c_str());

	// read dictionaries
	if (dict.read_state == 0) {
		if (dict.read(dictionary_fnames) < 0) {
			MERR("MedRepository: init: error: read dictionary failed\n");
			return -1;
		}
	}
	dict.read_state = 2;
	MLOG_D("MedRepository: read %d dictionary files\n", dictionary_fnames.size());

	// read signals
	if (signal_fnames.size() == 0) {
		MERR("MedRepository: read_all: error: no signals def file given, this is mandatory\n");
		return -1;
	}
	if (sigs.read(signal_fnames) < 0) {
		MERR("MedRepository: read_all: error: read signal files failed\n");
		return -1;
	}

	generate_fnames_for_prefix();
	MedTimer t("Rep Read Time");
	t.start();
	//MLOG("Reading Index Tables\n");
	vector<int> pids_sort_uniq = { -1 }, sids;
	if (read_index_tables(pids_sort_uniq, sids) < 0)
		return -1;
	if (!in_mem_mode_active() && read_pid_list() < 0)
		return -1;

	t.take_curr_time(); MLOG("Read data time %f seconds\n", t.diff_sec());
	return 0;
}

//-----------------------------------------------------------
int MedRepository::read_all(const string &conf_fname, const vector<int> &pids_to_take, const vector<int> &signals_to_take, int read_data_flag)
{
	if (dict.read_state != 1) {
		free_data();
		clear();
	}

	// read config
	if (read_config(conf_fname) < 0) {
		MERR("MedRepository: read_all: error: read_config %s failed\n", conf_fname.c_str());
		return -1;
	}

	MLOG_D("MedRepository: read config file %s\n", conf_fname.c_str());

	// read dictionaries
	if (dict.read_state == 0) {
		if (dict.read(dictionary_fnames) < 0) {
			MERR("MedRepository: read_all: error: read dictionary failed\n");
			return -1;
		}
	}
	dict.read_state = 2;
	MLOG_D("MedRepository: read %d dictionary files\n", dictionary_fnames.size());

	// read signals
	if (signal_fnames.size() == 0) {
		MERR("MedRepository: read_all: error: no signals def file given, this is mandatory\n");
		return -1;
	}
	if (sigs.read(signal_fnames) < 0) {
		MERR("MedRepository: read_all: error: read signal files failed\n");
		return -1;
	}

	if (rep_mode < 2) {
		MTHROW_AND_ERR("MedRepository: rep_mode %d is no longer supported\n", rep_mode);
	}
	else {
		// mode 2 TBD !!!
		generate_fnames_for_prefix();
	}

	MedTimer t("Rep Read Time");
	t.start();
	if (rep_mode < 3) {
		// read index
		if (read_index(pids_to_take, signals_to_take) < 0) {
			MERR("MedRepository: read_all: failed reading index files\n");
			return -1;
		}
		MLOG_D("MedRepository: read %d index files\n", index_fnames.size());
		t.take_curr_time(); MLOG("Read Index time %f seconds\n", t.diff_sec()); t.start();

		// read data - this will read the data that matches the sub index read
		if (read_data_flag == 1) {
			if (read_data(data_fnames) < 0) {
				MERR("MedRepository: read_all: failed reading data files\n");
				return -1;
			}
			MLOG_D("MedRepository: read %d data files\n", data_fnames.size());
		}
		if (read_data_flag == 2) {
			free_data();
			if (index.read_full_data(work_area, work_size, data_fnames) < 0) {
				MERR("MedRepository: read_all: failed reading data files in full mode\n");
				return -1;
			}
			MLOG_D("MedRepository: read %d data files in full mode\n", data_fnames.size());
		}
	}
	else {
		//MLOG("Reading Index Tables\n");
		vector<int> pids_sort_uniq = pids_to_take;
		sort(pids_sort_uniq.begin(), pids_sort_uniq.end());
		auto it = unique(pids_sort_uniq.begin(), pids_sort_uniq.end());
		pids_sort_uniq.resize(distance(pids_sort_uniq.begin(), it));

		if (read_index_tables(pids_sort_uniq, signals_to_take) < 0)
			return -1;

		if (read_pid_list() < 0)
			return -1;

	}

	t.take_curr_time(); MLOG("Read data time %f seconds\n", t.diff_sec());
	return 0;
}

//----------------------------------------------------------
SDateVal *MedRepository::get_before_date(int pid, int sid, int date, int &len)
{
	int orig_len;
	SDateVal *sdv = (SDateVal *)get(pid, sid, orig_len);

	if (orig_len == 0)
		sdv = NULL;
	len = 0;
	// ToDo:: Improve to binary serarch
	for (int i = 0; i < orig_len; i++) {
		if (sdv[i].date < date)
			len++;
		else
			break;
	}

	return sdv;
}

//----------------------------------------------------------
SDateVal *MedRepository::get_before_date(int pid, const string &sig_name, int date, int &len)
{
	int sid = sigs.sid(sig_name);
	if (sid < 0) {
		len = 0;
		return NULL;
	}
	return(get_before_date(pid, sid, date, len));
}

//----------------------------------------------------------
SDateVal *MedRepository::get_date(int pid, const string &sig_name, int date, const string &mode)
{
	int sid = sigs.sid(sig_name);
	if (sid < 0)
		return NULL;
	return(get_date(pid, sid, date, mode));
}

//----------------------------------------------------------
SDateVal *MedRepository::get_date(int pid, int sid, int date, const string &mode)
{
	SDateVal *sdv;

	int len;
	sdv = (SDateVal *)get(pid, sid, len);
	if (sdv == NULL)
		return NULL;

	SDateVal *sdv_res = NULL;
	if (mode.compare("==") == 0) {

		// ToDo:: move to binary search
		for (int i = 0; i < len; i++) {
			if (sdv[i].date == date)
				return(&sdv[i]);
			else if (sdv[i].date > date)
				return NULL;
		}

	}

	else if (mode.compare("<=") == 0) {

		// ToDo:: move to binary search
		for (int i = 0; i < len; i++) {
			if (sdv[i].date <= date)
				sdv_res = &sdv[i];
			else if (sdv[i].date > date)
				return sdv_res;
		}
	}

	else if (mode.compare(">=") == 0) {

		// ToDo:: move to binary search
		for (int i = len - 1; i > 0; i++) {
			if (sdv[i].date >= date)
				sdv_res = &sdv[i];
			else if (sdv[i].date < date)
				return sdv_res;
		}
	}

	else if (mode.compare("<") == 0) {

		// ToDo:: move to binary search
		for (int i = 0; i < len; i++) {
			if (sdv[i].date < date)
				sdv_res = &sdv[i];
			else if (sdv[i].date >= date)
				return sdv_res;
		}
	}

	else if (mode.compare(">") == 0) {

		// ToDo:: move to binary search
		for (int i = len - 1; i > 0; i++) {
			if (sdv[i].date > date)
				sdv_res = &sdv[i];
			else if (sdv[i].date <= date)
				return sdv_res;
		}
	}

	return NULL;
}

string MedRepository::convert_date(int d, int sid) {
	return med_time_converter.convert_times_S(
		sigs.Sid2Info[sid].time_unit, MedTime::DateTimeString, d);
}
void MedRepository::print_channel_helper(int sid, int channel, float val) {
	if (sigs.is_categorical_channel(sid, channel)) {
		MOUT(" %d ", (int)val);
		int section_id = dict.section_id(sigs.name(sid));
		int drug_sid = sigs.sid("Drug");
		int drugs_nice_names_section = dict.section_id("Drugs_nice_names");
		if (sid == drug_sid && drugs_nice_names_section != 0) {
			// ugly hack: Drug names do not specify ATC category, so we use Drugs_nice_names instead
			MLOG_D("replacing Drug=[%d] with Drugs_nice_names=[%d]\n", section_id, drugs_nice_names_section);
			section_id = drugs_nice_names_section;
		}
		int names_printed = 0;
		if (dict.dict(section_id)->Id2Names.find(val) != dict.dict(section_id)->Id2Names.end())
			for (int j = 0; j < dict.dict(section_id)->Id2Names[val].size(); j++) {
				string st = dict.dict(section_id)->Id2Names[val][j];
				MOUT("|%s", st.c_str());
				if (++names_printed == 3)
					break;
			}
	}
	else
		MOUT(" %f ", val);
}
string MedRepository::get_channel_info(int sid, int channel, float val) {
	string info;
	if (sigs.is_categorical_channel(sid, channel)) {
		info = to_string((int)val);
		int section_id = dict.section_id(sigs.name(sid));
		int drug_sid = sigs.sid("Drug");
		int drugs_nice_names_section = dict.section_id("Drugs_nice_names");
		if (sid == drug_sid && drugs_nice_names_section != 0) {
			// ugly hack: Drug names do not specify ATC category, so we use Drugs_nice_names instead
			MLOG_D("replacing Drug=[%d] with Drugs_nice_names=[%d]\n", section_id, drugs_nice_names_section);
			section_id = drugs_nice_names_section;
		}
		int names_printed = 0;
		if (dict.dict(section_id)->Id2Names.find(val) != dict.dict(section_id)->Id2Names.end())
			for (int j = 0; j < dict.dict(section_id)->Id2Names[val].size(); j++) {
				string st = dict.dict(section_id)->Id2Names[val][j];
				info += "|" + st;
				if (++names_printed == 3)
					break;
			}
	}
	else
		info = to_string(val);

	return info;
}
//-----------------------------------------------------------
void MedRepository::print_vec_dict(void *data, int len, int pid, int sid)
{
	MOUT("pid %d sid %d %s", pid, sid, sigs.name(sid).c_str());
	for (int channel = 0; channel < sigs.Sid2Info[sid].n_val_channels; channel++) {
		string u = sigs.unit_of_measurement(sid, 0);
		if (u != "")
			MOUT(" %s ", u.c_str());
	}
	MOUT(" ::");

	MOUT(" %d :: ", len);
	if (sigs.has_any_categorical_channel(sid))
		MOUT("\n");

	for (int i = 0; i < len; i++) {
		if (sigs.type(sid) == T_Value) {
			SVal *v = (SVal *)data;
			print_channel_helper(sid, 0, v[i].val);
		}
		else if (sigs.type(sid) == T_DateVal) {
			SDateVal *v = (SDateVal *)data;
			MOUT(" %s ", convert_date(v[i].date, sid).c_str());
			print_channel_helper(sid, 0, v[i].val);
		}
		else if (sigs.type(sid) == T_TimeVal) {
			STimeVal *v = (STimeVal *)data;
			MOUT(" %lld ", v[i].time);
			print_channel_helper(sid, 0, v[i].val);
		}
		else if (sigs.type(sid) == T_DateRangeVal) {
			SDateRangeVal *v = (SDateRangeVal *)data;
			MOUT(" %s %s ", convert_date(v[i].date_start, sid).c_str(), convert_date(v[i].date_end, sid).c_str());
			print_channel_helper(sid, 0, v[i].val);
		}
		else if (sigs.type(sid) == T_TimeStamp) {
			STimeStamp *v = (STimeStamp*)data;
			MOUT(" %s ", convert_date(v[i].time, sid).c_str());
		}
		else if (sigs.type(sid) == T_TimeRangeVal) {
			STimeRangeVal *v = (STimeRangeVal *)data;
			MOUT(" %lld - %lld  ", v[i].time_start, v[i].time_end);
			print_channel_helper(sid, 0, v[i].val);
		}
		else if (sigs.type(sid) == T_DateVal2) {
			SDateVal2 *v = (SDateVal2 *)data;
			MOUT(" %s ", convert_date(v[i].date, sid).c_str());
			print_channel_helper(sid, 0, v[i].val);
			print_channel_helper(sid, 1, v[i].val2);
		}

		else if (sigs.type(sid) == T_DateShort2) {
			SDateShort2 *v = (SDateShort2 *)data;
			MOUT(" %s ", convert_date(v[i].date, sid).c_str());
			print_channel_helper(sid, 0, v[i].val1);
			print_channel_helper(sid, 1, v[i].val2);
		}
		else if (sigs.type(sid) == T_ValShort2) {
			SValShort2 *v = (SValShort2 *)data;
			print_channel_helper(sid, 0, v[i].val1);
			print_channel_helper(sid, 1, v[i].val2);
		}

		else if (sigs.type(sid) == T_CompactDateVal) {
			SCompactDateVal *v = (SCompactDateVal *)data;
			MOUT(" %d ", compact_date_to_date(v[i].compact_date));
			print_channel_helper(sid, 0, v[i].val);
		}
		else if (sigs.type(sid) == T_DateRangeVal2) {
			SDateRangeVal2 *v = (SDateRangeVal2 *)data;
			MOUT(" %s %s ", convert_date(v[i].date_start, sid).c_str(), convert_date(v[i].date_end, sid).c_str());
			print_channel_helper(sid, 0, v[i].val);
			print_channel_helper(sid, 1, v[i].val2);
		}
		else if (sigs.type(sid) == T_DateFloat2) {
			SDateFloat2 *v = (SDateFloat2 *)data;
			MOUT(" %s ", convert_date(v[i].date, sid).c_str());
			print_channel_helper(sid, 0, v[i].val);
			print_channel_helper(sid, 1, v[i].val2);
		}
		else if (sigs.type(sid) == T_TimeShort4) {
			STimeShort4 *v = (STimeShort4 *)data;
			MOUT(" %lld ", v[i].time);
			print_channel_helper(sid, 0, v[i].val1);
			print_channel_helper(sid, 1, v[i].val2);
			print_channel_helper(sid, 2, v[i].val3);
			print_channel_helper(sid, 3, v[i].val4);
		}
		else if (sigs.type(sid) == T_Generic) {
			GenericSigVec v;
			v.init(sigs.Sid2Info[sid]);
			v.set_data(data, 1);
			for (int tchan = 0; tchan < v.n_time_channels(); ++tchan) {
				MOUT(" %s ", convert_date(v.Time(i, tchan), sid).c_str());
			}
			for (int vchan = 0; vchan < v.n_val_channels(); ++vchan) {
				print_channel_helper(sid, vchan, v.Val(i, vchan));
			}
		}

		if (sigs.has_any_categorical_channel(sid))
			MOUT(" :\n");
		else
			MOUT(" : ");
	}
	MOUT("\n");
}

//-----------------------------------------------------------
void MedRepository::long_print_vec_dict(void *data, int len, int pid, int sid)
{

	string signal_name = sigs.name(sid);

	for (int i = 0; i < len; i++) {
		string out = to_string(pid) + "\t" + to_string(sid) + "\t" + signal_name + "\t";
		if (sigs.type(sid) == T_Value) {
			SVal *v = (SVal *)data;
			out += "0 " + get_channel_info(sid, 0, v[i].val);
		}
		else if (sigs.type(sid) == T_DateVal) {
			SDateVal *v = (SDateVal *)data;
			out += convert_date(v[i].date, sid) + " " + get_channel_info(sid, 0, v[i].val);
		}
		else if (sigs.type(sid) == T_TimeVal) {
			STimeVal *v = (STimeVal *)data;
			out += to_string(v[i].time) + " " + get_channel_info(sid, 0, v[i].val);
		}
		else if (sigs.type(sid) == T_DateRangeVal) {
			SDateRangeVal *v = (SDateRangeVal *)data;
			out += convert_date(v[i].date_start, sid) + " " + convert_date(v[i].date_end, sid) + " " + get_channel_info(sid, 0, v[i].val);
		}
		else if (sigs.type(sid) == T_TimeStamp) {
			STimeStamp *v = (STimeStamp*)data;
			out += convert_date(v[i].time, sid);
		}
		else if (sigs.type(sid) == T_TimeRangeVal) {
			STimeRangeVal *v = (STimeRangeVal *)data;
			out += to_string(v[i].time_start) + " " + to_string(v[i].time_end) + " " + get_channel_info(sid, 0, v[i].val);
		}
		else if (sigs.type(sid) == T_DateVal2) {
			SDateVal2 *v = (SDateVal2 *)data;
			out += convert_date(v[i].date, sid) + " " + get_channel_info(sid, 0, v[i].val) + " " + get_channel_info(sid, 1, v[i].val2);
		}
		else if (sigs.type(sid) == T_DateShort2) {
			SDateShort2 *v = (SDateShort2 *)data;
			out += convert_date(v[i].date, sid) + " " + get_channel_info(sid, 0, v[i].val1) + " " + get_channel_info(sid, 1, v[i].val2);
		}
		else if (sigs.type(sid) == T_ValShort2) {
			SValShort2 *v = (SValShort2 *)data;
			out += "0 " + get_channel_info(sid, 0, v[i].val1) + " " + get_channel_info(sid, 1, v[i].val2);
		}

		else if (sigs.type(sid) == T_CompactDateVal) {
			SCompactDateVal *v = (SCompactDateVal *)data;
			out += to_string(compact_date_to_date(v[i].compact_date)) + " " + get_channel_info(sid, 1, v[i].val);
		}
		else if (sigs.type(sid) == T_DateRangeVal2) {
			SDateRangeVal2 *v = (SDateRangeVal2 *)data;
			out += convert_date(v[i].date_start, sid) + " " + convert_date(v[i].date_end, sid) + " " + get_channel_info(sid, 0, v[i].val) + " " + get_channel_info(sid, 1, v[i].val2);
		}
		else if (sigs.type(sid) == T_DateFloat2) {
			SDateFloat2 *v = (SDateFloat2 *)data;
			out += convert_date(v[i].date, sid) + " " + get_channel_info(sid, 0, v[i].val) + " " + get_channel_info(sid, 1, v[i].val2);
		}
		else if (sigs.type(sid) == T_TimeShort4) {
			STimeShort4 *v = (STimeShort4 *)data;
			out += to_string(v[i].time) + " "
				+ get_channel_info(sid, 0, v[i].val1) + " "
				+ get_channel_info(sid, 1, v[i].val2) + " "
				+ get_channel_info(sid, 2, v[i].val3) + " "
				+ get_channel_info(sid, 3, v[i].val4);
		}
		else if (sigs.type(sid) == T_Generic) {
			GenericSigVec v;
			v.init(sigs.Sid2Info[sid]);
			v.set_data(data, 1);
			bool is_first = true;
			for (int tchan = 0; tchan < v.n_time_channels(); ++tchan) {
				if (is_first) is_first = false;
				else out += " ";
				out += to_string(v.Time(i, tchan));
			}
			for (int vchan = 0; vchan < v.n_val_channels(); ++vchan) {
				if (is_first) is_first = false;
				else out += " ";
				out += get_channel_info(sid, vchan, v.Val(i, vchan));
			}
		}


		MOUT("%s\n", out.c_str());
	}

}

//-----------------------------------------------------------
void MedRepository::long_print_vec_dict(void *data, int len, int pid, int sid, int from, int to)
{

	string signal_name = sigs.name(sid);

	for (int i = 0; i < len; i++) {
		string out = to_string(pid) + "\t" + to_string(sid) + "\t" + signal_name + "\t";
		if (sigs.type(sid) == T_Value) {
			SVal *v = (SVal *)data;
			out += "0 " + get_channel_info(sid, 0, v[i].val);
		}
		else if (sigs.type(sid) == T_DateVal) {
			SDateVal *v = (SDateVal *)data;
			if (v[i].date < from || v[i].date > to) continue;
			out += convert_date(v[i].date, sid) + " " + get_channel_info(sid, 0, v[i].val);
		}
		else if (sigs.type(sid) == T_TimeVal) {
			STimeVal *v = (STimeVal *)data;
			out += to_string(v[i].time) + " " + get_channel_info(sid, 0, v[i].val);
		}
		else if (sigs.type(sid) == T_DateRangeVal) {
			SDateRangeVal *v = (SDateRangeVal *)data;
			if (v[i].date_start < from || v[i].date_start > to) continue;
			out += convert_date(v[i].date_start, sid) + " " + convert_date(v[i].date_end, sid) + " " + get_channel_info(sid, 0, v[i].val);
		}
		else if (sigs.type(sid) == T_TimeStamp) {
			STimeStamp *v = (STimeStamp*)data;
			out += convert_date(v[i].time, sid);
		}
		else if (sigs.type(sid) == T_TimeRangeVal) {
			STimeRangeVal *v = (STimeRangeVal *)data;
			if (v[i].time_start < from || v[i].time_start > to) continue;
			out += to_string(v[i].time_start) + " " + to_string(v[i].time_end) + " " + get_channel_info(sid, 0, v[i].val);
		}
		else if (sigs.type(sid) == T_DateVal2) {
			SDateVal2 *v = (SDateVal2 *)data;
			if (v[i].date < from || v[i].date > to) continue;
			out += convert_date(v[i].date, sid) + " " + get_channel_info(sid, 0, v[i].val) + " " + get_channel_info(sid, 1, v[i].val2);
		}
		else if (sigs.type(sid) == T_DateShort2) {
			SDateShort2 *v = (SDateShort2 *)data;
			if (v[i].date < from || v[i].date > to) continue;
			out += convert_date(v[i].date, sid) + " " + get_channel_info(sid, 0, v[i].val1) + " " + get_channel_info(sid, 1, v[i].val2);
		}
		else if (sigs.type(sid) == T_ValShort2) {
			SValShort2 *v = (SValShort2 *)data;
			out += "0 " + get_channel_info(sid, 0, v[i].val1) + " " + get_channel_info(sid, 1, v[i].val2);
		}

		else if (sigs.type(sid) == T_CompactDateVal) {
			SCompactDateVal *v = (SCompactDateVal *)data;
			if (v[i].compact_date < from || v[i].compact_date > to) continue;
			out += to_string(compact_date_to_date(v[i].compact_date)) + " " + get_channel_info(sid, 1, v[i].val);
		}
		else if (sigs.type(sid) == T_DateRangeVal2) {
			SDateRangeVal2 *v = (SDateRangeVal2 *)data;
			if (v[i].date_start < from || v[i].date_start > to) continue;
			out += convert_date(v[i].date_start, sid) + " " + convert_date(v[i].date_end, sid) + get_channel_info(sid, 0, v[i].val) + " " + get_channel_info(sid, 1, v[i].val2);
		}
		else if (sigs.type(sid) == T_DateFloat2) {
			SDateFloat2 *v = (SDateFloat2 *)data;
			if (v[i].date < from || v[i].date > to) continue;
			out += convert_date(v[i].date, sid).c_str() + get_channel_info(sid, 0, v[i].val) + " " + get_channel_info(sid, 1, v[i].val2);
		}
		else if (sigs.type(sid) == T_TimeShort4) {
			STimeShort4 *v = (STimeShort4 *)data;
			if (v[i].time < from || v[i].time > to) continue;
			out += to_string(v[i].time) + " "
				+ get_channel_info(sid, 0, v[i].val1) + " "
				+ get_channel_info(sid, 1, v[i].val2) + " "
				+ get_channel_info(sid, 2, v[i].val3) + " "
				+ get_channel_info(sid, 3, v[i].val4);
		}
		else if (sigs.type(sid) == T_Generic) {
			GenericSigVec v;
			v.init(sigs.Sid2Info[sid]);
			v.set_data(data, 1);
			if (v.Time(i) < from || v.Time(i) > to) continue;
			bool is_first = true;
			for (int tchan = 0; tchan < v.n_time_channels(); ++tchan) {
				if (is_first) is_first = false;
				else out += " ";
				out += to_string(v.Time(i, tchan));
			}
			for (int vchan = 0; vchan < v.n_val_channels(); ++vchan) {
				if (is_first) is_first = false;
				else out += " ";
				out += get_channel_info(sid, vchan, v.Val(i, vchan));
			}
		}
		MOUT("%s\n", out.c_str());
	}

}

//-----------------------------------------------------------
void MedRepository::long_print_vec_dict(void *data, int len, int pid, int sid, int i)
{

	string signal_name = sigs.name(sid);

	string out = to_string(pid) + "\t" + to_string(sid) + "\t" + signal_name + "\t";
	if (sigs.type(sid) == T_Value) {
		SVal *v = (SVal *)data;
		out += "0 " + get_channel_info(sid, 0, v[i].val);
	}
	else if (sigs.type(sid) == T_DateVal) {
		SDateVal *v = (SDateVal *)data;
		out += convert_date(v[i].date, sid) + " " + get_channel_info(sid, 0, v[i].val);
	}
	else if (sigs.type(sid) == T_TimeVal) {
		STimeVal *v = (STimeVal *)data;
		out += to_string(v[i].time) + " " + get_channel_info(sid, 0, v[i].val);
	}
	else if (sigs.type(sid) == T_DateRangeVal) {
		SDateRangeVal *v = (SDateRangeVal *)data;
		out += convert_date(v[i].date_start, sid) + " " + convert_date(v[i].date_end, sid) + " " + get_channel_info(sid, 0, v[i].val);
	}
	else if (sigs.type(sid) == T_TimeStamp) {
		STimeStamp *v = (STimeStamp*)data;
		out += convert_date(v[i].time, sid);
	}
	else if (sigs.type(sid) == T_TimeRangeVal) {
		STimeRangeVal *v = (STimeRangeVal *)data;
		out += to_string(v[i].time_start) + " " + to_string(v[i].time_end) + " " + get_channel_info(sid, 0, v[i].val);
	}
	else if (sigs.type(sid) == T_DateVal2) {
		SDateVal2 *v = (SDateVal2 *)data;
		out += convert_date(v[i].date, sid) + " " + get_channel_info(sid, 0, v[i].val) + " " + get_channel_info(sid, 1, v[i].val2);
	}
	else if (sigs.type(sid) == T_DateShort2) {
		SDateShort2 *v = (SDateShort2 *)data;
		out += convert_date(v[i].date, sid) + " " + get_channel_info(sid, 0, v[i].val1) + " " + get_channel_info(sid, 1, v[i].val2);
	}
	else if (sigs.type(sid) == T_ValShort2) {
		SValShort2 *v = (SValShort2 *)data;
		out += "0 " + get_channel_info(sid, 0, v[i].val1) + " " + get_channel_info(sid, 1, v[i].val2);
	}

	else if (sigs.type(sid) == T_CompactDateVal) {
		SCompactDateVal *v = (SCompactDateVal *)data;
		out += to_string(compact_date_to_date(v[i].compact_date)) + " " + get_channel_info(sid, 1, v[i].val);
	}
	else if (sigs.type(sid) == T_DateRangeVal2) {
		SDateRangeVal2 *v = (SDateRangeVal2 *)data;
		out += convert_date(v[i].date_start, sid) + " " + convert_date(v[i].date_end, sid) + get_channel_info(sid, 0, v[i].val) + " " + get_channel_info(sid, 1, v[i].val2);
	}
	else if (sigs.type(sid) == T_DateFloat2) {
		SDateFloat2 *v = (SDateFloat2 *)data;
		out += convert_date(v[i].date, sid).c_str() + get_channel_info(sid, 0, v[i].val) + " " + get_channel_info(sid, 1, v[i].val2);
	}
	else if (sigs.type(sid) == T_TimeShort4) {
		STimeShort4 *v = (STimeShort4 *)data;
		out += to_string(v[i].time) + " "
			+ get_channel_info(sid, 0, v[i].val1) + " "
			+ get_channel_info(sid, 1, v[i].val2) + " "
			+ get_channel_info(sid, 2, v[i].val3) + " "
			+ get_channel_info(sid, 3, v[i].val4);
	}
	else if (sigs.type(sid) == T_Generic) {
		GenericSigVec v;
		v.init(sigs.Sid2Info[sid]);
		v.set_data(data, 1);
		bool is_first = true;
		for (int tchan = 0; tchan < v.n_time_channels(); ++tchan) {
			if (is_first) is_first = false;
			else out += " ";
			out += to_string(v.Time(i, tchan));
		}
		for (int vchan = 0; vchan < v.n_val_channels(); ++vchan) {
			if (is_first) is_first = false;
			else out += " ";
			out += get_channel_info(sid, vchan, v.Val(i, vchan));
		}
	}

	MOUT("%s\n", out.c_str());

}

//-----------------------------------------------------------
void MedRepository::print_data_vec_dict(int pid, int sid)
{


	int len = 0;
	void *data = get(pid, sid, len);

	print_vec_dict(data, len, pid, sid);
}

//-----------------------------------------------------------
void MedRepository::long_print_data_vec_dict(int pid, int sid)
{


	int len = 0;
	void *data = get(pid, sid, len);

	long_print_vec_dict(data, len, pid, sid);
}

void MedRepository::long_print_data_vec_dict(int pid, int sid, int from, int to)
{


	int len = 0;
	void *data = get(pid, sid, len);

	long_print_vec_dict(data, len, pid, sid, from, to);
}

//-----------------------------------------------------------
void MedRepository::print_csv_vec(void *data, int len, int pid, int sid, bool dict_val = false)
{
	int section_id = 0;
	//	int drug_sid;
	if (dict_val) {
		section_id = dict.section_id(sigs.name(sid));
		if (sigs.name(sid) == "Drug") {
			section_id = dict.section_id("Drugs_nice_names");
		}
	}
	for (int i = 0; i < len; i++) {
		int val;
		MOUT("%d,%d,%s,%d,%d,%d,", pid, sid, sigs.name(sid).c_str(), sigs.type(sid), len, i);
		if (sigs.type(sid) == T_Value) {
			SVal *v = (SVal *)data;
			if (dict_val)
				MOUT("%d,", val = (int)v[i].val);
			else MOUT("%f,", v[i].val);
			MOUT("0,0,0,0,0,");
		}
		else if (sigs.type(sid) == T_DateVal) {
			SDateVal *v = (SDateVal *)data;
			if (dict_val)
				MOUT("%d,", val = (int)v[i].val);
			else MOUT("%f,", v[i].val);
			MOUT("0,%d,0,0,0,", v[i].date);
		}
		else if (sigs.type(sid) == T_DateRangeVal) {
			SDateRangeVal *v = (SDateRangeVal *)data;
			if (dict_val)
				MOUT("%d,", val = (int)v[i].val);
			else MOUT("%f,", v[i].val);
			MOUT("0,%d,%d,0,0,", v[i].date_start, v[i].date_end);
		}
		else if (sigs.type(sid) == T_DateRangeVal2) {
			SDateRangeVal2 *v = (SDateRangeVal2 *)data;
			if (dict_val)
				MOUT("%d,", val = (int)v[i].val);
			else MOUT("%f,", v[i].val);
			MOUT("%f,%d,%d,0,0,", v[i].val2, v[i].date_start, v[i].date_end);
		}
		else if (sigs.type(sid) == T_DateFloat2) {
			SDateFloat2 *v = (SDateFloat2 *)data;
			if (dict_val)
				MOUT("%d,", val = (int)v[i].val);
			else MOUT("%f,", v[i].val);
			MOUT("%f,%d,0,0,0,", v[i].val2, v[i].date);
		}
		else if (sigs.type(sid) == T_DateVal2) {
			SDateVal2 *v = (SDateVal2 *)data;
			if (dict_val)
				MOUT("%d,", val = (int)v[i].val);
			else MOUT("%f,", v[i].val);
			MOUT("%d,%d,0,0,0,", v[i].val2, v[i].date);
		}
		else if (sigs.type(sid) == T_DateShort2) {
			SDateShort2 *v = (SDateShort2 *)data;
			val = v[i].val1;
			MOUT("%d,%d,%d,0,0,0,", v[i].val1, v[i].val2, v[i].date);
		}
		else if (sigs.type(sid) == T_ValShort2) {
			SValShort2 *v = (SValShort2 *)data;
			val = v[i].val1;
			MOUT("%d,%d,0,0,0,0,", v[i].val1, v[i].val2);
		}
		else if (sigs.type(sid) == T_ValShort4) {
			SValShort4 *v = (SValShort4 *)data;
			val = v[i].val1;
			MOUT("%d,%d,0,0,%d,%d,", v[i].val1, v[i].val2, v[i].val3, v[i].val4);
		}
		else if (sigs.type(sid) == T_CompactDateVal) {
			SCompactDateVal *v = (SCompactDateVal *)data;
			val = v[i].val;
			MOUT("%d,0,%d,0,0,0,", v[i].val, compact_date_to_date(v[i].compact_date));
		}
		else if (sigs.type(sid) == T_TimeVal) {
			STimeVal *v = (STimeVal *)data;
			if (dict_val)
				MOUT("%d,", val = (int)v[i].val);
			else MOUT("%f,", v[i].val);
			MOUT("0,%d,0,0,0,0,", v[i].time);
		}
		else if (sigs.type(sid) == T_TimeRangeVal) {
			STimeRangeVal *v = (STimeRangeVal *)data;
			if (dict_val)
				MOUT("%d,", val = (int)v[i].val);
			else MOUT("%f,", v[i].val);
			MOUT("0,%lld,%lld,0,0,", v[i].time_start, v[i].time_end, v[i].val);
		}
		else if (sigs.type(sid) == T_TimeShort4) {
			STimeShort4 *v = (STimeShort4 *)data;
			//All values are always int 
			MOUT("%d,%d,%lld,%d,%d,%d,", val = (int)v[i].val1, (int)v[i].val2, v[i].time, 0, (int)v[i].val3, (int)v[i].val4);
		}
		else if (sigs.type(sid) == T_Generic) {
			// The format is:
			// (pid,sid,sname,stype,len,i,) ...  val,val2,date1,date2,val3,val4,desc1,desc2,extra1,extra2
			GenericSigVec v;
			v.init(sigs.Sid2Info[sid]);
			v.set_data(data, 1);
			// val chan #1 (.. or dict val or 0)
			if (v.n_val_channels() > 0) {
				if (dict_val) {
					val = v.Val<int>(i, 0);
					MOUT("%d,", v.Val<int>(i, 0));
				}
				else MOUT("%f,", v.Val(i, 0));
			}
			else MOUT("0,");
			//val chan #2
			if (v.n_val_channels() > 1) {
				MOUT("%f,", v.Val(i, 1));
			}
			else MOUT("0,");
			//time chan #1
			if (v.n_time_channels() > 0) {
				MOUT("%d,", v.Time(i, 0));
			}
			else MOUT("0,");
			//time chan #2
			if (v.n_time_channels() > 1) {
				MOUT("%d,", v.Time(i, 1));
			}
			else MOUT("0,");
			// val chan #3
			if (v.n_val_channels() > 2) {
				MOUT("%f,", v.Val(i, 2));
			}
			else MOUT("0,");
			//val chan #4
			if (v.n_val_channels() > 3) {
				MOUT("%f,", v.Val(i, 3));
			}
			else MOUT("0,");
		}
		if (dict_val && dict.dict(section_id)->Id2Names.find(val) != dict.dict(section_id)->Id2Names.end()) {
			for (size_t j = 0; j < 2; j++)
			{
				if (j < dict.dict(section_id)->Id2Names[val].size()) {
					string st = dict.dict(section_id)->Id2Names[val][j];
					MOUT("\"%s\",", st.c_str());
				}
			}
		}
		MOUT("\n");
	}


}

//-----------------------------------------------------------
void MedRepository::convert_pid_sigs(const UniversalSigVec &usv, vector<pair<vector<string>, vector<string>>> &pid_result, const string &sig_name, int sig_id, int limit_count) {
	pid_result.resize(usv.len);
	char buffer[500];
	int section_id = dict.section_id(sig_name);

	for (int l = 0; l < usv.len; ++l) {
		vector<string> &rec_times = pid_result[l].first;
		vector<string> &rec_vals = pid_result[l].second;
		rec_times.resize(usv.n_time_channels());
		rec_vals.resize(usv.n_val_channels());

		for (int time_ch = 0; time_ch < usv.n_time_channels(); ++time_ch) {
			if (usv.time_channel_types[time_ch] != GenericSigVec::type_enc::INT32) {
				stringstream str_buff;
				if (usv.time_channel_types[time_ch] == GenericSigVec::type_enc::INT64)
					str_buff << usv.Time<long long>(l, time_ch);
				else if (usv.time_channel_types[time_ch] == GenericSigVec::type_enc::UINT64)
					str_buff << usv.Time<unsigned long long>(l, time_ch);
				else if (usv.time_channel_types[time_ch] == GenericSigVec::type_enc::UINT32)
					str_buff << usv.Time<unsigned int>(l, time_ch);
				else if (usv.time_channel_types[time_ch] == GenericSigVec::type_enc::FLOAT32)
					str_buff << usv.Time<float>(l, time_ch);
				else if (usv.time_channel_types[time_ch] == GenericSigVec::type_enc::FLOAT64)
					str_buff << usv.Time<double>(l, time_ch);
				else if (usv.time_channel_types[time_ch] == GenericSigVec::type_enc::FLOAT80)
					str_buff << usv.Time<long double>(l, time_ch);
				else
					str_buff << usv.Time(l, time_ch);
				rec_times[time_ch] = str_buff.str();
			}
			else {
				int raw_time = usv.Time(l, time_ch);
				string format_time = med_time_converter.convert_times_S(
					sigs.Sid2Info[sig_id].time_unit, MedTime::DateTimeString, raw_time);
				if (sigs.Sid2Info[sig_id].time_unit != MedTime::Date)
					snprintf(buffer, sizeof(buffer), "%d|%s", raw_time, format_time.c_str());
				else
					snprintf(buffer, sizeof(buffer), "%d", raw_time);
				format_time = string(buffer);

				rec_times[time_ch] = format_time;
			}
		}

		for (int val_ch = 0; val_ch < usv.n_val_channels(); ++val_ch) {
			bool is_categ = sigs.is_categorical_channel(sig_id, val_ch);
			if (usv.val_channel_types[val_ch] != GenericSigVec::type_enc::INT32
				&& usv.val_channel_types[val_ch] != GenericSigVec::type_enc::FLOAT32
				&& usv.val_channel_types[val_ch] != GenericSigVec::type_enc::INT16
				&& usv.val_channel_types[val_ch] != GenericSigVec::type_enc::INT8
				&& usv.val_channel_types[val_ch] != GenericSigVec::type_enc::UINT16
				&& usv.val_channel_types[val_ch] != GenericSigVec::type_enc::UINT8) {

				stringstream str_buff;
				if (usv.val_channel_types[val_ch] == GenericSigVec::type_enc::FLOAT64)
					str_buff << usv.Val<double>(l, val_ch);
				else if (usv.val_channel_types[val_ch] == GenericSigVec::type_enc::FLOAT80)
					str_buff << usv.Val<long double>(l, val_ch);
				else if (usv.val_channel_types[val_ch] == GenericSigVec::type_enc::INT64)
					str_buff << usv.Val<long long>(l, val_ch);
				else if (usv.val_channel_types[val_ch] == GenericSigVec::type_enc::UINT64)
					str_buff << usv.Val<unsigned long long>(l, val_ch);
				else if (usv.val_channel_types[val_ch] == GenericSigVec::type_enc::UINT32)
					str_buff << usv.Val<unsigned int>(l, val_ch);
				else
					str_buff << usv.Val(l, val_ch);
				rec_vals[val_ch] = str_buff.str();
			}
			else {
				int raw_val = usv.Val<int>(l, val_ch);
				//is_categ |= section_id > 0 && abs(raw_val - usv.Val(l, val_ch)) < 1e-5 && dict.dict(section_id)->Id2Names.find(raw_val) != dict.dict(section_id)->Id2Names.end();
				if (!is_categ) {
					stringstream str_buff;
					if (usv.val_channel_types[val_ch] == GenericSigVec::type_enc::INT32)
						str_buff << usv.Val<int>(l, val_ch);
					else
						str_buff << usv.Val(l, val_ch);
					rec_vals[val_ch] = str_buff.str();
				}
				else {
					stringstream str_buff;
					str_buff << raw_val;
					int printed_cnt = 0;
					if (dict.dict(section_id)->Id2Names.find(raw_val) != dict.dict(section_id)->Id2Names.end())
						for (int j = 0; j < dict.dict(section_id)->Id2Names[raw_val].size()
							&& printed_cnt < limit_count; j++) {
							string st = dict.dict(section_id)->Id2Names[raw_val][j];
							str_buff << "|" << st;
							++printed_cnt;
						}

					rec_vals[val_ch] = str_buff.str();
				}
			}
		}
	}
}

//-----------------------------------------------------------
void MedRepository::print_pid_sig(int pid, const string &sig_name, const vector<pair<vector<string>, vector<string>>> &usv) {
	for (int l = 0; l < usv.size(); ++l) {
		cout << pid << '\t' << sig_name;
		for (int time_ch = 0; time_ch < usv[l].first.size(); ++time_ch)
			cout << '\t' << "Time_ch_" << time_ch << '\t' << usv[l].first[time_ch];
		for (int val_ch = 0; val_ch < usv[l].second.size(); ++val_ch)
			cout << '\t' << "Val_ch_" << val_ch << '\t' << usv[l].second[val_ch];
		cout << "\n";
	}
}

//----------------------------------------------------------------------------------------------
int MedRepository::get_dates_with_signal(int pid, vector<string> &sig_names, vector<int> &dates)
{
	int i, j;

	dates.clear();
	vector<UniversalSigVec> usvs(sig_names.size());
	for (i = 0; i < sig_names.size(); i++) {
		uget(pid, sig_names[i], usvs[i]);
		for (j = 0; j < usvs[i].len; j++)
			dates.push_back(usvs[i].Time(j));
	}

	sort(dates.begin(), dates.end());

	vector<int>::iterator it;
	it = unique(dates.begin(), dates.end());
	dates.resize(distance(dates.begin(), it));

	return 0;
}

//--------------------------------------------------------------------------------------
int MedRepository::generate_fnames_for_prefix()
{
	data_fnames.clear();
	index_fnames.clear();
	//vector<int> ii;
	//for (int i=0; i<sigs.signals_names.size(); i++) ii.push_back(i); // this is just for debugging, making sure we are order independent
	//random_shuffle(ii.begin(), ii.end());
	for (int i = 0; i < sigs.signals_names.size(); i++) {
		string fixed_sig_name = sigs.signals_names[i];
		boost::replace_all(fixed_sig_name, "/", "_div_");
		boost::replace_all(fixed_sig_name, ":", "_over_");
		boost::replace_all(fixed_sig_name, "%", "_percent_");
		string name = path + "/" + rep_files_prefix + "_" + fixed_sig_name; //sigs.signals_names[i];
		data_fnames.push_back(name + ".data");
		index_fnames.push_back(name + ".idx");
		int sid = sigs.sid(sigs.signals_names[i]);
		sigs.Sid2Info[sid].fno = i;
		//MLOG("i %d name %s sid %d\n", i, name.c_str(), sid);
	}

	return 0;
}

//--------------------------------------------------------------------------------------
// mode 3 index and data preps - looping over all sids
int MedRepository::read_index_tables(const vector<int> &pids_to_take, const vector<int> &signals_to_take)
{
	vector<int> sig_ids = signals_to_take;

	if (sig_ids.size() == 0)
		sig_ids = sigs.signals_ids;

	int nerr = 0;

	index.index_table.resize(MAX_SID_NUMBER);
	if (pids_to_take.size() == 0 || pids_to_take[0] >= 0) { // this allows for reading an empty repository by signing it with a negative value for the first pid
		//#pragma omp parallel for num_threads(2)
#pragma omp parallel for
		for (int i = 0; i < sig_ids.size(); i++) {
			int sid = sig_ids[i];
			int fno = sigs.Sid2Info[sid].fno;
			MLOG_D("Reading: sid %d index %s data %s\n", sid, index_fnames[fno].c_str(), data_fnames[fno].c_str());
			if (index.read_index_table_and_data(sid, index_fnames[fno], data_fnames[fno], pids_to_take) < 0) {
				MERR("MedRepository:read_index_tables: FAILED at sid %d index file %s\n", sid, index_fnames[fno].c_str());
				nerr++;
			}
		}
	}
	if (nerr > 0) return -1;
	index.update_pids();
	pids = index.pids;

	unsigned long long all_index_size = 0;
	unsigned long long all_data_size = 0;
	for (int i = 0; i < sig_ids.size(); i++) {
		int sid = sig_ids[i];
		all_index_size += index.index_table[sid].get_size();
		all_data_size += index.index_table[sid].w_size;
	}

	double idx_gb, data_gb, tot_gb;
	idx_gb = (double)all_index_size / (double)(1 << 30);
	data_gb = (double)all_data_size / (double)(1 << 30);
	tot_gb = idx_gb + data_gb;

	MLOG("Read %d signals, %d pids :: data %6.3fGB :: idx %6.3fGB :: tot %6.3fGB\n", index.signals.size(), index.pids.size(), data_gb, idx_gb, tot_gb);
	return 0;
}

//--------------------------------------------------------------------------------------
int MedRepository::load(const string &sig_name)
{
	int sid;
	if ((sid = sigs.sid(sig_name)) < 0)
		return -1;
	return load(sid);
}

//--------------------------------------------------------------------------------------
int MedRepository::load(const int sid)
{
	vector<int> pids;
	return load(sid, pids);
}

//--------------------------------------------------------------------------------------
int MedRepository::load(const vector<string> &sig_names)
{
	int rc = 0;
	for (auto &sname : sig_names)
		rc += load(sname);
	return rc;
}

//--------------------------------------------------------------------------------------
int MedRepository::load(const vector<int> &sids)
{
	int rc = 0;
	for (int sid : sids)
		rc += load(sid);
	return rc;
}

//--------------------------------------------------------------------------------------
int MedRepository::load(const string &sig_name, vector<int> &pids_to_take)
{
	int sid;
	if ((sid = sigs.sid(sig_name)) < 0)
		return -1;
	return load(sid, pids_to_take);
}

//--------------------------------------------------------------------------------------
int MedRepository::load(const vector<string> &sig_names, vector<int> &pids_to_take)
{
	unordered_set<string> sig_set;
	for (auto &s : sig_names) sig_set.insert(s);
	vector<int> sids;
	for (auto &s : sig_set) {
		int sid = sigs.sid(s);
		if (sid > 0)
			sids.push_back(sid);
		else
			MERR("MedRepository::load() : ERROR: signal %s is illegal\n", s.c_str());
	}

	return load(sids, pids_to_take);

}

//--------------------------------------------------------------------------------------
int MedRepository::load(const vector<int> &sids, vector<int> &pids_to_take)
{
	unordered_set<int> sid_set;
	for (auto s : sids) sid_set.insert(s);

	vector<int> pids_sort_uniq = pids_to_take;
	sort(pids_sort_uniq.begin(), pids_sort_uniq.end());
	auto it = unique(pids_sort_uniq.begin(), pids_sort_uniq.end());
	pids_sort_uniq.resize(distance(pids_sort_uniq.begin(), it));

	vector<int> usids;
	copy(sid_set.begin(), sid_set.end(), back_inserter(usids));

	int rc = 0;
#pragma omp parallel for
	for (int j = 0; j < usids.size(); j++) {
		int sid = usids[j];
		int local_rc = load_pids_sorted(sid, pids_sort_uniq);
		if (local_rc < 0) MERR("MedRepository::load() : ERROR: Failed reading sid %d (%s)\n", sid, sigs.Sid2Name[sid].c_str());
#pragma omp critical
		rc += local_rc;
	}
	return rc;

}

//--------------------------------------------------------------------------------------
// rc -2: not loaded due to locking
int MedRepository::load(const int sid, vector<int> &pids_to_take)
{
	if (sid < 0 || sid >= index.index_table.size())
		return -1;

	if (index.index_table[sid].full_load)
		return 0; // nothing to do already fully loaded

	if (index.index_table[sid].is_locked)
		return -2; // need to load but can't since it is locked

	int fno = sigs.Sid2Info[sid].fno;

	vector<int> pids_sort_uniq = pids_to_take;
	sort(pids_sort_uniq.begin(), pids_sort_uniq.end());
	auto it = unique(pids_sort_uniq.begin(), pids_sort_uniq.end());
	pids_sort_uniq.resize(distance(pids_sort_uniq.begin(), it));

	if (index.index_table[sid].read_index_and_data(index_fnames[fno], data_fnames[fno], pids_sort_uniq) < -1)
		return -1;

	return 0;
}

//--------------------------------------------------------------------------------------
// rc -2: not loaded due to locking
int MedRepository::load_pids_sorted(const int sid, vector<int> &pids_sort_uniq)
{
	if (sid < 0 || sid >= index.index_table.size())
		return -1;

	if (index.index_table[sid].full_load)
		return 0; // nothing to do already fully loaded

	if (index.index_table[sid].is_locked)
		return -2; // need to load but can't since it is locked

	int fno = sigs.Sid2Info[sid].fno;


	if (index.index_table[sid].read_index_and_data(index_fnames[fno], data_fnames[fno], pids_sort_uniq) < -1)
		return -1;

	return 0;
}


//--------------------------------------------------------------------------------------
int MedRepository::lock_all_sigs()
{
	return lock(sigs.signals_ids);
}

//--------------------------------------------------------------------------------------
int MedRepository::lock(const string &sig_name)
{
	int sid;
	if ((sid = sigs.sid(sig_name)) < 0)
		return -1;
	return lock(sid);
}

//--------------------------------------------------------------------------------------
int MedRepository::lock(const int sid)
{
	if (sid < 0 || sid >= index.index_table.size())
		return -1;

	index.index_table[sid].lock();
	return 0;
}

//--------------------------------------------------------------------------------------
int MedRepository::lock(const vector<string> &sig_names)
{
	int rc = 0;
	for (auto &sname : sig_names)
		rc += lock(sname);
	return rc;
}

//--------------------------------------------------------------------------------------
int MedRepository::lock(const vector<int> &sids)
{
	int rc = 0;
	for (int sid : sids)
		rc += lock(sid);
	return rc;
}

//--------------------------------------------------------------------------------------
int MedRepository::unlock_all_sigs()
{
	return unlock(sigs.signals_ids);
}

//--------------------------------------------------------------------------------------
int MedRepository::unlock(const string &sig_name)
{
	int sid;
	if ((sid = sigs.sid(sig_name)) < 0)
		return -1;
	return unlock(sid);
}

//--------------------------------------------------------------------------------------
int MedRepository::unlock(const int sid)
{
	if (sid < 0 || sid >= index.index_table.size())
		return -1;

	index.index_table[sid].unlock();
	return 0;
}

//--------------------------------------------------------------------------------------
int MedRepository::unlock(const vector<string> &sig_names)
{
	int rc = 0;
	for (auto &sname : sig_names)
		rc += unlock(sname);
	return rc;
}

//--------------------------------------------------------------------------------------
int MedRepository::unlock(const vector<int> &sids)
{
	int rc = 0;
	for (int sid : sids)
		rc += unlock(sid);
	return rc;
}

//--------------------------------------------------------------------------------------
int MedRepository::free_all_sigs()
{
	return this->free(sigs.signals_ids);
}

//--------------------------------------------------------------------------------------
int MedRepository::free(const string &sig_name)
{
	int sid;
	if ((sid = sigs.sid(sig_name)) < 0)
		return -1;
	return free(sid);
}

//--------------------------------------------------------------------------------------
int MedRepository::free(const int sid)
{
	if (sid < 0 || sid >= index.index_table.size())
		return -1;

	// freeing the index table matching sid
	index.index_table[sid].clear();
	return 0;
}

//--------------------------------------------------------------------------------------
int MedRepository::free(const vector<string> &sig_names)
{
	int rc = 0;

	for (auto &sname : sig_names)
		rc += free(sname);

	return rc;
}

//--------------------------------------------------------------------------------------
int MedRepository::free(const vector<int> &sids)
{
	int rc = 0;
	for (int sid : sids)
		rc += free(sid);
	return rc;
}

//--------------------------------------------------------------------------------------
int MedRepository::get_pids_with_sig(const string &sig_name, vector<int> &in_pids)
{
	if (rep_mode < 3) {
		MERR("get_pids_with_sig():: supporting mode 3 and up only\n");
		return -1;
	}

	int sid = sigs.sid(sig_name);
	in_pids.clear();
	if (index.index_table[sid].is_loaded) {
		vector<unsigned int> _in_pids;
		index.index_table[sid].sv.get_all_keys(_in_pids);
		copy(_in_pids.begin(), _in_pids.end(), back_inserter(in_pids));
	}

	return 0;
}

//--------------------------------------------------------------------------------------
int IndexTable::insert(unsigned int pid, unsigned int delta, int len)
{
	if (sv.data.size() <= 1)
		sv.set_min(pid);
	last_len = len;
	acc = delta;
	return sv.insert(pid, delta);
}

//--------------------------------------------------------------------------------------
int IndexTable::get_len(const unsigned int pid)
{
	unsigned int ind;

	ind = sv.get_ind(pid);
	if (ind == 0)
		return 0;

	if (ind < sv.data.size() - 1)
		return (int)(sv.data[ind + 1] - sv.data[ind]);

	return last_len;

}

//--------------------------------------------------------------------------------------
int IndexTable::get(const unsigned int pid, unsigned long long &pos, int &len)
{
	pos = 0;
	len = 0;
	unsigned int ind;

	ind = sv.get_ind(pid);
	if (ind == 0)
		return -1;

	pos = sv.data[ind];

	/*
	if (factor == 8)
		pos = pos << 3; // saving time, 8 is the most common factor (after that 4).
	else if (factor == 4)
		pos = pos << 2;
	else if (factor == 16)
		pos = (pos << 4);
	else if (factor == 12)
		pos = (pos << 3) + (pos << 2);
	else if (factor == 20)
		pos = (pos << 4) + (pos << 2);
	else*/
	pos *= factor;

	pos += base;

	if (ind < sv.data.size() - 1)
		len = sv.data[ind + 1] - sv.data[ind];
	else
		len = last_len;

	//MLOG("IndexTable get: sid %d pid %d pos %lld base %lld factor %d ind %d len %d\n", sid, pid, pos, base, factor, ind, len);

	return 0;
}

//--------------------------------------------------------------------------------------
unsigned long long IndexTable::get_data_size()
{
	lock_guard<mutex> guard(index_table_locks[sid]);

	if (sv.data.size() < 2)
		return 0;

	unsigned long long size = sv.data.back() - sv.data[1];
	size += last_len;
	size *= factor;
	return size;
}

//--------------------------------------------------------------------------------------
size_t IndexTable::get_size()
{

	lock_guard<mutex> guard(index_table_locks[sid]);
	size_t size = 0;

	size += sizeof(int); // last_len
	size += sizeof(int); // sid
	size += sizeof(unsigned long long); // base
	size += sizeof(unsigned int); // factor
	size += sv.get_size(); // sparse vec

	return size;
}

//--------------------------------------------------------------------------------------
size_t IndexTable::serialize(unsigned char *blob)
{

	lock_guard<mutex> guard(index_table_locks[sid]);
	unsigned char *curr = blob;
	size_t len = 0;

	((int *)curr)[0] = last_len; curr += sizeof(int);		// last_len
	((int *)curr)[0] = sid; curr += sizeof(int);				// sid
	((unsigned long long *)curr)[0] = base; curr += sizeof(unsigned long long);				// base
	((unsigned int *)curr)[0] = factor; curr += sizeof(unsigned int);				// factor
	len += curr - blob;
	len += sv.serialize(curr);

	return len;
}

//--------------------------------------------------------------------------------------
size_t IndexTable::deserialize(unsigned char *blob)
{
	lock_guard<mutex> guard(index_table_locks[sid]);
	unsigned char *curr = blob;
	size_t len = 0;

	last_len = ((int *)curr)[0]; curr += sizeof(int);
	sid = ((int *)curr)[0]; curr += sizeof(int);
	base = ((unsigned long long *)curr)[0]; curr += sizeof(unsigned long long);
	factor = ((unsigned int *)curr)[0]; curr += sizeof(unsigned int);
	len += curr - blob;
	len += sv.deserialize(curr);

	return len;
}

//--------------------------------------------------------------------------------------
int IndexTable::write_to_file(string &fname)
{

	// get size & allocate
	size_t size = get_size();
	unsigned char *data = new unsigned char[size];

	// serialize
	size_t ser_size;
	ser_size = serialize(data);
	MLOG("Serializing index for file %s : size %lld (%lld): sid %d : base %lld : factor %d : acc %d : last_len %d\n",
		fname.c_str(), ser_size, size, sid, base, factor, acc, last_len);

	// write
	int rc;
	{
		lock_guard<mutex> guard(index_table_locks[sid]);
		rc = write_bin_file_IM(fname, data, size);
	}

	// dealloc
	delete[] data;

	return rc;
}

//--------------------------------------------------------------------------------------
int IndexTable::read_from_file(string &fname)
{

	unsigned char *data = NULL;
	unsigned long long size = 0;
	int rc;

	{
		lock_guard<mutex> guard(index_table_locks[sid]);
		//		rc = read_bin_file_IM_parallel(fname, data, size);
		//		MLOG("IndexTable::read_from_file fname is %s\n", fname.c_str());
		rc = read_bin_file_IM(fname, data, size);
		if (rc != 0 && file_exists_IM(fname))
			MTHROW_AND_ERR("IndexTable::read_from_file could not read from [%s]\n", fname.c_str());
	}
	if (size > 0) {
		deserialize(data);
		//		MLOG("IndexTable::read_from_file fname [%s] rc %d before delete[] data %d\n", fname.c_str(), rc, data);
		if (data != NULL)
			delete[] data;
		//		MLOG("IndexTable::read_from_file fname [%s] rc %d after delete[] data %d\n", fname.c_str(), rc, data);
	}

	//if (rc < 0 && !file_exists_IM(fname)) rc = 0;

	return rc;
}


//--------------------------------------------------------------------------------------
void IndexTable::clear()
{
	if (is_locked || sid == 0)
		return;

	lock_guard<mutex> guard(index_table_locks[sid]);

	if (work_area_allocated) {
		if (work_area != NULL) {
			//MLOG("IndexTable::clear before work_area [%d]\n", work_area);
			delete[] work_area;
			//MLOG("IndexTable::clear after work_area [%d]\n", work_area);
		}
	}
	work_area = NULL;
	init();
	try {
		sv.clear();
	}
	catch (...) {
		MERR("Error IndexTable::clear - not enough memory\n");
	}

}

//--------------------------------------------------------------------------------------
void IndexTable::lock()
{
	lock_guard<mutex> guard(index_table_locks[sid]);
	is_locked = 1;
}

//--------------------------------------------------------------------------------------
void IndexTable::unlock()
{
	lock_guard<mutex> guard(index_table_locks[sid]);
	is_locked = 0;
}

//--------------------------------------------------------------------------------------
int IndexTable::read_index_and_data(string &idx_fname, string &data_fname)
{
	vector<int> take_pids;
	return read_index_and_data(idx_fname, data_fname, take_pids);
}

//--------------------------------------------------------------------------------------
int IndexTable::read_index_and_data(string &idx_fname, string &data_fname, const vector<int> &pids_to_include)
{
	string prefix = "IndexTable::read_index_and_data :: ";

	// threads locking: making sure only a single thread at a time can load this specific sid
	// using the index_table_read_locks set of locks, since the inner functions use the index_table_locks mechanism
	// MLOG("in read_index_and_data with sid %d\n", sid);
	lock_guard<mutex> guard(index_table_read_locks[sid]);

	if (is_loaded && is_locked) {
		MTHROW_AND_ERR("%s ERROR: failed reading index table %s since this signal is LOCKED\n", prefix.c_str(), idx_fname.c_str());
		return -2;
	}

	int keep_lock = is_locked;
	clear(); // making sure we free all previous loads
	is_locked = keep_lock;

	if (pids_to_include.size() == 0) {

		// full mode size

		//MLOG("in here:: idx_fname %s data_fname %s\n", idx_fname.c_str(), data_fname.c_str());
		// read index table
		if (read_from_file(idx_fname) < 0) {
			MTHROW_AND_ERR("%s ERROR: failed reading index table %s\n", prefix.c_str(), idx_fname.c_str());
			return -1;
		}

		if (file_exists_IM(data_fname)) {
			//if (read_bin_file_IM_parallel(data_fname, work_area, w_size) < 0) {
			if (read_bin_file_IM(data_fname, work_area, w_size) < 0) {
				MTHROW_AND_ERR("%s ERROR: can't read or allocate for file %s\n", prefix.c_str(), data_fname.c_str());
				return -1;
			}
		}
		else {
			MERR("%s ERROR: file %s not found\n", prefix.c_str(), data_fname.c_str());
			w_size = 0;
			work_area = NULL;
		}

		if (w_size > 0) work_area_allocated = 1;
		is_loaded = 1;
		full_load = 1;
		//MLOG("sid %d : %s index size %ld data_size %ld %ld\n", sid, idx_fname.c_str(), get_size(), d_size, w_size);

	}
	else {
		// partial mode - read only selected pids

		// this is more complex, we need to open the data file and read only the pids we need
		// while creating an index for them

		// first we read the full index for the sid

		//MLOG("Reading Index Table for %d pids from %s (sid %d factor %d) (thread %d)\n", pids_to_include.size(), idx_fname.c_str(), sid, factor, omp_get_thread_num());

		IndexTable idx;
		if (idx.read_from_file(idx_fname) < 0) {
			MTHROW_AND_ERR("%s ERROR: failed reading index table %s\n", prefix.c_str(), idx_fname.c_str());
			return -1;
		}

		// from here on we assume pids_to_include is sorted and uniq.
		// this will allow us a faster code
		unsigned long long d_size = 0;

		vector<int> keys;
		vector<int> inds;
		vector<int> lens;
		idx.sv.get_all_intersected_keys(pids_to_include, keys, inds);
		lens.resize(inds.size());

		if (inds.size() > 0) {
			int ii = 0;
			for (ii = 0; ii < inds.size() - 1; ii++) {
				d_size += (lens[ii] = idx.sv.data[inds[ii] + 1] - idx.sv.data[inds[ii]]);
			}
			if (inds[ii] < idx.sv.data.size() - 1)
				d_size += (lens[ii] = idx.sv.data[inds[ii] + 1] - idx.sv.data[inds[ii]]);
			else
				d_size += (lens[ii] = idx.last_len);
		}

		d_size *= idx.factor;
		//MLOG("sid = %d , d_size = %d\n", sid, d_size);



		// initializing our index_table 
		base = 0;
		factor = idx.factor;
		acc = 0;
		last_len = 0;
		sid = idx.sid;

		if (d_size == 0)
			return 0;

		//MLOG("before allocating work_area d_size %d\n", d_size);
		work_area = new unsigned char[d_size];
		if (work_area == NULL)
			MTHROW_AND_ERR("IndexTable::read_index_and_data could not allocate work_area\n")
			w_size = d_size;
		work_area_allocated = 1;

		// now reading with a buffer one by one , copying and updating index_table
		MedBufferedFile mbf;

		if (mbf.open(data_fname) < 0) {
			MTHROW_AND_ERR("ERROR: Can't open data file %s\n", data_fname.c_str());
			return -1;
		}

		unsigned char *buf = work_area;
		//MLOG("pid_to_read : %d pids factor %d sid %d\n", pids_to_include.size(),factor,sid);

		for (int i = 0; i < inds.size(); i++) {
			unsigned long long pos = idx.base + (unsigned long long)idx.sv.data[inds[i]] * (unsigned long long)idx.factor;
			int blen = lens[i] * factor;
			//MLOG("before pid %d sid %d pos %lld factor %d %d blen %d i %d idx.base %d\n", inds[i], sid, pos, factor, idx.factor, blen, i, idx.base);
			mbf.read(buf, pos, blen);
			//MLOG("after pid %d sid %d pos %lld factor %d %d blen %d i %d idx.base %d\n", inds[i], sid, pos, factor, idx.factor, blen, i, idx.base);
			buf += blen;
			if (insert(keys[i], lens[i]) < 0)
				MTHROW_AND_ERR("IndexTable::read_index_and_data could not insert(keys[i], lens[i])\n");
		}

		mbf.close();
		is_loaded = 1;
		full_load = 0;
	}
	tot_size = get_size() + get_data_size();
	tot_size_gb = (double)tot_size / (double)(1 << 30);

	return 0;
}

//======================================================================================================================
int UsvsIterator::init(MedRepository *_rep, const vector<string> &_sig_names, const vector<UniversalSigVec *> &_usvs)
{
	rep = _rep;
	sig_names = _sig_names;
	sids.clear();
	for (auto sig : sig_names) {
		sids.push_back(rep->sigs.sid(sig));
		if (sids.back() <= 0) {
			MERR("ERROR: no susch sig %s in UsvsIterator init()\n", sig.c_str());
			return -1;
		}
	}
	usvs = _usvs;
	if (usvs.size() != sids.size()) {
		MERR("ERROR: non matching lengths: %d usvs, %d sigs, in UsvsIterator init()\n", usvs.size(), sids.size());
		return -1;
	}
	return 0;
}

//======================================================================================================================
int UsvsIterator::read_pid(int pid)
{
	for (int i = 0; i < usvs.size(); i++)
		rep->uget(pid, sids[i], *usvs[i]);
	return 0;
}

//======================================================================================================================
int UsvsIterator::read_pid(int pid, const vector<UniversalSigVec *> &_usvs)
{
	if (usvs.size() != sids.size()) {
		MERR("ERROR: non matching lengths: %d usvs, %d sigs, in UsvsIterator read_pid()\n", _usvs.size(), sids.size());
		return -1;
	}
	for (int i = 0; i < _usvs.size(); i++)
		rep->uget(pid, sids[i], *_usvs[i]);
	return 0;
}


float medial::repository::DateDiff(int refDate, int dateSample) {
	return float((med_time_converter.convert_times(global_default_time_unit, MedTime::Days, dateSample) -
		med_time_converter.convert_times(global_default_time_unit, MedTime::Days, refDate)) / 365.0);
}

int medial::repository::DateAdd(int refDate, int daysAdd) {
	if (global_default_time_unit == MedTime::Date) {
		int days_val = med_time_converter.convert_times(global_default_time_unit, MedTime::Days, refDate) + daysAdd;
		if (days_val > 0)
			return med_time_converter.convert_days(global_default_time_unit, days_val);
		else
			return 0;
	}
	else
		return refDate + daysAdd;
}

int medial::repository::get_value(MedRepository &rep, int pid, int sigCode) {
	int gend = -1;
	UniversalSigVec usv;

	rep.uget(pid, sigCode, usv);

	if (usv.len > 0) {
		if (usv.n_val_channels() > 0)
			gend = (int)usv.Val(0);
		else if (usv.n_time_channels() > 0)
			gend = usv.Time(0);
	}
	return gend;
}

int medial::repository::get_value(PidRec &rep, int sigCode) {
	int gend = -1;
	UniversalSigVec usv;
	rep.uget(sigCode, usv);

	if (usv.len > 0) {
		if (usv.n_val_channels() > 0)
			gend = (int)usv.Val(0);
		else if (usv.n_time_channels() > 0)
			gend = usv.Time(0);
	}
	return gend;
}


template<class T>int medial::repository::fetch_next_date(vector<T> &patientFile, vector<int> &signalPointers) {
	int minDate = -1, minDate_index = -1;
	for (size_t i = 0; i < patientFile.size(); ++i)
	{
		UniversalSigVec &data = patientFile[i];
		if (signalPointers[i] >= data.len)
			continue; //already reached the end for this signal
		if (data.get_type() == T_Value || data.n_time_channels() == 0) {
			if (minDate_index == -1 || data.Val(signalPointers[i]) < minDate) {
				minDate = (int)data.Val(signalPointers[i]);
				minDate_index = (int)i;
			}
		}
		else {
			if (minDate_index == -1 || data.Time(signalPointers[i]) < minDate) {
				minDate = data.Time(signalPointers[i]);
				minDate_index = (int)i;
			}
		}
	}
	if (minDate_index >= 0)
		++signalPointers[minDate_index];
	return minDate_index;
}
template int medial::repository::fetch_next_date<UniversalSigVec_mem>(vector<UniversalSigVec_mem> &patientFile, vector<int> &signalPointers);
template int medial::repository::fetch_next_date<UniversalSigVec>(vector<UniversalSigVec> &patientFile, vector<int> &signalPointers);


string medial::signal_hierarchy::filter_code_hierarchy(const vector<string> &vec, const string &signalHirerchyType) {
	if (vec.empty())
		return "";
	return vec.front(); //always first is the coded
}
vector<int> medial::signal_hierarchy::parents_code_hierarchy(MedDictionarySections &dict,
	const string &group, const string &signalHirerchyType, int depth, int max_nodes) {
	int sectionId = 0;
	if (dict.SectionName2Id.find(signalHirerchyType) == dict.SectionName2Id.end())
		MTHROW_AND_ERR("Signal_Hirerchy_Type not suppoted %s. please select dictionary section: \n",
			signalHirerchyType.c_str());
	sectionId = dict.SectionName2Id.at(signalHirerchyType);

	vector<int> parents;
	vector<int> last_parents = { dict.id(sectionId, group) };
	if (last_parents.front() < 0)
		return parents; //no parents

	for (size_t k = 0; k < depth; ++k) {
		vector<int> tmp_par, new_layer;
		for (int par : last_parents)
		{
			dict.get_member_sets(sectionId, par, tmp_par);
			new_layer.insert(new_layer.end(), tmp_par.begin(), tmp_par.end());
			parents.insert(parents.end(), tmp_par.begin(), tmp_par.end()); //aggregate all parents
		}
		if (parents.size() >= max_nodes)
			break;
		new_layer.swap(last_parents);
		if (last_parents.empty())
			break; //no more parents to loop up
	}

	return parents;
}
vector<int> medial::signal_hierarchy::sons_code_hierarchy(MedDictionarySections &dict, const string &group, const string &signalHirerchyType) {
	int sectionId = 0;
	if (dict.SectionName2Id.find(signalHirerchyType) == dict.SectionName2Id.end())
		MTHROW_AND_ERR("Signal_Hirerchy_Type not suppoted %s. please select dictionary section: \n",
			signalHirerchyType.c_str());
	sectionId = dict.SectionName2Id.at(signalHirerchyType);
	vector<int> sons;
	dict.get_set_members(sectionId, group, sons);
	return sons;
}

void medial::signal_hierarchy::sons_code_hierarchy_recursive(MedDictionarySections &dict,
	const string &signalName, const string &code, vector<int> &flat_all_sons, int max_depth) {
	int sectionId = dict.section_id(signalName);
	//includes the code itself!

	int code_id = dict.id(sectionId, code);
	if (code_id < 0)
		MTHROW_AND_ERR("Error in medial::signal_hierarchy::sons_code_hierarchy_recursive - can't find code %s in signal %s\n",
			code.c_str(), signalName.c_str());

	vector<int> last_layer_new_ids;
	unordered_set<int> full_list;
	full_list.insert(code_id);
	last_layer_new_ids.push_back(code_id);
	int depth_cnt = 0;
	while (!last_layer_new_ids.empty()) {
		vector<int> new_ids;
		for (int id_in : last_layer_new_ids)
		{
			vector<int> fetch_ids;
			dict.get_set_members(sectionId, id_in, fetch_ids);
			for (int fetch_id : fetch_ids)
			{
				//not seen id - new one
				if (full_list.find(fetch_id) == full_list.end()) {
					full_list.insert(fetch_id);
					new_ids.push_back(fetch_id);
				}
			}
		}
		last_layer_new_ids = move(new_ids);
		++depth_cnt;
		if (max_depth > 0 && depth_cnt >= max_depth)
			break;
	}

	//return all ids in full_list
	flat_all_sons.clear();
	flat_all_sons.insert(flat_all_sons.end(), full_list.begin(), full_list.end());
}

string medial::signal_hierarchy::get_readcode_code(MedDictionarySections &dict, int id, const string &signalHirerchyType) {
	int sectionId = 0;
	if (dict.SectionName2Id.find(signalHirerchyType) == dict.SectionName2Id.end())
		MTHROW_AND_ERR("Signal_Hirerchy_Type not suppoted %s. please select dictionary section: \n",
			signalHirerchyType.c_str());
	sectionId = dict.SectionName2Id.at(signalHirerchyType);

	string res = "";
	vector<int> sets;
	dict.get_member_sets(sectionId, id, sets);
	if (dict.dicts[sectionId].Id2Names.find(id) != dict.dicts[sectionId].Id2Names.end()) {
		vector<string> nms = dict.dicts[sectionId].Id2Names.at(id);
		string filterd = filter_code_hierarchy(nms, signalHirerchyType);
		if (!filterd.empty()) {
			res = filterd;
			return res;
		}
	}
	for (int s : sets)
	{
		vector<string> names;
		if (dict.dicts[sectionId].Id2Names.find(s) != dict.dicts[sectionId].Id2Names.end()) {
			names = dict.dicts[sectionId].Id2Names.at(s);
		}
		string name = filter_code_hierarchy(names, signalHirerchyType);
		if (!name.empty() && !res.empty()) {
			//more than one option - for now not supported:
			MTHROW_AND_ERR("not supported code %d for signal=%s\n",
				id, signalHirerchyType.c_str());
		}
		if (!name.empty())
			res = name;
	}

	return res;
}

bool medial::repository::fix_contradictions(UniversalSigVec &s, fix_method method,
	UniversalSigVec_mem &edited) {
	edited.set(s);
	edited.manage = false;
	if (s.n_time_channels() == 0 || s.get_type() == T_TimeStamp || method == fix_method::none)
		return false;

	vector<vector<char *>> data_by_dates;
	char *p = (char *)s.data;
	int st_size = (int)s.size();
	if (s.len > 0)
		data_by_dates.push_back({ p });
	p += st_size;
	for (int i = 1; i < s.len; ++i)
	{
		if (s.Time(i) == s.Time(i - 1) && (data_by_dates.back().size() > 1 || s.Val(i) != s.Val(i - 1)))
			data_by_dates.back().push_back(p);
		else
			data_by_dates.push_back({ p });
		p += st_size;
	}

	bool changed = false;
	vector<int> to_remove;
	int curr_pos = 0;
	float mean_val = 0;
	for (size_t i = 0; i < data_by_dates.size(); ++i)
	{
		vector<char *> &data_group = data_by_dates[i];
		if (data_group.size() == 1) {
			curr_pos += (int)data_group.size();
			continue; //all ok - do nothing
		}

		changed = true;
		switch (method)
		{
		case medial::repository::drop:
			for (int k = 0; k < data_group.size(); ++k)
				to_remove.push_back(curr_pos + k);
			break;
		case medial::repository::take_first:
			for (int k = 1; k < data_group.size(); ++k)
				to_remove.push_back(curr_pos + k);
			break;
		case medial::repository::take_last:
			for (int k = 0; k < (int)data_group.size() - 1; ++k)
				to_remove.push_back(curr_pos + k);
			break;
		case medial::repository::take_mean:
			for (int k = 1; k < data_group.size(); ++k)
				to_remove.push_back(curr_pos + k);
			//change first to be mean value:
			mean_val = 0;
			for (int k = 0; k < data_group.size(); ++k)
				mean_val += s.Val(curr_pos + k);
			mean_val /= data_group.size();
			s.SetVal_ch_vec(curr_pos, 0, mean_val, s.data);
			break;
		case medial::repository::take_min:
			//lets find min;
			mean_val = s.Val(curr_pos);
			for (int k = 1; k < data_group.size(); ++k) {
				if (s.Val(curr_pos + k) < mean_val)
					mean_val = s.Val(curr_pos + k);
				to_remove.push_back(curr_pos + k);
			}
			s.SetVal_ch_vec(curr_pos, 0, mean_val, s.data);
			break;
		case medial::repository::take_max:
			//lets find min;
			mean_val = s.Val(curr_pos);
			for (int k = 1; k < data_group.size(); ++k) {
				if (s.Val(curr_pos + k) > mean_val)
					mean_val = s.Val(curr_pos + k);
				to_remove.push_back(curr_pos + k);
			}
			s.SetVal_ch_vec(curr_pos, 0, mean_val, s.data);
			break;
		default:
			MTHROW_AND_ERR("Error in fix_contradictions - Not Implemented\n");
		}

		curr_pos += (int)data_group.size();
	}

	if (!to_remove.empty()) {
		changed = true;
		int rem_idx = 0, curr_i = 0;
		p = (char *)s.data;
		char *after_filter = new char[(s.len - (int)to_remove.size()) * st_size];
		for (size_t i = 0; i < s.len; ++i)
		{
			if (rem_idx < to_remove.size() && i == to_remove[rem_idx]) {
				++rem_idx;
				continue;
				//remove current
			}
			memcpy(after_filter + curr_i * st_size, p + i * st_size, st_size);
			++curr_i;
		}

		s.data = after_filter;
		s.len -= (int)to_remove.size();

		edited.set(s);
		edited.manage = true;
	}
	return changed;
}

void medial::repository::set_global_time_unit(const string &repository_path) {
	MedRepository temp_rep;
	if (temp_rep.read_config(repository_path) < 0)
		MTHROW_AND_ERR("Can't read repository %s\n", repository_path.c_str());
	global_default_time_unit = temp_rep.time_unit;
	global_default_windows_time_unit = MedTime::Days;
	if (global_default_time_unit != MedTime::Date) {
		global_default_windows_time_unit = MedTime::Minutes;
		MLOG("model_runner: running on ICU repository\n");
	}
}