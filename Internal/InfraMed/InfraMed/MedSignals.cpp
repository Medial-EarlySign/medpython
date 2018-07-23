//
// MedSignals.c
//

#ifndef _SCL_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS
#endif

#include "MedSignals.h"
#include "InfraMed.h"
#include "Logger/Logger/Logger.h"
#include <fstream>
#include <algorithm>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <mutex>

mutex insert_signal_mutex;

#define LOCAL_SECTION LOG_SIG
#define LOCAL_LEVEL LOG_DEF_LEVEL
extern MedLogger global_logger;

using namespace std;
using namespace boost;

//-----------------------------------------------------------------------------------------------
int MedRep::get_type_size(SigType t)
{
	switch (t) {

	case T_Value:
		return ((int)sizeof(SVal));

	case T_DateVal:
		return ((int)sizeof(SDateVal));

	case T_TimeVal:
		return ((int)sizeof(STimeVal));

	case T_DateRangeVal:
		return ((int)sizeof(SDateRangeVal));

	case T_TimeRangeVal:
		return ((int) sizeof(STimeRangeVal));

	case T_TimeStamp:
		return ((int) sizeof(STimeStamp));

	case T_DateVal2:
		return ((int)sizeof(SDateVal2));

	case T_TimeLongVal:
		return ((int)sizeof(STimeLongVal));

	case T_DateShort2:
		return ((int)sizeof(SDateShort2));

	case T_ValShort2:
		return ((int)sizeof(SValShort2));

	case T_ValShort4:
		return ((int)sizeof(SValShort4));

	case T_CompactDateVal:
		return ((int)sizeof(SCompactDateVal));

	default:
		MERR("Cannot get size of signal type %d\n", t);
		return 0;
	}

	return 0;
}

//-----------------------------------------------------------------------------------------------
int MedRep::get_type_channels(SigType t, int &time_unit, int &n_time_chans, int &n_val_chans)
{
	switch (t) {

	case T_Value:
		return MedRep::get_type_channels_info<SVal>(time_unit, n_time_chans, n_val_chans);

	case T_DateVal:
		return MedRep::get_type_channels_info<SDateVal>(time_unit, n_time_chans, n_val_chans);

	case T_TimeVal:
		return MedRep::get_type_channels_info<STimeVal>(time_unit, n_time_chans, n_val_chans);

	case T_DateRangeVal:
		return MedRep::get_type_channels_info<SDateRangeVal>(time_unit, n_time_chans, n_val_chans);

	case T_TimeRangeVal:
		return MedRep::get_type_channels_info<STimeRangeVal>(time_unit, n_time_chans, n_val_chans);

	case T_TimeStamp:
		return MedRep::get_type_channels_info<STimeStamp>(time_unit, n_time_chans, n_val_chans);

	case T_DateVal2:
		return MedRep::get_type_channels_info<SDateVal2>(time_unit, n_time_chans, n_val_chans);

	case T_TimeLongVal:
		return MedRep::get_type_channels_info<STimeLongVal>(time_unit, n_time_chans, n_val_chans);

	case T_DateShort2:
		return MedRep::get_type_channels_info<SDateShort2>(time_unit, n_time_chans, n_val_chans);

	case T_ValShort2:
		return MedRep::get_type_channels_info<SValShort2>(time_unit, n_time_chans, n_val_chans);

	case T_ValShort4:
		return MedRep::get_type_channels_info<SValShort4>(time_unit, n_time_chans, n_val_chans);

	case T_CompactDateVal:
		return MedRep::get_type_channels_info<SCompactDateVal>(time_unit, n_time_chans, n_val_chans);

	default:
		MERR("Cannot get channels for signal type %d\n", t);
		return 0;
	}

	return 0;
}


//-----------------------------------------------------------------------------------------------
int MedSignals::read(vector<string> &sfnames)
{
	int rc = 0;

	for (int i = 0; i < sfnames.size(); i++) {
		rc += read(sfnames[i]);
	}

	return rc;
}

//-----------------------------------------------------------------------------------------------
int MedSignals::read(string path, vector<string> &sfnames)
{
	int rc = 0;
	for (int i = 0; i < sfnames.size(); i++) {
		string fname = (path == "") ? sfnames[i] : path + "/" + sfnames[i];
		try {
			rc += read(fname);
		}
		catch (...) {
			MERR("Error in reading Signal %s\n", fname.c_str());
			throw;
		}
	}
	return rc;
}

//-----------------------------------------------------------------------------------------------
int MedSignals::read(const string &fname)
{
	lock_guard<mutex> guard(insert_signal_mutex);

	ifstream inf(fname);

	if (!inf) {
		MERR("MedSignals: read: Can't open file %s\n", fname.c_str());
		return -1;
	}

	fnames.push_back(fname); // TBD : check that we didn't already load this file
	string curr_line;
	MLOG_D("Working on signals file %s\n", fname.c_str());
	while (getline(inf, curr_line)) {
		if ((curr_line.size() > 1) && (curr_line[0] != '#')) {
			vector<string> fields;
			split(fields, curr_line, boost::is_any_of("\t"));
			MLOG_D("MedSignals: read: file %s line %s\n", fname.c_str(), curr_line.c_str());
			if (fields.size() >= 4 && fields.size() <= 5 && fields[0].compare(0, 6, "SIGNAL") == 0) {
				// line format: SIGNAL <name> <signal id> <signal type num> <description>

				int sid = stoi(fields[2]);
				if (Name2Sid.find(fields[1]) != Name2Sid.end() || Sid2Name.find(sid) != Sid2Name.end()) {
					MERR("MedSignals: read: name or id already used in line %s\n", curr_line.c_str());
					return -1;
				}
				else {
					Name2Sid[fields[1]] = sid;
					Sid2Name[sid] = fields[1];
					signals_names.push_back(fields[1]);
					signals_ids.push_back(sid);
					int type = stoi(fields[3]);
					if (type<0 || type>(int)T_Last) {
						MERR("MedSignals: read: type %d not recognized in line %s\n", curr_line.c_str());
						return -1;
					}
					SignalInfo info;
					info.sid = sid;
					info.name = fields[1];
					info.type = type;
					info.bytes_len = MedRep::get_type_size((SigType)type);
					if (fields.size() == 4)
						info.description = "";
					else
						info.description = fields[4];
					// default time_units and channels ATM, time_unit may be optional as a parameter in the sig file in the future.
					MedRep::get_type_channels((SigType)type, info.time_unit, info.n_time_channels, info.n_val_channels);
					_Sid2Info[sid] = info;
					if (sid >= Sid2Info.size()) {
						SignalInfo si;
						si.sid = -1;
						Sid2Info.resize(sid + 1, si);
					}
					if (my_repo != NULL)
						info.time_unit = my_repo->time_unit;
					Sid2Info[sid] = info;
				}

			}
			else {
				MWARN("MedSignals: read: can't parse line: %s (%d)\n", curr_line.c_str(), fields.size());
			}

		}
	}

	sort(signals_ids.begin(), signals_ids.end());
	int max_sid = *max_element(signals_ids.begin(), signals_ids.end());
	sid2serial.clear();
	sid2serial.resize(max_sid + 1, -1);
	for (int i = 0; i < signals_ids.size(); i++)
		sid2serial[signals_ids[i]] = i;


	inf.close();
	fnames.push_back(fname);
	MLOG_D("Finished reading signals file %s\n", fname.c_str());
	return 0;
}

//-----------------------------------------------------------------------------------------------
int MedSignals::read_sfile(const string &fname)
{
	ifstream inf(fname);

	if (!inf) {
		MERR("MedSignals: read: Can't open file %s\n", fname.c_str());
		return -1;
	}

	fnames.push_back(fname); // TBD : check that we didn't already load this file
	string curr_line;
	MLOG_D("Working on signals file %s\n", fname.c_str());
	while (getline(inf, curr_line)) {
		if ((curr_line.size() > 1) && (curr_line[0] != '#')) {
			vector<string> fields;
			split(fields, curr_line, boost::is_any_of(" \t"));

			if (fields.size() == 2) {
				if (Name2Sid.find(fields[1]) != Name2Sid.end()) {
					int sid = Name2Sid[fields[1]];
					_Sid2Info[sid].fno = stoi(fields[0]);
					if (sid >= Sid2Info.size()) {
						SignalInfo si;
						si.sid = -1;
						Sid2Info.resize(sid + 1, si);
					}
					Sid2Info[sid].fno = stoi(fields[0]);
					Sid2Info[sid].sid = sid;
					if (my_repo != NULL)
						Sid2Info[sid].time_unit = my_repo->time_unit;
					//MLOG("read_sfile: %s\n", curr_line.c_str());
					//MLOG("read_sfile: sid %d fno %d\n", sid, _Sid2Info[sid].fno);
				}
				else {
					MERR("MedSignals:: ERROR: Unrecognized signal %s in file %s\n", fields[1].c_str(), fname.c_str());
				}
			}
		}
	}

	signals_to_files = fname;
	return 0;
}

//-----------------------------------------------------------------------------------------------
int MedSignals::read_sfile(const string &path, const string &sfile_name)
{
	string fname;

	if (path == "")
		fname = sfile_name;
	else
		fname = path + "/" + sfile_name;

	return(read_sfile(fname));
}

//-----------------------------------------------------------------------------------------------
string MedSignals::name(int sid)
{
	if (Sid2Name.find(sid) == Sid2Name.end())
		return string("");
	return Sid2Name[sid];
}

//-----------------------------------------------------------------------------------------------
int MedSignals::type(const string &name)
{
	if (Name2Sid.find(name) == Name2Sid.end())
		return -1;
	return type(Name2Sid[name]);
}

//-----------------------------------------------------------------------------------------------
int MedSignals::type(int sid)
{
	if (sid >= Sid2Info.size() || Sid2Info[sid].sid == -1)
		return -1;
	return Sid2Info[sid].type;
}

//-----------------------------------------------------------------------------------------------
string MedSignals::desc(const string &name)
{
	if (Name2Sid.find(name) == Name2Sid.end())
		return string("");
	return desc(Name2Sid[name]);
}

//-----------------------------------------------------------------------------------------------
string MedSignals::desc(int sid)
{
	if (sid >= Sid2Info.size() || Sid2Info[sid].sid == -1)
		return string("");
	return Sid2Info[sid].description;
}

//-----------------------------------------------------------------------------------------------
int MedSignals::fno(const string &sig_name)
{
	if (Name2Sid.find(sig_name) == Name2Sid.end())
		return -1;
	return fno(Name2Sid[sig_name]);
}


//-----------------------------------------------------------------------------------------------
int MedSignals::fno(int sid)
{
	if (sid >= Sid2Info.size() || Sid2Info[sid].sid == -1)
		return -1;
	return Sid2Info[sid].fno;
}


//-----------------------------------------------------------------------------------------------
int MedSignals::insert_virtual_signal(const string &sig_name, int type)
{
	// lock to allow concurrency
	lock_guard<mutex> guard(insert_signal_mutex);

	if (Name2Sid.find(sig_name) != Name2Sid.end()) {
		MERR("MedSignals: ERROR: Can't insert %s as virtual , it already exists\n", sig_name.c_str());
		return -1;
	}

	if (type<0 || type>(int)T_Last) {
		MERR("MedSignals: ERROR: Can't insert virtual signal %s : type %d not recognized\n", sig_name.c_str(), type);
		return -1;
	}

	// get_max_sid used currently, we will enter our sig_name as this + 1
	int max_sid = Sid2Name.rbegin()->first;
	int new_sid = max_sid + 1;

	// take care of all basic tables: Name2Sid , Sid2Name, signal_names, signal_ids, Sid2Info
	Name2Sid[sig_name] = new_sid;
	Sid2Name[new_sid] = sig_name;
	signals_names.push_back(sig_name);
	signals_ids.push_back(new_sid);


	SignalInfo info;
	info.sid = new_sid;
	info.name = sig_name;
	info.type = type;
	info.bytes_len = MedRep::get_type_size((SigType)type);
	info.description = "Virtual Signal";
	info.virtual_sig = 1;
	if (my_repo != NULL)
		info.time_unit = my_repo->time_unit;
	// default time_units and channels ATM, time_unit may be optional as a parameter in the sig file in the future.
	MedRep::get_type_channels((SigType)type, info.time_unit, info.n_time_channels, info.n_val_channels);
	if (new_sid >= Sid2Info.size()) {
		SignalInfo si;
		si.sid = -1;
		Sid2Info.resize(new_sid + 1, si);
	}
	_Sid2Info[new_sid] = info;
	Sid2Info[new_sid] = info;

	// take care of sid2serial
	sid2serial.resize(new_sid + 1, -1); // resize to include current new_sid
	sid2serial[new_sid] = (int)signals_ids.size() - 1; // -1 since serials start at 0


	return new_sid; // returning the new sid (always positive) as the rc
}


//================================================================================================
// UniversalSigVec
//================================================================================================
void UniversalSigVec::init(SigType _type)
{
	if (_type == type) return; // no need to init, same type as initiated

	type = _type;
	switch (_type) {
	case T_Value: set_funcs<SVal>(); return;
	case T_DateVal: set_funcs<SDateVal>(); return;
	case T_TimeVal: set_funcs<STimeVal>(); return;
	case T_DateRangeVal: set_funcs<SDateRangeVal>(); return;
	case T_TimeStamp: set_funcs<STimeStamp>(); return;
	case T_TimeRangeVal: set_funcs<STimeRangeVal>(); return;
	case T_DateVal2: set_funcs<SDateVal2>(); return;
		//case T_TimeLongVal: set_funcs<STimeLongVal>(); return;
	case T_DateShort2: set_funcs<SDateShort2>(); return;
	case T_ValShort2: set_funcs<SValShort2>(); return;
	case T_ValShort4: set_funcs<SValShort4>(); return;
		//case T_CompactDateVal: set_funcs<SCompactDateVal>(); return;
	}

	type = T_Last;

}


//-------------------------------------------------------------------------------------------------------------------
int MedSignalsSingleElemFill(int sig_type, char *buf, int *time_data, float *val_data)
{
	switch ((SigType)sig_type) {

	case T_Value:				SetSignalElement<SVal>(buf, time_data, val_data);				break;
	case T_DateVal:				SetSignalElement<SDateVal>(buf, time_data, val_data);			break;
	case T_TimeVal:				SetSignalElement<STimeVal>(buf, time_data, val_data);			break;
	case T_DateRangeVal:		SetSignalElement<SDateRangeVal>(buf, time_data, val_data);		break;
	case T_TimeStamp:			SetSignalElement<STimeStamp>(buf, time_data, val_data);			break;
	case T_TimeRangeVal:		SetSignalElement<STimeRangeVal>(buf, time_data, val_data);		break;
	case T_DateVal2:			SetSignalElement<SDateVal2>(buf, time_data, val_data);			break;
	case T_TimeLongVal:			SetSignalElement<STimeLongVal>(buf, time_data, val_data);		break;
	case T_DateShort2:			SetSignalElement<SDateShort2>(buf, time_data, val_data);		break;
	case T_ValShort2:			SetSignalElement<SValShort2>(buf, time_data, val_data);			break;
	case T_ValShort4:			SetSignalElement<SValShort4>(buf, time_data, val_data);			break;
		//case T_CompactDateVal:		SetSignalElement<SCompactDateVal>(buf, time_data, val_data);	break; // not fully supported yet
	default: MERR("ERROR:MedSignalsSingleElemFill Unknown sig_type %d\n", sig_type);
		return -1;

	}
	return 0;
}

int MedSignalsPrintVecByType(ostream &os, int sig_type, void* vec, int len_bytes)
{
	switch ((SigType)sig_type) {

	case T_Value:				MedSignalsPrintVec<SVal>(os, (SVal *)vec, len_bytes / sizeof(SVal));									break;
	case T_DateVal:				MedSignalsPrintVec<SDateVal>(os, (SDateVal *)vec, len_bytes / sizeof(SDateVal));						break;
	case T_TimeVal:				MedSignalsPrintVec<STimeVal>(os, (STimeVal *)vec, len_bytes / sizeof(STimeVal));						break;
	case T_DateRangeVal:		MedSignalsPrintVec<SDateRangeVal>(os, (SDateRangeVal *)vec, len_bytes / sizeof(SDateRangeVal));		break;
	case T_TimeStamp:			MedSignalsPrintVec<STimeStamp>(os, (STimeStamp *)vec, len_bytes / sizeof(STimeStamp));				break;
	case T_TimeRangeVal:		MedSignalsPrintVec<STimeRangeVal>(os, (STimeRangeVal *)vec, len_bytes / sizeof(STimeRangeVal));		break;
	case T_DateVal2:			MedSignalsPrintVec<SDateVal2>(os, (SDateVal2 *)vec, len_bytes / sizeof(SDateVal2));					break;
	case T_TimeLongVal:			MedSignalsPrintVec<STimeLongVal>(os, (STimeLongVal *)vec, len_bytes / sizeof(STimeLongVal));			break;
	case T_DateShort2:			MedSignalsPrintVec<SDateShort2>(os, (SDateShort2 *)vec, len_bytes / sizeof(SDateShort2));				break;
	case T_ValShort2:			MedSignalsPrintVec<SValShort2>(os, (SValShort2 *)vec, len_bytes / sizeof(SValShort2));				break;
	case T_ValShort4:			MedSignalsPrintVec<SValShort4>(os, (SValShort4 *)vec, len_bytes / sizeof(SValShort4));				break;
	default: MERR("ERROR: MedSignalsPrintVecByType Unknown sig_type %d\n", sig_type);
		return -1;

	}
	return 0;
}
