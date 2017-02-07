//
// MedSignals.c
//

#ifndef _SCL_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS
#endif

#include "MedSignals.h"
#include "Logger/Logger/Logger.h"
#include <fstream>
#include <algorithm>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>


#define LOCAL_SECTION LOG_SIG
#define LOCAL_LEVEL LOG_DEF_LEVEL
extern MedLogger global_logger;

using namespace std;
using namespace boost;

//-----------------------------------------------------------------------------------------------
int get_type_size(SigType t)
{
	switch (t) {

		case T_Value: 
			return ((int)sizeof(SVal));

		case T_DateVal:
			return ((int)sizeof(SDateVal)) ;

		case T_TimeVal:
			return ((int)sizeof(STimeVal));

		case T_DateRangeVal:
			return ((int)sizeof(SDateRangeVal));
		
		case T_TimeRangeVal:
			return ((int) sizeof(STimeRangeVal));

		case T_TimeStamp:
			return ((int) sizeof(STimeStamp)) ;

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
			MERR("Cannot get size of signal type %d\n",t) ;
			return 0;
	}

	return 0;
}

//-----------------------------------------------------------------------------------------------
int MedSignals::read(vector<string> &sfnames)
{
	int rc = 0;

	for (int i=0; i<sfnames.size(); i++) {
		rc += read(sfnames[i]);
	}

	return rc;
}

//-----------------------------------------------------------------------------------------------
int MedSignals::read(string path, vector<string> &sfnames)
{
	int rc = 0;

	for (int i=0; i<sfnames.size(); i++) {
		string fname = (path == "") ? sfnames[i] : path + "/" + sfnames[i] ;
		rc += read(fname);
	}

	return rc;
}

//-----------------------------------------------------------------------------------------------
int MedSignals::read(const string &fname)
{

	ifstream inf(fname);

	if (!inf) {
		MERR("MedSignals: read: Can't open file %s\n",fname.c_str());
		return -1;
	}

	fnames.push_back(fname); // TBD : check that we didn't already load this file
	string curr_line;
	MLOG_D("Working on signals file %s\n",fname.c_str());
	while (getline(inf,curr_line)) {
		if ((curr_line.size() > 1) && (curr_line[0] != '#')) {
			vector<string> fields;
			split(fields, curr_line, boost::is_any_of("\t"));
			MLOG_D("MedSignals: read: file %s line %s\n",fname.c_str(),curr_line.c_str());
			if (fields.size() >= 4 && fields.size() <= 5 && fields[0].compare(0,6,"SIGNAL")==0) {
				// line format: SIGNAL <name> <signal id> <signal type num> <description>

				int sid = stoi(fields[2]);
				if (Name2Sid.find(fields[1]) != Name2Sid.end() || Sid2Name.find(sid) != Sid2Name.end()) {
					MERR("MedSignals: read: name or id already used in line %s\n",curr_line.c_str());
					return -1;
				} else {
					Name2Sid[fields[1]] = sid;
					Sid2Name[sid] = fields[1];
					signals_names.push_back(fields[1]);
					signals_ids.push_back(sid);
					int type = stoi(fields[3]);
					if (type<0 || type>(int)T_Last) {
						MERR("MedSignals: read: type %d not recognized in line %s\n",curr_line.c_str());
						return -1;
					}
					SignalInfo info;
					info.sid = sid;
					info.name = fields[1];
					info.type = type;
					info.bytes_len = get_type_size((SigType)type);
					if (fields.size() == 4)
						info.description = "";
					else
						info.description = fields[4];
					Sid2Info[sid] = info;
				}

			} else {
				MWARN("MedSignals: read: can't parse line: %s (%d)\n",curr_line.c_str(),fields.size());
			}

		}
	}

	sort(signals_ids.begin(), signals_ids.end());
	int max_sid = *max_element(signals_ids.begin(), signals_ids.end());
	sid2serial.clear();
	sid2serial.resize(max_sid+1, -1);
	for (int i=0; i<signals_ids.size(); i++)
		sid2serial[signals_ids[i]] = i;


	inf.close();
	fnames.push_back(fname);
	MLOG_D("Finished reading signals file %s\n",fname.c_str());
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
					Sid2Info[sid].fno = stoi(fields[0]);
					//MLOG("read_sfile: %s\n", curr_line.c_str());
					//MLOG("read_sfile: sid %d fno %d\n", sid, Sid2Info[sid].fno);
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
	if (Sid2Info.find(sid) == Sid2Info.end())
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
	if (Sid2Info.find(sid) == Sid2Info.end())
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
	if (Sid2Info.find(sid) == Sid2Info.end())
		return -1;
	return Sid2Info[sid].fno;
}

//================================================================================================
// UniversalSigVec
//================================================================================================

//-----------------------------------------------------------------------------------------------
//float UniversalSigVec::val(int idx)
//{
//	// idx is assumed to be < len ! (saving the time for checking it)
//	return ((float *)&data[idx*size_of_element + value_channels_offsets[0]])[0];
//}

//int val_int(int idx);
//long long val_long(int idx);
//int time_days(int idx);
//int time_minutes(int idx);

//float val(int idx, int channel);

//int time_days(int idx, int channel);
//int time_minutes(int idx, int channel);

//-----------------------------------------------------------------------------------------------
int UniversalSigVec::init(SigType _type)
{
	switch (_type) {
	case T_DateVal:
		type = T_DateVal;
		n_time_channels = 1;
		n_val_channels = 1;
		size_of_element = (int)sizeof(SDateVal);
		time_channels_offsets = { offsetof(SDateVal, date) };
		value_channels_offsets = { offsetof(SDateVal, val) };
		break;

	case T_TimeVal:
		type = T_TimeVal;
		n_time_channels = 1;
		n_val_channels = 1;
		size_of_element = (int)sizeof(STimeVal);
		time_channels_offsets ={ offsetof(STimeVal, time) };
		value_channels_offsets ={ offsetof(SDateVal, val) };
		break;
	}

	return 0;
}
