//
// MedSignals.h : Signal types definitions
//
#define __INFRAMED_DLL
#ifndef __MEDSIGNALS__H__
#define __MEDSIGNALS__H__

#include <string>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;

#define N_SignalTypes

enum SigType {T_Value = 0,		// 0 :: single float Value
			  T_DateVal,		// 1 :: date (32 bit space yyyymmdd reccomended) , float value (MOST COMMON !)
			  T_TimeVal,		// 2 :: date-time tag (64 bit space), float value
			  T_DateRangeVal,	// 3 :: date start, date end, float value
			  T_TimeStamp,		// 4 :: 64 bits of data (mainly for time stamps)
			  T_TimeRangeVal,	// 5 :: time-time + value
			  T_DateVal2,		// 6 :: date, float value, unsigned short additional value (specially tailored to drug code + drug period - to save a lot of space)	
			  T_TimeLongVal,	// 7 :: date-time (64 bit) + long long value
			  T_DateShort2,		// 8 :: date (32 bits) + 2 short values (perfect for BP for example).
			  T_ValShort2,		// 9 :: 2 short values
			  T_ValShort4,		// 10 :: 4 short values
			  T_CompactDateVal,	// 11 :: 2 unsigned shorts - first is a compact date (in 16 bits), second in an unsigned short value
			  T_Last};			// 12 :: next free slot for type id

int get_type_size(SigType t);

//====================================================================
// UniversalSigVec :
// -----------------
// A unified wrapper for signals to allow getting times and values
// from signals in a unified API.
//====================================================================
class UniversalSigVec {
	public:
		char *data;
		int len;		// type len (not bytes len)

		inline float val(int idx);
		int val_int(int idx);
		//inline float val_float(int idx) { return val(idx); }
		double val_double(int idx);
		long long val_long(int idx);

		int date(int idx);
		int days(int idx);
		int minutes(int idx);

		int time_days(int idx);
		int time_minutes(int idx);

		float val(int idx, int channel);

		int date(int idx, int channel);
		int days(int idx, int channel);
		int minutes(int idx, int channel);

		int init(SigType _type);

	private:
		// these will init when init() is called.
		SigType type; // type of the embedded signal
		int n_time_channels;
		int n_val_channels;
		int size_of_element;
		vector<int> time_channels_offsets;
		vector<int> value_channels_offsets;

		int val_cast_to_float_flag;
		int val_cast_to_double_flag;
		int val_cast_to__flag;

		// SDateVal funcs
		//inline float SDateVal_val_float(int idx, int chan) { return ((SDateVal *)data)[idx].val; }

};

//===========================================
// UnifiedSig - unifiying API's for signals
// This has only virtual functions
//===========================================
class UnifiedSig {
public:

};

//===================================
// SVal
//===================================
class SVal {
	public:
		float val;
};

//===================================
// SDateVal
//===================================
class SDateVal {
	public:
		int date;
		float val;

		inline int n_time_channels() { return 1; }
		inline int n_val_channels() { return 1; }
		inline float val_float(int chan) { return val; }
//		inline int date(int chan) { return date; }
//		inline int days(int chan) { return date_to_days_IM(date); }
};

//===================================
// STimeVal
//===================================
class STimeVal {
	public:
		long long time;
		float val;
};

//===================================
// SDateRangeVal
//===================================
class SDateRangeVal {
	public:
		int date_start;
		int date_end;
		float val;
};

//===================================
// STimeRangeVal
//===================================
class STimeRangeVal {
	public:
		long long time_start;
		long long time_end;
		float val;
};

//===================================
// STimeStamp
//===================================
class STimeStamp {
	public:
		long long time;
};

//===================================
// SDateVal2
//===================================
class SDateVal2 {
	public:
		int date;
		float val;
		unsigned short val2;
};

//===================================
// STimeLongVal
//===================================
class STimeLongVal {
	public:
		long long time ;
		long long val ;
};

//===================================
// SDateShort2
//===================================
class SDateShort2 {
public:
	int date;
	short val1;
	short val2;
};

//===================================
// SValShort2
//===================================
class SValShort2 {
public:
	short val1;
	short val2;
};


//===================================
// SValShort4
//===================================
class SValShort4 {
public:
	short val1;
	short val2;
	short val3;
	short val4;
};

//===================================
// SCompactDateVal
//===================================
class SCompactDateVal {
public:
	unsigned short compact_date;		// kept as: top 7 bits cY, then 4 bits M, then 5 bits day. Year is cY+1923 (format lasts until 2050)
	unsigned short val;
};

//==========================================
// Compact date to normal date conversions.
//==========================================
inline unsigned short date_to_compact_date(int date) {
	int d = date % 100; int m = (date-d)/100 % 100; int y=max(date/10000, 1923); y = y-1923; unsigned short cd = (y<<9)+(m<<5)+d; return cd;
}

inline int compact_date_to_date(unsigned short cd) {
	int d = cd & 0x1f; int m = (cd >> 5) & 0xf; int y = (m>>9); int date=y*10000+m*100+d; return date;
}

//=============================================================================================================
class SignalInfo {
	public:
		int sid;
		string name;
		int type;
		int bytes_len;
		string description;
		int fno; // currently each signal is in a single data and index file. This helps make things faster and is doable.
		int shift;
		float factor;

		SignalInfo() { fno = -1; };
};

//===================================================================
// Signals file handler
//===================================================================
class MedSignals {
	public:
		vector<string> fnames;
		string signals_to_files;
		map<string, int> Name2Sid;
		map<int, string> Sid2Name;
		map<int, SignalInfo> Sid2Info;
		vector<string> signals_names;
		vector<int>	signals_ids;
		vector<int> sid2serial; // inverse of signal_ids, -1: empty slots

		void clear() { fnames.clear(); Name2Sid.clear(); Sid2Name.clear(); signals_names.clear(); signals_ids.clear(); }

		int read(const string &fname);
		int read(vector<string> &sfnames);
		int read(string path, vector<string> &sfnames);

		int read_sfile(const string &fname);
		int read_sfile(const string &path, const string &sfile_name);

		inline int sid(const string &name);
		string name(int sid);
		int type(const string &name);
		int type(int sid);
		string desc(const string &name);
		string desc(int sid);
		int fno(const string &sig_name);
		int fno(int sid);

	private:

};

//=============================================================================================
// Inline functions
//=============================================================================================
inline int MedSignals::sid(const string &name) 
{ 
	if (Name2Sid.find(name) == Name2Sid.end()) 
		return -1; return 
	Name2Sid[name];
};

#endif