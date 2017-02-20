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
#include <MedTime/MedTime/MedTime.h>

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

namespace MedRep {
	int get_type_size(SigType t);
	int get_type_channels(SigType t, int &time_unit, int &n_time_chans, int &n_val_chans);
	template <class T> int get_type_channels_info(int &time_unit, int &n_time_chans, int &n_val_chans) {
		time_unit = T::time_unit();
		n_time_chans = T::n_time_channels();
		n_val_chans = T::n_val_channels();		
		return 0;
	}
}

//======================================================================
// UnifiedSig - unifiying API's for signals
// This has only virtual functions and functions built with them
// Never add a data member to UnifiedSig !
//======================================================================
class UnifiedSig {
public:

	// channels numbers
	inline int n_time_channels() { return 0; }
	inline int n_val_channels() { return 0; }
	
	// time unit & unitless time
	inline int time_unit() { return 0; }
	inline int Time(int chan) { return 0; }

	// value channels float
	inline float Val(int chan) { return 0; }
	inline float SetVal(int chan) { return 0; }

	// Following functions are implemented based on the functions above (and save lots of coding hence)
	// time channels int
	inline int Date(int chan) { return med_time_converter.convert_times(time_unit(), MedTime::Date, Time(chan)); }
	inline int Years(int chan) { return med_time_converter.convert_times(time_unit(), MedTime::Years, Time(chan)); }
	inline int Months(int chan) { return med_time_converter.convert_times(time_unit(), MedTime::Months, Time(chan)); }
	inline int Days(int chan) { return med_time_converter.convert_times(time_unit(), MedTime::Days, Time(chan)); }
	inline int Hours(int chan) { return med_time_converter.convert_times(time_unit(), MedTime::Hours, Time(chan)); }
	inline int Minutes(int chan) { return med_time_converter.convert_times(time_unit(), MedTime::Minutes, Time(chan)); }


	// channel 0 easy access
	inline int Date() { return Date(0); }
	inline int Years() { return Years(0); }
	inline int Months() { return Months(0); }
	inline int Days() { return Days(0); }
	inline int Hours() { return Hours(0); }
	inline int Minutes() { return Minutes(0); }
	inline float Val() { return Val(0); }
};

//===================================
// SVal
//===================================
class SVal : UnifiedSig {
	public:
		float val;

		// unified API extension
		static inline int n_time_channels() { return 0; }
		static inline int n_val_channels() { return 1; }
		static inline int time_unit() { return 0; }
		inline int Time(int chan) { return 0; }
		inline float Val(int chan) { return val; }
		inline void SetVal(int chan, float val) { return; };
		
};

//===================================
// SDateVal
//===================================
class SDateVal : UnifiedSig {
	public:
		int date;
		float val;

		// unified API extention
		static inline int n_time_channels() { return 1; }
		static inline int n_val_channels() { return 1; }
		static inline int time_unit() { return MedTime::Date; }
		inline int Time(int chan) { return date; }
		inline float Val(int chan) { return val; }
		inline void SetVal(int chan, float _val) { val = _val; };
};

//===================================
// STimeVal
//===================================
class STimeVal : UnifiedSig {
	public:
		long long time;
		float val;

		// unified API extention
		static inline int n_time_channels() { return 1; }
		static inline int n_val_channels() { return 1; }
		static inline int time_unit() { return MedTime::Minutes; }
		inline int Time(int chan) { return (int)time; } // assuming minutes span are within the size of an int
		inline float Val(int chan) { return val; }
		inline void SetVal(int chan, float _val) { val = _val; };
};

//===================================
// SDateRangeVal
//===================================
class SDateRangeVal : UnifiedSig {
	public:
		int date_start;
		int date_end;
		float val;

		// unified API extention
		static inline int n_time_channels() { return 2; }
		static inline int n_val_channels() { return 1; }
		static inline int time_unit() { return MedTime::Date; }
		inline int Time(int chan) { return ((chan) ? (date_end) : (date_start)); } // assuming minutes span are within the size of an int
		inline float Val(int chan) { return val; }
		inline void SetVal(int chan, float _val) { val = _val; };
};

//===================================
// STimeRangeVal
//===================================
class STimeRangeVal : UnifiedSig {
	public:
		long long time_start;
		long long time_end;
		float val;

		// unified API extention
		static inline int n_time_channels() { return 2; }
		static inline int n_val_channels() { return 1; }
		static inline int time_unit() { return MedTime::Minutes; }
		inline int Time(int chan) { return ((chan) ? ((int)time_end) : ((int)time_start)); } // assuming minutes span are within the size of an int
		inline float Val(int chan) { return val; }
		inline void SetVal(int chan, float _val) { val = _val; };
};

//===================================
// STimeStamp
//===================================
class STimeStamp : UnifiedSig {
	public:
		long long time;

		// unified API extention
		static inline int n_time_channels() { return 1; }
		static inline int n_val_channels() { return 0; }
		static inline int time_unit() { return MedTime::Minutes; }
		inline int Time(int chan) { return (int)time; } // assuming minutes span are within the size of an int
		inline float Val(int chan) { return 0; }
		inline void SetVal(int chan, float _val) { return; };

};

//===================================
// SDateVal2
//===================================
class SDateVal2 : UnifiedSig {
	public:
		int date;
		float val;
		unsigned short val2;

		// unified API extention
		static inline int n_time_channels() { return 1; }
		static inline int n_val_channels() { return 2; }
		static inline int time_unit() { return MedTime::Date; }
		inline int Time(int chan) { return date; } // assuming minutes span are within the size of an int
		inline float Val(int chan) { return ((chan) ? (float)val2 : (float)val) ; }
		inline void SetVal(int chan, float _val) { (chan) ? val2 = (unsigned short)_val : val = _val; };

};

//===================================
// STimeLongVal
//===================================
class STimeLongVal : UnifiedSig {
	public:
		long long time ;
		long long val ;

		// unified API extention
		static inline int n_time_channels() { return 1; }
		static inline int n_val_channels() { return 1; }
		static inline int time_unit() { return MedTime::Minutes; }
		inline int Time(int chan) { return (int)time; } // assuming minutes span are within the size of an int
		inline float Val(int chan) { return (float)val; }
		inline void SetVal(int chan, float _val) { val = (long long)_val; };

		// Waiting with unified here until we support long version of values.
};

//===================================
// SDateShort2
//===================================
class SDateShort2 : UnifiedSig {
public:
	int date;
	short val1;
	short val2;

	// unified API extention
	static inline int n_time_channels() { return 1; }
	static inline int n_val_channels() { return 2; }
	static inline int time_unit() { return MedTime::Date; }
	inline int Time(int chan) { return date; } // assuming minutes span are within the size of an int
	inline float Val(int chan) { return ((chan) ? (float)val2 : (float)val1); }
	inline void SetVal(int chan, float _val) { (chan) ? val2 = (short)_val : val1 = (short)_val; };

};

//===================================
// SValShort2
//===================================
class SValShort2 : UnifiedSig {
public:
	short val1;
	short val2;
	// unified API extention
	static inline int n_time_channels() { return 0; }
	static inline int n_val_channels() { return 2; }
	static inline int time_unit() { return 0; }
	inline int Time(int chan) { return 0; }
	inline float Val(int chan) { return ((chan) ? (float)val2 : (float)val1); }
	inline void SetVal(int chan, float _val) { (chan) ? val2 = (short)_val : val1 = (short)_val; };
};


//===================================
// SValShort4
//===================================
class SValShort4 : UnifiedSig {
public:
	short val1;
	short val2;
	short val3;
	short val4;

	// unified API extention
	static inline int n_time_channels() { return 0; }
	static inline int n_val_channels() { return 4; }
	static inline int time_unit() { return 0; }
	inline int Time(int chan) { return 0; }

	inline float Val(int chan) { 
		switch (chan) {
			case 0: return val1;
			case 1: return val2;
			case 2: return val3;
			case 3: return val4;
		}
		return 0;
	}
	inline void SetVal(int chan, float _val) { 
		switch (chan) {
			case 0: val1 = (short)_val; return;
			case 1: val2 = (short)_val; return;
			case 2: val3 = (short)_val; return;
			case 3: val4 = (short)_val; return;
		}
	};

};

//===================================
// SCompactDateVal
//===================================
class SCompactDateVal {
public:
	unsigned short compact_date;		// kept as: top 7 bits cY, then 4 bits M, then 5 bits day. Year is cY+1923 (format lasts until 2050)
	unsigned short val;

	// No unified support until we support compact_date as a date in MedTime
	// unified API extention
	static inline int n_time_channels() { return 1; }
	static inline int n_val_channels() { return 1; }
	static inline int time_unit() { return MedTime::Date; }
	inline int Time(int chan) { return (int)compact_date; } // assuming minutes span are within the size of an int
	inline float Val(int chan) { return (float)val; }
	inline void SetVal(int chan, float _val) { val = (unsigned short)_val; };
};


//====================================================================
// UniversalSigVec :
// -----------------
// A unified wrapper for signals to allow getting times and values
// from signals in a unified API.
//====================================================================
template <class T> class UnifiedSignalsAPIs {
public:

	static inline int Time_ch_vec(int idx, int chan, void *data) { return ((T *)data)[idx].Time(chan); }
	static inline float Val_ch_vec(int idx, int chan, void *data) { return ((T *)data)[idx].Val(chan); }
	static inline void SetVal_ch_vec(int idx, int chan, float _val, void *data) { ((T *)data)[idx].SetVal(chan, _val); }
	static inline size_t size() { return sizeof(T); }
};


class UniversalSigVec {
public:
	void *data;
	int len;		// type len (not bytes len)

	//--------------------------------------------------------------------------------------
	// function pointers - to be set before using the relevant type (use the init function)
	//--------------------------------------------------------------------------------------
	// channels numbers
	int (*n_time_channels)();
	int (*n_val_channels)();

	// time unit & unitless time
	int (*time_unit)();
	int(*Time_ch_vec)(int, int, void *); // Time(idx,chan)

	// value channels float
	float (*Val_ch_vec)(int, int, void *);
	void (*SetVal_ch_vec)(int, int, float, void *);

	size_t (*size)();

	// init function : call before using a certain type
	void init(SigType _type);
	void init(int _type) { return init((SigType)_type); }

	//--------------------------------------------------------------------------------------
	// Following are based on the pointed functions above
	//--------------------------------------------------------------------------------------
	// Following functions are implemented based on the functions above (and save lots of coding hence)
	// time channels int
	inline int Time(int idx, int chan) { return Time_ch_vec(idx, chan, data); }
	inline float Val(int idx, int chan) { return Val_ch_vec(idx, chan, data); }

	// channel 0 easy API
	inline int Time(int idx) { return Time(idx, 0); }
	inline float Val(int idx) { return Val(idx, 0); }

	inline int TimeU(int idx, int to_time_unit) { return med_time_converter.convert_times(time_unit(), to_time_unit, Time(idx)); }
	inline int Date(int idx) { return med_time_converter.convert_times(time_unit(), MedTime::Date, Time(idx)); }
	inline int Years(int idx) { return med_time_converter.convert_times(time_unit(), MedTime::Years, Time(idx)); }
	inline int Months(int idx) { return med_time_converter.convert_times(time_unit(), MedTime::Months, Time(idx)); }
	inline int Days(int idx) { return med_time_converter.convert_times(time_unit(), MedTime::Days, Time(idx)); }
	inline int Hours(int idx) { return med_time_converter.convert_times(time_unit(), MedTime::Hours, Time(idx)); }
	inline int Minutes(int idx) { return med_time_converter.convert_times(time_unit(), MedTime::Minutes, Time(idx)); }

	// general channel API
	inline int TimeU(int idx, int chan, int to_time_unit) { return med_time_converter.convert_times(time_unit(), to_time_unit, Time(idx, chan)); }
	inline int Date(int idx, int chan) { return med_time_converter.convert_times(time_unit(), MedTime::Date, Time(idx, chan)); }
	inline int Years(int idx, int chan) { return med_time_converter.convert_times(time_unit(), MedTime::Years, Time(idx, chan)); }
	inline int Months(int idx, int chan) { return med_time_converter.convert_times(time_unit(), MedTime::Months, Time(idx, chan)); }
	inline int Days(int idx, int chan) { return med_time_converter.convert_times(time_unit(), MedTime::Days, Time(idx, chan)); }
	inline int Hours(int idx, int chan) { return med_time_converter.convert_times(time_unit(), MedTime::Hours, Time(idx, chan)); }
	inline int Minutes(int idx, int chan) { return med_time_converter.convert_times(time_unit(), MedTime::Minutes, Time(idx, chan)); }

	template <class S> void set_funcs() {
		n_time_channels = &S::n_time_channels;
		n_val_channels = &S::n_val_channels;
		time_unit = &S::time_unit;
		Time_ch_vec = &UnifiedSignalsAPIs<S>::Time_ch_vec;
		Val_ch_vec = &UnifiedSignalsAPIs<S>::Val_ch_vec;
		SetVal_ch_vec = &UnifiedSignalsAPIs<S>::SetVal_ch_vec;
		size = &UnifiedSignalsAPIs<S>::size;
	}

	SigType get_type() { return type; }

private:
	SigType type = T_Last; // type of the embedded signal

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
		int time_unit;
		int n_time_channels;
		int n_val_channels;

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

//-----------------------------------------------------------------------------------------------
//template <class T> int MedRep:get_type_channels_info(int &time_unit, int &n_time_chans, int &n_val_chans)
//{
//	T t;
//
//	time_unit = t.time_unit();
//	n_time_chans = t.n_time_channels();
//	n_val_chans = t.n_val_channels();
//
//	return 0;
//}

#endif