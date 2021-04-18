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
#include <iostream>
#include <memory>
#include <cstring>
#include <array>

using namespace std;

#define N_SignalTypes
#define GENERIC_SIG_VEC_MAX_CHANNELS 30


#ifndef USE_LEGACY_USV
class GenericSigVec;
class GenericSigVec_mem;
typedef class GenericSigVec UniversalSigVec;
typedef class GenericSigVec_mem UniversalSigVec_mem;
#else // legacy USV
class UniversalSigVec_legacy;
class UniversalSigVec_mem_legacy;
typedef class UniversalSigVec_legacy UniversalSigVec;
typedef class UniversalSigVec_mem_legacy UniversalSigVec_mem;
#endif


enum SigType {
	T_Value = 0,		// 0 :: single float Value
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
	T_DateRangeVal2,	// 12 :: date start, date end, 2 float values
	T_DateFloat2,	// 13 :: date + 2 float values
	T_TimeRange,	// 14 :: time-time
	T_TimeShort4,   // 15 :: time + 4 shorts
	T_Generic,	// 16 :: generic signal
	T_Last
};		//    :: next free slot for type id

namespace MedRep {
	int get_type_size(SigType t);
	int get_type_channels(SigType t, int &time_unit, int &n_time_chans, int &n_val_chans);
	int get_type_channels(const string& sigSpec, int &time_unit, int &n_time_chans, int &n_val_chans);
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
	inline void SetVal(int chan, float _val) {}

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

//=============================================================================================
// General Fill in function
//=============================================================================================
int MedSignalsSingleElemFill(int sig_type, char *buf, int *time_data, float *val_data);


//===================================
// SVal
//===================================
class SVal : public UnifiedSig {
public:
	float val;

	// unified API extension
	static inline int n_time_channels() { return 0; }
	static inline int n_val_channels() { return 1; }
	static inline int time_unit() { return 0; }
	inline int Time(int chan) { return 0; }
	inline float Val(int chan) { return val; }
	inline void SetVal(int chan, float _val) { val = _val; };

	inline void Set(float _val) { val = _val; }
	inline void Set(int *times, float *vals) { val = vals[0]; }

	bool operator<(const SVal& s) { return (this->val < s.val); }
	bool operator==(const SVal& s) { return (this->val == s.val); }

	friend ostream& operator<<(ostream& os, const SVal& s) { os << s.val; return os; }
};


//===================================
// SDateVal
//===================================
class SDateVal : public UnifiedSig {
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

	inline void Set(int _date, float _val) { date = _date; val = _val; }
	inline void Set(int *times, float *vals) { date = times[0]; val = vals[0]; }

	bool operator<(const SDateVal& s) { if (this->date < s.date) return true; if (this->date > s.date) return false; return (this->val < s.val); }
	bool operator==(const SDateVal& s) { return (this->val == s.val && this->date == s.date); }

	friend ostream& operator<<(ostream& os, const SDateVal& s) { os << s.date << ":" << s.val; return os; }
};

//===================================
// STimeVal
//===================================
class STimeVal : public UnifiedSig {
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

	inline void Set(long long _time, float _val) { time = _time; val = _val; }
	inline void Set(int *times, float *vals) { time = (long long)times[0]; val = vals[0]; }

	bool operator<(const STimeVal& s) { if (this->time < s.time) return true; if (this->time > s.time) return false; return (this->val < s.val); }
	bool operator==(const STimeVal& s) { return (this->val == s.val && this->time == s.time); }

	friend ostream& operator<<(ostream& os, const STimeVal& s) { os << s.time << ":" << s.val; return os; }

};

//===================================
// SDateRangeVal
//===================================
class SDateRangeVal : public UnifiedSig {
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

	inline void Set(int _date_start, int _date_end, float _val) { date_start = _date_start; date_end = _date_end; val = _val; }
	inline void Set(int *times, float *vals) { date_start = times[0]; date_end = times[1]; val = vals[0]; }

	bool operator<(const SDateRangeVal& s) {
		if (this->date_start < s.date_start) return true;
		if (this->date_start > s.date_start) return false;
		if (this->date_end < s.date_end) return true;
		if (this->date_end > s.date_end) return false;
		return (this->val < s.val);
	}
	bool operator==(const SDateRangeVal& s) { return (this->val == s.val && this->date_start == s.date_start && this->date_end == s.date_end); }

	friend ostream& operator<<(ostream& os, const SDateRangeVal& s) { os << s.date_start << "-" << s.date_end << ":" << s.val; return os; }

};

//===================================
// SDateRangeVal2
//===================================
class SDateRangeVal2 : public UnifiedSig {
public:
	int date_start;
	int date_end;
	float val;
	float val2;

	// unified API extention
	static inline int n_time_channels() { return 2; }
	static inline int n_val_channels() { return 2; }
	static inline int time_unit() { return MedTime::Date; }
	inline int Time(int chan) { return ((chan) ? (date_end) : (date_start)); } // assuming minutes span are within the size of an int
	inline float Val(int chan) { return ((chan) ? (float)val2 : (float)val); }
	inline void SetVal(int chan, float _val) { (chan) ? val2 = _val : val = _val; };

	inline void Set(int _date_start, int _date_end, float _val, float _val2) { date_start = _date_start; date_end = _date_end; val = _val; val2 = _val2; }
	inline void Set(int *times, float *vals) { date_start = times[0]; date_end = times[1]; val = vals[0]; val2 = vals[1]; }

	bool operator<(const SDateRangeVal2& s) {
		if (this->date_start < s.date_start) return true;
		if (this->date_start > s.date_start) return false;
		if (this->date_end < s.date_end) return true;
		if (this->date_end > s.date_end) return false;
		if (this->val < s.val) return true;
		if (this->val > s.val) return false;
		return (this->val2 < s.val2);
	}
	bool operator==(const SDateRangeVal2& s) { return (this->val == s.val && this->val2 == s.val2 && this->date_start == s.date_start && this->date_end == s.date_end); }

	friend ostream& operator<<(ostream& os, const SDateRangeVal2& s) { os << s.date_start << "-" << s.date_end << ":" << s.val << "," << s.val2; return os; }

};

//===================================
// SDateFloat2
//===================================
class SDateFloat2 : public UnifiedSig {
public:
	int date;
	float val;
	float val2;

	// unified API extention
	static inline int n_time_channels() { return 1; }
	static inline int n_val_channels() { return 2; }
	static inline int time_unit() { return MedTime::Date; }
	inline int Time(int chan) { return date; } // assuming minutes span are within the size of an int
	inline float Val(int chan) { return ((chan) ? (float)val2 : (float)val); }
	inline void SetVal(int chan, float _val) { (chan) ? val2 = _val : val = _val; };

	inline void Set(int _date, float _val, float _val2) { date = _date; val = _val; val2 = _val2; }
	inline void Set(int *times, float *vals) { date = times[0]; val = vals[0]; val2 = vals[1]; }

	bool operator<(const SDateFloat2& s) {
		if (this->date < s.date) return true;
		if (this->date > s.date) return false;
		if (this->val < s.val) return true;
		if (this->val > s.val) return false;
		return (this->val2 < s.val2);
	}
	bool operator==(const SDateFloat2& s) { return (this->val == s.val && this->val2 == s.val2 && this->date == s.date); }

	friend ostream& operator<<(ostream& os, const SDateFloat2& s) { os << s.date << ":" << s.val << "," << s.val2; return os; }

};

//===================================
// STimeRangeVal
//===================================
class STimeRangeVal : public UnifiedSig {
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

	inline void Set(long long _time_start, long long _time_end, float _val) { time_start = _time_start; time_end = _time_end; val = _val; }
	inline void Set(int *times, float *vals) { time_start = (long long)times[0]; time_end = (long long)times[1]; val = vals[0]; }

	bool operator<(const STimeRangeVal& s) {
		if (this->time_start < s.time_start) return true;
		if (this->time_start > s.time_start) return false;
		if (this->time_end < s.time_end) return true;
		if (this->time_end > s.time_end) return false;
		return (this->val < s.val);
	}
	bool operator==(const STimeRangeVal& s) { return (this->val == s.val && this->time_start == s.time_start && this->time_end == s.time_end); }

	friend ostream& operator<<(ostream& os, const STimeRangeVal& s) { os << s.time_start << "-" << s.time_end << ":" << s.val; return os; }

};

//===================================
// STimeRange - 14
//===================================
class STimeRange : public UnifiedSig {
public:
	long long time_start;
	long long time_end;

	// unified API extention
	static inline int n_time_channels() { return 2; }
	static inline int n_val_channels() { return 0; }
	static inline int time_unit() { return MedTime::Minutes; }
	inline int Time(int chan) { return ((chan) ? ((int)time_end) : ((int)time_start)); } // assuming minutes span are within the size of an int
	inline float Val(int chan) { return 0; }
	inline void SetVal(int chan, float _val) { return; };

	inline void Set(long long _time_start, long long _time_end) { time_start = _time_start; time_end = _time_end; }
	inline void Set(int *times, float *vals) { time_start = (long long)times[0]; time_end = (long long)times[1]; }

	bool operator<(const STimeRange& s) {
		if (this->time_start < s.time_start) return true;
		if (this->time_start > s.time_start) return false;
		if (this->time_end < s.time_end) return true;
		return false;
	}
	bool operator==(const STimeRange& s) { return (this->time_start == s.time_start && this->time_end == s.time_end); }

	friend ostream& operator<<(ostream& os, const STimeRange& s) { os << s.time_start << "-" << s.time_end; return os; }

};

//===================================
// SDateShort4
//===================================
class STimeShort4 : public UnifiedSig {
public:
	long long time;
	unsigned short val1;
	unsigned short val2;
	unsigned short val3;
	unsigned short val4;

	// unified API extention
	static inline int n_time_channels() { return 1; }
	static inline int n_val_channels() { return 4; }
	static inline int time_unit() { return MedTime::Minutes; }
	inline int Time(int chan) { return (int)time; } // assuming minutes span are within the size of an int
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
		case 0: val1 = (unsigned short)_val; return;
		case 1: val2 = (unsigned short)_val; return;
		case 2: val3 = (unsigned short)_val; return;
		case 3: val4 = (unsigned short)_val; return;
		}
	};

	inline void Set(int _time, short _val1, short _val2, short _val3, short _val4) { time = _time; val1 = _val1; val2 = _val2; val3 = _val3; val4 = _val4; }
	inline void Set(int *times, float *vals) {
		time = times[0];
		val1 = (unsigned short)vals[0];
		val2 = (unsigned short)vals[1];
		val3 = (unsigned short)vals[2];
		val4 = (unsigned short)vals[3];
	}

	bool operator<(const STimeShort4& s) {
		if (this->time < s.time) return true;
		if (this->time > s.time) return false;
		if (this->val1 < s.val1) return true;
		if (this->val1 > s.val1) return false;
		if (this->val2 < s.val2) return true;
		if (this->val2 > s.val2) return false;
		if (this->val3 < s.val3) return true;
		if (this->val3 > s.val3) return false;
		return (this->val4 < s.val4);
	}
	bool operator==(const STimeShort4& s) {
		return (this->time == s.time
			&& this->val1 == s.val1
			&& this->val2 == s.val2
			&& this->val3 == s.val3
			&& this->val4 == s.val4);
	}

	friend ostream& operator<<(ostream& os, const STimeShort4& s) {
		os << s.time << ":"
			<< s.val1 << ","
			<< s.val2 << ","
			<< s.val3 << ","
			<< s.val4;
		return os;
	}
};



//===================================
// STimeStamp
//===================================
class STimeStamp : public UnifiedSig {
public:
	long long time;

	// unified API extention
	static inline int n_time_channels() { return 1; }
	static inline int n_val_channels() { return 0; }
	static inline int time_unit() { return MedTime::Minutes; }
	inline int Time(int chan) { return (int)time; } // assuming minutes span are within the size of an int
	inline float Val(int chan) { return 0; }
	inline void SetVal(int chan, float _val) { return; };

	inline void Set(long long _time) { time = _time; }
	inline void Set(int *times, float *vals) { time = (long long)times[0]; }

	bool operator<(const STimeStamp& s) { if (this->time < s.time) return true; return false; }
	bool operator==(const STimeStamp& s) { return (this->time == s.time); }

	friend ostream& operator<<(ostream& os, const STimeStamp& s) { os << s.time; return os; }

};

//===================================
// SDateVal2
//===================================
class SDateVal2 : public UnifiedSig {
public:
	int date;
	float val;
	unsigned short val2;

	// unified API extention
	static inline int n_time_channels() { return 1; }
	static inline int n_val_channels() { return 2; }
	static inline int time_unit() { return MedTime::Date; }
	inline int Time(int chan) { return date; } // assuming minutes span are within the size of an int
	inline float Val(int chan) { return ((chan) ? (float)val2 : (float)val); }
	inline void SetVal(int chan, float _val) { (chan) ? val2 = (unsigned short)_val : val = _val; };

	inline void Set(int _date, float _val, unsigned short _val2) { date = _date; val = _val; val2 = _val2; }
	inline void Set(int *times, float *vals) { date = times[0]; val = vals[0]; val2 = (unsigned short)vals[1]; }

	bool operator<(const SDateVal2& s) {
		if (this->date < s.date) return true;
		if (this->date > s.date) return false;
		if (this->val < s.val) return true;
		if (this->val > s.val) return false;
		return (this->val2 < s.val2);
	}
	bool operator==(const SDateVal2& s) { return (this->date == s.date && this->val == s.val && this->val2 == s.val2); }

	friend ostream& operator<<(ostream& os, const SDateVal2& s) { os << s.date << ":" << s.val << "," << s.val2; return os; }

};

//===================================
// STimeLongVal
//===================================
class STimeLongVal : public UnifiedSig {
public:
	long long time;
	long long val;

	// unified API extention
	static inline int n_time_channels() { return 1; }
	static inline int n_val_channels() { return 1; }
	static inline int time_unit() { return MedTime::Minutes; }
	inline int Time(int chan) { return (int)time; } // assuming minutes span are within the size of an int
	inline float Val(int chan) { return (float)val; }
	inline void SetVal(int chan, float _val) { val = (long long)_val; };

	inline void Set(long long _time, long long _val) { time = _time; val = _val; }
	// Waiting with unified here until we support long version of values.
	inline void Set(int *times, float *vals) { time = (long long)times[0]; val = (long long)vals[0]; }

	bool operator<(const STimeLongVal& s) {
		if (this->time < s.time) return true;
		if (this->time > s.time) return false;
		return (this->val < s.val);
	}
	bool operator==(const STimeLongVal& s) { return (this->time == s.time && this->val == s.val); }

	friend ostream& operator<<(ostream& os, const STimeLongVal& s) { os << s.time << ":" << s.val; return os; }

};

//===================================
// SDateShort2
//===================================
class SDateShort2 : public UnifiedSig {
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

	inline void Set(int _date, short _val1, short _val2) { date = _date; val1 = _val1; val2 = _val2; }
	inline void Set(int *times, float *vals) { date = times[0]; val1 = (short)vals[0]; val2 = (short)vals[1]; }

	bool operator<(const SDateShort2& s) {
		if (this->date < s.date) return true;
		if (this->date > s.date) return false;
		if (this->val1 < s.val1) return true;
		if (this->val1 > s.val1) return false;
		return (this->val2 < s.val2);
	}
	bool operator==(const SDateShort2& s) { return (this->date == s.date && this->val1 == s.val1 && this->val2 == s.val2); }

	friend ostream& operator<<(ostream& os, const SDateShort2& s) { os << s.date << ":" << s.val1 << "," << s.val2; return os; }
};

//===================================
// SValShort2
//===================================
class SValShort2 : public UnifiedSig {
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

	inline void Set(short _val1, short _val2) { val1 = _val1; val2 = _val2; }
	inline void Set(int *times, float *vals) { val1 = (short)vals[0]; val2 = (short)vals[1]; }

	bool operator<(const SValShort2& s) {
		if (this->val1 < s.val1) return true;
		if (this->val1 > s.val1) return false;
		return (this->val2 < s.val2);
	}
	bool operator==(const SValShort2& s) { return (this->val1 == s.val1 && this->val2 == s.val2); }

	friend ostream& operator<<(ostream& os, const SValShort2& s) { os << s.val1 << "," << s.val2; return os; }

};


//===================================
// SValShort4
//===================================
class SValShort4 : public UnifiedSig {
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

	inline void Set(short _val1, short _val2, short _val3, short _val4) { val1 = _val1; val2 = _val2; val3 = _val3; val4 = _val4; }
	inline void Set(int *times, float *vals) { val1 = (short)vals[0]; val2 = (short)vals[1]; val3 = (short)vals[2]; val4 = (short)vals[3]; }

	bool operator<(const SValShort4& s) {
		if (this->val1 < s.val1) return true;
		if (this->val1 > s.val1) return false;
		if (this->val2 < s.val2) return true;
		if (this->val1 > s.val1) return false;
		if (this->val3 < s.val3) return true;
		if (this->val1 > s.val1) return false;
		return (this->val4 < s.val4);
	}
	bool operator==(const SValShort4& s) { return (this->val1 == s.val1 && this->val2 == s.val2 && this->val3 == s.val3 && this->val4 == s.val4); }
	friend ostream& operator<<(ostream& os, const SValShort4& s) { os << s.val1 << "," << s.val2 << "," << s.val3 << "," << s.val4; return os; }

};


template <class T> void SetSignalElement(void *elem_buf, int *times, float *vals) { (*(T *)elem_buf).Set(times, vals); }
template <class T> int MedSignalsCompareSig(const void *a, const void *b) { if (*(T*)a < *(T*)b) return -1; if (*(T*)a == *(T*)b) return 0; return 1; }
template <class T> int MedSignalsPrintVec(ostream& os, T *vec, int n_elem);
int MedSignalsPrintVecByType(ostream &os, int sig_type, void* vec, int len);

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
	static inline void Set(int idx, int *_times, float *_vals, void *data) { ((T *)data)[idx].Set(_times, _vals); }
	static inline size_t size() { return sizeof(T); }
};

//==========================================
// Compact date to normal date conversions.
//==========================================
inline unsigned short date_to_compact_date(int date) {
	int d = date % 100; int m = (date - d) / 100 % 100; int y = max(date / 10000, 1923); y = y - 1923; unsigned short cd = (y << 9) + (m << 5) + d; return cd;
}

inline int compact_date_to_date(unsigned short cd) {
	int d = cd & 0x1f; int m = (cd >> 5) & 0xf; int y = (m >> 9); int date = y * 10000 + m * 100 + d; return date;
}

//=============================================================================================================
class SignalInfo {
public:
	std::array<int, GENERIC_SIG_VEC_MAX_CHANNELS> time_channel_offsets;
	std::array<int, GENERIC_SIG_VEC_MAX_CHANNELS> val_channel_offsets;
	std::array<unsigned char, GENERIC_SIG_VEC_MAX_CHANNELS> time_channel_types;
	std::array<unsigned char, GENERIC_SIG_VEC_MAX_CHANNELS> val_channel_types;

	int sid;
	string name;
	int type;
	int bytes_len;
	string description;
	string generic_signal_spec;
	int fno; // currently each signal is in a single data and index file. This helps make things faster and is doable.
	int shift;
	float factor;
	int time_unit;
	int n_time_channels;
	int n_val_channels;
	int virtual_sig = 0; // flag to tell if the signal was defined in the signals files OR if it was defined as a virtual signal
	std::array<int, GENERIC_SIG_VEC_MAX_CHANNELS> is_categorical_per_val_channel ; // when 1, channel doens't hold numerical values but rather pointers to a dict
	std::array<string, GENERIC_SIG_VEC_MAX_CHANNELS> unit_of_measurement_per_val_channel;

	SignalInfo() { fno = -1; time_channel_types.fill(0); val_channel_types.fill(0); time_channel_offsets.fill(0); val_channel_offsets.fill(0); is_categorical_per_val_channel.fill(0); unit_of_measurement_per_val_channel.fill(""); };

	void set_gsv_spec(const string &gsv_spec_str);

};
class MedRepository;
//===================================================================
// Signals file handler
//===================================================================
class MedSignals {
	int _allocate_new_signal(const string &sig_name);
public:
	vector<string> fnames;
	string signals_to_files;
	map<string, int> Name2Sid;
	map<int, string> Sid2Name;
	vector<SignalInfo> Sid2Info;
	vector<string> signals_names;
	vector<int>	signals_ids;
	vector<int> sid2serial; // inverse of signal_ids, -1: empty slots
	MedRepository* my_repo = NULL; // backward pointer to the owning repo

	void clear() { fnames.clear(); Name2Sid.clear(); Sid2Name.clear(); signals_names.clear(); signals_ids.clear(); }

	int read(const string &fname);
	int read(vector<string> &sfnames);
	int read(string path, vector<string> &sfnames);

	inline int sid(const string &name);
	string name(int sid);
	int type(const string &name);
	int type(int sid);
	string desc(const string &name);
	string desc(int sid);
	int fno(const string &sig_name);
	int fno(int sid);
	int has_any_categorical_channel(const string &sig_name);
	int has_any_categorical_channel(int sid);
	int is_categorical_channel(const string &sig_name, int val_channel);
	int is_categorical_channel(int sid, int val_channel);
	string unit_of_measurement(const string &name, int val_channel);
	string unit_of_measurement(int sid, int val_channel);
	// this option allows adding new signals definitions to the class, that were not defined in the files.
	// this is useful when using repositories to calculate new features, etc.
	int insert_virtual_signal(const string &sig_name, int type);
	int insert_virtual_signal(const string &sig_name, const string& signalSpec);
	int get_sids(vector<string> &sigs, vector<int> &sids);
};


class UniversalSigVec_legacy {
public:
	void *data;
	int len;		// type len (not bytes len)

					//--------------------------------------------------------------------------------------
					// function pointers - to be set before using the relevant type (use the init function)
					//--------------------------------------------------------------------------------------
					// channels numbers
	int(*n_time_channels)();
	int(*n_val_channels)();

	// time unit & unitless time
	int time_unit() const { return _time_unit; };
	int(*Time_ch_vec)(int, int, void *); // Time(idx,chan)

										 // value channels float
	float(*Val_ch_vec)(int, int, void *);
	void(*SetVal_ch_vec)(int, int, float, void *);

	// Set() - setting a specific index in data, given all its time channels and val channels
	// channels are given as an array (of the proper length) or NULL for 0 length channels.
	void(*Set)(int, int *, float *, void *);

	size_t(*size)();

	// init function : call before using a certain type
	void init(const SignalInfo &info);
	//void init(int _type) { return init((SigType)_type); }

	//--------------------------------------------------------------------------------------
	// Following are based on the pointed functions above
	//--------------------------------------------------------------------------------------
	// Following functions are implemented based on the functions above (and save lots of coding hence)
	// time channels int
	inline int Time(int idx, int chan) const { return Time_ch_vec(idx, chan, data); }
	inline float Val(int idx, int chan) const { return Val_ch_vec(idx, chan, data); }

	// channel 0 easy API
	inline int Time(int idx) const { return Time(idx, 0); }
	inline float Val(int idx) const { return Val(idx, 0); }

	inline int TimeU(int idx, int to_time_unit) const { return med_time_converter.convert_times(time_unit(), to_time_unit, Time(idx)); }
	inline int Date(int idx) const { return med_time_converter.convert_times(time_unit(), MedTime::Date, Time(idx)); }
	inline int Years(int idx) const { return med_time_converter.convert_times(time_unit(), MedTime::Years, Time(idx)); }
	inline int Months(int idx) const { return med_time_converter.convert_times(time_unit(), MedTime::Months, Time(idx)); }
	inline int Days(int idx) const { return med_time_converter.convert_times(time_unit(), MedTime::Days, Time(idx)); }
	inline int Hours(int idx) const { return med_time_converter.convert_times(time_unit(), MedTime::Hours, Time(idx)); }
	inline int Minutes(int idx) const { return med_time_converter.convert_times(time_unit(), MedTime::Minutes, Time(idx)); }

	// general channel API
	inline int TimeU(int idx, int chan, int to_time_unit) const { return med_time_converter.convert_times(time_unit(), to_time_unit, Time(idx, chan)); }
	inline int Date(int idx, int chan) const { return med_time_converter.convert_times(time_unit(), MedTime::Date, Time(idx, chan)); }
	inline int Years(int idx, int chan) const { return med_time_converter.convert_times(time_unit(), MedTime::Years, Time(idx, chan)); }
	inline int Months(int idx, int chan) const { return med_time_converter.convert_times(time_unit(), MedTime::Months, Time(idx, chan)); }
	inline int Days(int idx, int chan) const { return med_time_converter.convert_times(time_unit(), MedTime::Days, Time(idx, chan)); }
	inline int Hours(int idx, int chan) const { return med_time_converter.convert_times(time_unit(), MedTime::Hours, Time(idx, chan)); }
	inline int Minutes(int idx, int chan) const { return med_time_converter.convert_times(time_unit(), MedTime::Minutes, Time(idx, chan)); }

	template <class S> void set_funcs() {
		n_time_channels = &S::n_time_channels;
		n_val_channels = &S::n_val_channels;
		Time_ch_vec = &UnifiedSignalsAPIs<S>::Time_ch_vec;
		Val_ch_vec = &UnifiedSignalsAPIs<S>::Val_ch_vec;
		SetVal_ch_vec = &UnifiedSignalsAPIs<S>::SetVal_ch_vec;
		Set = &UnifiedSignalsAPIs<S>::Set;
		size = &UnifiedSignalsAPIs<S>::size;
	}

	SigType get_type() const { return type; }

	// helper functions for common operations
	//----------------------------------------
	// returns the first index i in the usv that has Time(i, time_chan) > time_bound, if none : return -1
	int get_index_gt_time_bound(int time_chan, int time_bound);
	// returns the first index i in the usv that has Time(i, time_chan) >= time_bound, if none : return -1
	int get_index_ge_time_bound(int time_chan, int time_bound);

protected:
	SigType type = T_Last; // type of the embedded signal
	int _time_unit = MedTime::Undefined;
};

/**
* Managed memory version of UniversalSigVec.
* for example when allocating new virtual records
*/
class UniversalSigVec_mem_legacy : public UniversalSigVec_legacy {
public:
	bool manage;
	UniversalSigVec_mem_legacy() {
		manage = false;
	}

	~UniversalSigVec_mem_legacy() {
		if (len > 0 && data != NULL && manage) {
			delete[](char *)data;
			data = NULL;
			len = 0;
		}
	}

	void set(const UniversalSigVec_legacy &s) {
		data = s.data;
		len = s.len;
		manage = false;
		n_time_channels = s.n_time_channels;
		n_val_channels = s.n_val_channels;
		Set = s.Set;
		SetVal_ch_vec = s.SetVal_ch_vec;
		size = s.size;
		Time_ch_vec = s.Time_ch_vec;
		_time_unit = s.time_unit();
		Val_ch_vec = s.Val_ch_vec;
		type = s.get_type();
	}
};


// small helper class to help in the sorting
struct GSVElement {
	vector<unsigned char> data;
};


//=============================================================================================
// Inline/templated functions
//=============================================================================================
inline int MedSignals::sid(const string &name)
{
	if (Name2Sid.find(name) == Name2Sid.end())
		return -1;
	return Name2Sid[name];
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


template <class T> int MedSignalsPrintVec(ostream& os, T *vec, int n_elem)
{
	for (int i = 0; i < n_elem; i++)
		os << vec[i] << " ";
	return 0;
}

class GenericSigVec {
public:
	class type_enc
	{
	public:
		static constexpr const unsigned char UNDEFINED = 0b00000000;
		static constexpr const unsigned char SIGNED = 0b00000001;
		static constexpr const unsigned char UINT8 = 0b00000010;  //unsigned char
		static constexpr const unsigned char UINT16 = 0b00000100;  //unsigned short
		static constexpr const unsigned char UINT32 = 0b00001000;  //unsigned int
		static constexpr const unsigned char UINT64 = 0b00010000;  //unsigned long long
		static constexpr const unsigned char INT8 = 0b00000011;  //char
		static constexpr const unsigned char INT16 = 0b00000101;  //short
		static constexpr const unsigned char INT32 = 0b00001001;  //int
		static constexpr const unsigned char INT64 = 0b00010001;  //long long
		static constexpr const unsigned char FLOAT32 = 0b00100000;  //float
		static constexpr const unsigned char FLOAT64 = 0b01000000;  //double
		static constexpr const unsigned char FLOAT80 = 0b10000000;  //long double
		static unsigned char encode(char c, bool isSigned = false) {
			unsigned char _is_signed = isSigned ? type_enc::SIGNED : 0;
			switch (c) {
			case 'c': return _is_signed | type_enc::UINT8;
			case 's': return _is_signed | type_enc::UINT16;
			case 'i': return _is_signed | type_enc::UINT32;
			case 'l': return _is_signed | type_enc::UINT64;
			case 'f': return type_enc::FLOAT32;
			case 'd': return type_enc::FLOAT64;
			case 'D': return type_enc::FLOAT80;
			}
			return type_enc::UNDEFINED;
		}
		static char decode(unsigned char c, bool isSigned = false) {
			switch (c) {
			case  type_enc::INT8:
			case  type_enc::UINT8: return 'c';
			case  type_enc::INT16:
			case  type_enc::UINT16: return 's';
			case type_enc::INT32:
			case type_enc::UINT32: return 'i';
			case  type_enc::INT64:
			case  type_enc::UINT64: return 'l';
			case  type_enc::FLOAT32: return 'f';
			case type_enc::FLOAT64: return 'd';
			case type_enc::FLOAT80: return 'D';
			}
			return 0;
		}
		static int bytes_len(unsigned char enct) {
			switch (enct) {
			case type_enc::INT32:   return sizeof(int);
			case type_enc::INT64:   return sizeof(long long);
			case type_enc::UINT16:  return sizeof(unsigned short);
			case type_enc::UINT8:   return sizeof(unsigned char);
			case type_enc::UINT32:  return sizeof(unsigned int);
			case type_enc::UINT64:  return sizeof(unsigned long long);
			case type_enc::INT8:    return sizeof(char);
			case type_enc::INT16:   return sizeof(short);
			case type_enc::FLOAT32: return sizeof(float);
			case type_enc::FLOAT64: return sizeof(double);
			case type_enc::FLOAT80: return sizeof(long double);
			}
			return 0;
		}
		static bool is_signed(unsigned char c) {
			return c == FLOAT32 || c == FLOAT64 || c == FLOAT80 || c & SIGNED;
		}
	};
	void *data;
	int len;		// type len (not bytes len)

	int n_time_channels() const { return n_time; };
	int n_val_channels() const { return n_val; };

	// time unit & unitless time
	int time_unit() const { return _time_unit; };
	int Time_ch_vec(int idx, int chan, void * data_) const { return Time<int>(idx, chan, data); } // Time(idx,chan)
	void SetVal_ch_vec(int idx, int chan, float _val, void *data_) { setVal(idx, chan, _val, data_); };

	inline void Set(int idx, int *times, float *vals) { Set(idx, times, vals, data); }

	inline void Set(int idx, int *times, float *vals, void* data_)
	{
		for (int chan = 0; chan < n_time; ++chan)
			setTime(idx, chan, times[chan], (char*)data_);

		for (int chan = 0; chan < n_val; ++chan)
			setVal(idx, chan, vals[chan], (char*)data_);
	}

	size_t size() const { return struct_size; }

	void init(const SignalInfo &info) {
		if (sid == info.sid) return;
		_time_unit = info.time_unit;
		time_channel_offsets = info.time_channel_offsets;
		val_channel_offsets = info.val_channel_offsets;
		time_channel_types = info.time_channel_types;
		val_channel_types = info.val_channel_types;

		struct_size = info.bytes_len;
		n_time = info.n_time_channels;
		n_val = info.n_val_channels;

		sid = info.sid;
	}

	int struct_size;
	int n_time;
	int n_val;

	std::array<int, GENERIC_SIG_VEC_MAX_CHANNELS> time_channel_offsets;
	std::array<int, GENERIC_SIG_VEC_MAX_CHANNELS> val_channel_offsets;
	std::array<unsigned char, GENERIC_SIG_VEC_MAX_CHANNELS> time_channel_types;
	std::array<unsigned char, GENERIC_SIG_VEC_MAX_CHANNELS> val_channel_types;

	int sid;

	void set_data(void* _data, int _len) {
		data = _data;
		len = _len;
	}
	GenericSigVec() : data(nullptr), sid(-1), len(0), struct_size(0), n_time(0), n_val(0) { time_channel_offsets.fill(0); val_channel_offsets.fill(0); time_channel_types.fill(0); val_channel_types.fill(0); }
	GenericSigVec(const string& signalSpec, int time_unit = MedTime::Undefined) : GenericSigVec() { _time_unit = time_unit; init_from_spec(signalSpec); }
	GenericSigVec(SigType sigtype, int time_unit = MedTime::Undefined) : GenericSigVec() { _time_unit = time_unit; init_from_sigtype(sigtype); }
	GenericSigVec(const GenericSigVec& other) { *this = other; }

	void copy_signal_metadata(const GenericSigVec& other) {
		struct_size = other.struct_size;
		sid = other.sid;
		n_time = other.n_time;
		n_val = other.n_val;
		_time_unit = other._time_unit;
		time_channel_offsets = other.time_channel_offsets;
		time_channel_types = other.time_channel_types;
		val_channel_offsets = other.val_channel_offsets;
		val_channel_types = other.val_channel_types;
	}

	GenericSigVec& operator=(const GenericSigVec& other) {
		data = other.data;
		len = other.len;
		copy_signal_metadata(other);
		return *this;
	}

	void init_from_spec(const string& signalSpec);
	void init_from_sigtype(SigType sigtype);
	void init_from_repo(MedRepository& repo, int sid);

	template<typename T = int>
	T Time(int idx, int chan) const { return Time<T>(idx, chan, data); }

	template<typename T = int>
	T Time(int idx, int chan, const void* data_) const {
		auto field_ptr = ((const char*)data_) + idx * struct_size + time_channel_offsets[chan];
		switch (time_channel_types[chan]) {
		case type_enc::INT32:   return (T)(*(int*)(field_ptr));
		case type_enc::INT64:   return (T)(*(long long*)(field_ptr));
		case type_enc::UINT16:  return (T)(*(unsigned short*)(field_ptr));
		case type_enc::UINT8:   return (T)(*(unsigned char*)(field_ptr));
		case type_enc::UINT32:  return (T)(*(unsigned int*)(field_ptr));
		case type_enc::UINT64:  return (T)(*(unsigned long long*)(field_ptr));
		case type_enc::INT8:    return (T)(*(char*)(field_ptr));
		case type_enc::INT16:   return (T)(*(short*)(field_ptr));
		case type_enc::FLOAT32: return (T)(*(float*)(field_ptr));
		case type_enc::FLOAT64: return (T)(*(double*)(field_ptr));
		case type_enc::FLOAT80: return (T)(*(long double*)(field_ptr));
		}
		return 0;
	}

	template<typename T = int>
	void setTime(int idx, int chan, T new_val) { setTime<T>(idx, chan, new_val, data); }

	template<typename T = int>
	void setTime(int idx, int chan, T new_val, void* data_) {
		auto field_ptr = ((char*)data_) + idx * struct_size + time_channel_offsets[chan];
		switch (time_channel_types[chan]) {
		case type_enc::INT32:   (*(int*)(field_ptr)) = (int)new_val; return;
		case type_enc::INT64:   (*(long long*)(field_ptr)) = (long long)new_val; return;
		case type_enc::UINT16:  (*(unsigned short*)(field_ptr)) = (unsigned short)new_val; return;
		case type_enc::UINT8:   (*(unsigned char*)(field_ptr)) = (unsigned char)new_val; return;
		case type_enc::UINT32:  (*(unsigned int*)(field_ptr)) = (unsigned int)new_val; return;
		case type_enc::UINT64:  (*(unsigned long long*)(field_ptr)) = (unsigned long long)new_val; return;
		case type_enc::INT8:    (*(char*)(field_ptr)) = (char)new_val; return;
		case type_enc::INT16:   (*(short*)(field_ptr)) = (short)new_val; return;
		case type_enc::FLOAT32: (*(float*)(field_ptr)) = (float)new_val; return;
		case type_enc::FLOAT64: (*(double*)(field_ptr)) = (double)new_val; return;
		case type_enc::FLOAT80: (*(long double*)(field_ptr)) = (long double)new_val; return;
		}
	}

	template<typename T = float>
	T Val(int idx, int chan) const { return Val<T>(idx, chan, data); }

	template<typename T = float>
	T Val(int idx) const { return Val<T>(idx, 0, data); }

	template<typename T = float>
	T Val(int idx, int chan, const void* data_) const {
		auto field_ptr = ((char*)data_) + idx * struct_size + val_channel_offsets[chan];
		switch (val_channel_types[chan]) {
		case type_enc::FLOAT32: return (T)(*(float*)(field_ptr));
		case type_enc::INT16:   return (T)(*(short*)(field_ptr));
		case type_enc::UINT16:  return (T)(*(unsigned short*)(field_ptr));
		case type_enc::UINT8:   return (T)(*(unsigned char*)(field_ptr));
		case type_enc::UINT32:  return (T)(*(unsigned int*)(field_ptr));
		case type_enc::UINT64:  return (T)(*(unsigned long long*)(field_ptr));
		case type_enc::INT8:    return (T)(*(char*)(field_ptr));
		case type_enc::INT32:   return (T)(*(int*)(field_ptr));
		case type_enc::INT64:   return (T)(*(long long*)(field_ptr));
		case type_enc::FLOAT64: return (T)(*(double*)(field_ptr));
		case type_enc::FLOAT80: return (T)(*(long double*)(field_ptr));
		}
		return 0;
	}

	template<typename T = float>
	void setVal(int idx, int chan, T new_val) { setVal<T>(idx, chan, new_val, data); }

	template<typename T = float>
	void setVal(int idx, int chan, T new_val, void* data_) {
		auto field_ptr = ((char*)data_) + idx * struct_size + val_channel_offsets[chan];
		switch (val_channel_types[chan]) {
		case type_enc::FLOAT32: (*(float*)(field_ptr)) = (float)new_val; return;
		case type_enc::INT16:   (*(short*)(field_ptr)) = (short)new_val; return;
		case type_enc::UINT16:  (*(unsigned short*)(field_ptr)) = (unsigned short)new_val; return;
		case type_enc::UINT8:   (*(unsigned char*)(field_ptr)) = (unsigned char)new_val; return;
		case type_enc::UINT32:  (*(unsigned int*)(field_ptr)) = (unsigned int)new_val; return;
		case type_enc::UINT64:  (*(unsigned long long*)(field_ptr)) = (unsigned long long)new_val; return;
		case type_enc::INT8:    (*(char*)(field_ptr)) = (char)new_val; return;
		case type_enc::INT32:   (*(int*)(field_ptr)) = (int)new_val; return;
		case type_enc::INT64:   (*(long long*)(field_ptr)) = (long long)new_val; return;
		case type_enc::FLOAT64: (*(double*)(field_ptr)) = (double)new_val; return;
		case type_enc::FLOAT80: (*(long double*)(field_ptr)) = (long double)new_val; return;
		}
	}

	// channel 0 easy API
	inline int Time(int idx) const { return Time(idx, 0); }
	inline float Val(int idx) const { return Val(idx, 0); }

	inline int TimeU(int idx, int to_time_unit) const { return med_time_converter.convert_times(time_unit(), to_time_unit, Time(idx)); }
	inline int Date(int idx) const { return med_time_converter.convert_times(time_unit(), MedTime::Date, Time(idx)); }
	inline int Years(int idx) const { return med_time_converter.convert_times(time_unit(), MedTime::Years, Time(idx)); }
	inline int Months(int idx) const { return med_time_converter.convert_times(time_unit(), MedTime::Months, Time(idx)); }
	inline int Days(int idx) const { return med_time_converter.convert_times(time_unit(), MedTime::Days, Time(idx)); }
	inline int Hours(int idx) const { return med_time_converter.convert_times(time_unit(), MedTime::Hours, Time(idx)); }
	inline int Minutes(int idx) const { return med_time_converter.convert_times(time_unit(), MedTime::Minutes, Time(idx)); }

	// general channel API
	inline int TimeU(int idx, int chan, int to_time_unit) const { return med_time_converter.convert_times(time_unit(), to_time_unit, Time(idx, chan)); }
	inline int Date(int idx, int chan) const { return med_time_converter.convert_times(time_unit(), MedTime::Date, Time(idx, chan)); }
	inline int Years(int idx, int chan) const { return med_time_converter.convert_times(time_unit(), MedTime::Years, Time(idx, chan)); }
	inline int Months(int idx, int chan) const { return med_time_converter.convert_times(time_unit(), MedTime::Months, Time(idx, chan)); }
	inline int Days(int idx, int chan) const { return med_time_converter.convert_times(time_unit(), MedTime::Days, Time(idx, chan)); }
	inline int Hours(int idx, int chan) const { return med_time_converter.convert_times(time_unit(), MedTime::Hours, Time(idx, chan)); }
	inline int Minutes(int idx, int chan) const { return med_time_converter.convert_times(time_unit(), MedTime::Minutes, Time(idx, chan)); }

	template <class S> void set_funcs() const { /* Do nothing */ }

	SigType get_type() const { return type; }

	// helper functions for common operations
	//----------------------------------------
	// returns the first index i in the usv that has Time(i, time_chan) > time_bound, if none : return -1
	int get_index_gt_time_bound(int time_chan, int time_bound);
	// returns the first index i in the usv that has Time(i, time_chan) >= time_bound, if none : return -1
	int get_index_ge_time_bound(int time_chan, int time_bound);

	bool compareTimeLt(const void* data1, int idx1, const void* data2, int idx2) const {
		for (int tchan = 0; tchan < n_time; tchan++) {
			if (this->Time(idx1, tchan, data1) > this->Time(idx2, tchan, data2))
				return false;
			if (this->Time(idx1, tchan, data1) < this->Time(idx2, tchan, data2))
				return true;
		}
		for (int vchan = 0; vchan < n_val; vchan++) {
			if (this->Val(idx1, vchan, data1) > this->Val(idx2, vchan, data2))
				return false;
			if (this->Val(idx1, vchan, data1) < this->Val(idx2, vchan, data2))
				return true;
		}
		return false;
	}

	int compare_elements(const void* elem1, const void* elem2) {
		if (compareTimeLt(elem1, 0, elem2, 0)) return - 1;
		return 0;
	}


	bool compare_gsv_elements(const GSVElement &e1, const GSVElement &e2) {
		return compareTimeLt(&e1.data[0], 0, &e2.data[0], 0);
	}

	int inplace_sort_data(const void *data, int nelem) {
		vector<GSVElement> dc(nelem);
		unsigned char* p_data = (unsigned char *)data;
		int k = 0;
		for (int i = 0; i < nelem; i++) {
			dc[i].data.resize(this->struct_size);
			for (int j = 0; j < this->struct_size; j++)
				dc[i].data[j] = p_data[k++];
		}
		
		std::sort(dc.begin(), dc.end(), [this](GSVElement l, GSVElement r) {return compare_gsv_elements(l, r); });

		k = 0;
		for (int i = 0; i < nelem; i++) {
			dc[i].data.resize(this->struct_size);
			for (int j = 0; j < this->struct_size; j++)
				p_data[k++] = dc[i].data[j];
		}

		return 0;
		//qsort(data, nelem, this->struct_size, &this->compare_elements);
	}

	bool compareData(int idx, const GenericSigVec& other_gsv, int other_idx) const {
		return std::memcmp(((const char*)data) + idx * struct_size, ((const char*)other_gsv.data) + other_idx * struct_size, struct_size) == 0;
	}

	static string get_type_generic_spec(SigType t);

	string get_signal_generic_spec() const;

protected:
	const SigType type = T_Generic; // type of the generic signal
	int _time_unit = MedTime::Undefined;
};

/**
* Managed memory version of UniversalSigVec.
* for example when allocating new virtual records
*/
class GenericSigVec_mem : public GenericSigVec {
public:
	bool manage; /// if true manages the memory
	string signal_spec; ///< stores the signal spec
	GenericSigVec_mem() {
		manage = false;
	}

	~GenericSigVec_mem() {
		if (len > 0 && data != NULL && manage) {
			delete[](char *)data;
			data = NULL;
			len = 0;
		}
	}

	void set(const GenericSigVec &s) {
		manage = false;
		GenericSigVec::operator=(s);
	}
};

#endif
