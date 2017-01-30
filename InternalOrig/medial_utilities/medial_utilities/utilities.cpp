//General functions

#define _CRT_SECURE_NO_WARNINGS

#include "medial_utilities.h"
using namespace std;

// File utilities
bool file_exists (const char *filename) {
	FILE *fp = fopen(filename,"r");
	if (fp != NULL) {
		fclose(fp) ;
		return true ;
	}

	return false ;
}

FILE* safe_fopen(const char* filename, const char* mode, bool exit_on_fail) {
	string fn;
	fix_path(string(filename), fn);
	FILE* fp = fopen(fn.c_str(), mode);
	if (fp == NULL) {
		fprintf(stderr, "Failed to open file %s in mode %s, error code is %d (%s)\n", filename, mode, errno, strerror(errno));
		if (exit_on_fail) exit(EXIT_FAILURE);
	}

	return(fp);
}

gzFile safe_gzopen(const char* filename, const char* mode, bool print_msg, bool exit_on_fail) {
	string fn;
	fix_path(string(filename), fn);

	string fm = mode;
	if (fm == "r") fm = "rb";
	if (fm.substr(0, 1) == "w" && (fm.length() == 1 || fm.substr(1,1) != "b")) {
		fm = (fm.length() == 1) ? "wb" : ("wb" + fm.substr(1)); // handle correctly "w", "wb", "wb8", "wT", "wbH" 
	}

	// prefer opening for reading a '.gz' version
	if (fm.substr(0, 1) == "r" && (file_exists((fn + ".gz").c_str()))) {
		fn += ".gz";
	}

	if (print_msg) {
		fprintf(stderr, "Going to open file %s in mode %s\n", fn.c_str(), fm.c_str());
		fflush(stderr);
	}
	gzFile fp = gzopen(fn.c_str(), fm.c_str());
	if (fp == NULL) {
		int err;
        const char * error_string;
        error_string = gzerror(fp, &err);
		fprintf(stderr, "Failed to gzopen file %s in mode %s, error code is %d (\'%s\')\n", fn.c_str(), fm.c_str(), err, error_string);
		fflush(stderr);
		if (exit_on_fail) exit(EXIT_FAILURE);
	}

	return(fp);
}
		
char * gzGetLine(gzFile file, string& str) {
	char txtline[MAX_STRING_LEN];

	char * res;
	if ((res = gzgets(file, txtline, MAX_STRING_LEN)) == NULL) {
		return(NULL);
	}

	size_t txtlen = strnlen(txtline, MAX_STRING_LEN);
	assert (txtlen > 0);
	assert (txtline[txtlen - 1] == '\n'); // verify full line
	if (txtlen >= 2 && txtline[txtlen - 2] == '\r') txtlen --; // remove CR

	str = string(txtline, txtline + (txtlen - 1));

	return(res);
}

// Time Utilities

int days2month[] = {0,31,59,90,120,151,181,212,243,273,304,334,365} ;

// hours from 01/01/1900
double get_hour(const char *time, int format) {
			
	int yr,mon,dy,hr,min  ;

	int cnt ;
	if (format == 0) {
		cnt = sscanf(time,"%d/%d/%d %d:%d",&dy,&mon,&yr,&hr,&min) ;
	} else if (format == 1) {
		cnt = sscanf( time, "%d-%d-%d %d:%d", &yr, &mon, &dy, &hr, &min);
	} else {
		fprintf(stderr,"unknown format %d\n",format) ;
		return -1 ;
	}

	if (cnt != 5 && cnt != 3)
		return -1.0 ;

	if (cnt == 3)
		hr=min=0 ;

	if (mon < 1 || mon > 12 || yr < 1900 || dy < 0 || hr < 0 || hr > 23 || min < 0 || min > 59)
		return -1.0 ;

	// Full years
	int days = 365 * (yr-1900) ;
	days += (yr-1897)/4 ;
	days -= (yr-1801)/100 ;
	days += (yr-1601)/400 ;

	// Full Months
	days += days2month[mon-1] ;
	if (mon>2 && (yr%4)==0 && ((yr%100)!=0 || (yr%400)==0))
		days ++ ;
	days += (dy-1) ;

	double hours = days*24 + hr + min/60.0 ;

	return ((double) hours) ;
}

void hours2time (double hours, char *time) {

	int year = 1900 + ((int)hours)/24/365 - 10 ;
	sprintf(time,"1/1/%d 0:0",year+1) ;

	while (get_hour(time) < hours) {
		year ++ ;
		sprintf(time,"1/1/%d 0:0",year+1) ;
	}

	int month = 1 ;
	sprintf(time,"1/%d/%d 0:0",month+1,year) ;

	while (month <= 12 &&get_hour(time) < hours) {
		month ++ ;
		sprintf(time,"1/%d/%d 0:0",month+1,year) ;
	}

	sprintf(time,"1/%d/%d 0:0",month,year) ;
	hours -= get_hour(time) ;

	int day = (int) hours/24 ;
	hours -= 24.0*day ;

	sprintf(time,"%d_%d_%d_%d_%d",day,month,year,(int) hours/60, ((int) hours)%60) ;
	return ;
}

// days from 01/01/1900
int get_day(int val) {
		
	int year = val/100/100 ;
	int month = (val/100)%100 ;
	int day = val%100 ;

	if (month < 1 || month > 12 || year < 1900 || day < 0)
		return -1;

	// Full years
	int days = 365 * (year-1900) ;
	days += (year-1897)/4 ;
	days -= (year-1801)/100 ;
	days += (year-1601)/400 ;

	// Full Months
	days += days2month[month-1] ;
	if (month>2 && (year%4)==0 && ((year%100)!=0 || (year%400)==0))
		days ++ ;
	days += (day-1) ;

	return days;
}

double get_day(const char *time, int format) {
			
	int yr,mon,dy  ;

	int cnt ;
	if (format == 0) {
		cnt = sscanf(time,"%d/%d/%d",&dy,&mon,&yr) ;
	} else if (format == 1) {
		cnt = sscanf( time, "%d-%d-%d", &yr, &mon, &dy);
	} else {
		fprintf(stderr,"unknown format %d\n",format) ;
		return -1 ;
	}

	if (cnt != 3)
		return -1.0 ;

	if (mon < 1 || mon > 12 || yr < 1900 || dy < 0)
		return -1.0 ;

	// Full years
	int days = 365 * (yr-1900) ;
	days += (yr-1897)/4 ;
	days -= (yr-1801)/100 ;
	days += (yr-1601)/400 ;

	// Full Months
	days += days2month[mon-1] ;
	if (mon>2 && (yr%4)==0 && ((yr%100)!=0 || (yr%400)==0))
		days ++ ;

	days += (dy-1) ;

	return ((double) days) ;
}

// Days -> Date
int get_date(int days) {

	// Full Years
	int year = 1900 + days/365 ;
	days %= 365 ;

	days -= (year-1897)/4 ;
	days += (year-1801)/100 ;
	days -= (year-1601)/400 ;

	if (days < 0) {
		year -- ;
		days += 365 ;
		if ((year%4)==0 && ((year%100)!=0 || (year%400)==0)) {
			days ++ ;
			if (days == 366) {
				days = 0 ;
				year ++ ;
			}
		}
	}

	// Full Months
	bool leap_year = ((year%4)==0 && ((year%100)!=0 || (year%400)==0)) ;
	int month ;
	for (int i = 1; i <= 12; i ++) {
		int mdays = days2month[i] + ((leap_year && i > 1) ? 1 : 0) ;
		if (days < mdays) {
			month = i ;
			days -= (days2month[i-1] + ((leap_year && i > 2) ? 1 : 0)) ;
			break ;
		}
	}

	days ++ ;
	return days + 100*month + 10000*year ;
}


// minutes from 01/01/2000
double get_min(char *time, int format) {
			
	int yr,mon,dy,hr,min  ;

	int cnt ;
	if (format == 0) {
		cnt = sscanf(time,"%d/%d/%d %d:%d",&dy,&mon,&yr,&hr,&min) ;
	} else if (format == 1) {
		cnt = sscanf( time, "%d-%d-%d %d:%d", &yr, &mon, &dy, &hr, &min);
	} else {
		fprintf(stderr,"unknown format %d\n",format) ;
		return -1 ;
	}

	if (cnt != 5)
		return -1.0 ;

	if (mon < 1 || mon > 12 || yr < 2000 || dy < 0 || hr < 0 || hr > 23 || min<0 || min>59)
		return -1.0 ;

	// Full years
	int days = 365 * (yr-2000) ;
	days += (yr-1997)/4 ;
	days -= (yr-1901)/100 ;
	days += (yr-1601)/400 ;

	// Full Months
	days += days2month[mon-1] ;
	if (mon>2 && (yr%4)==0 && ((yr%100)!=0 || (yr%400)==0))
		days ++ ;

	days += (dy-1) ;

	double mins = 60*(days*24 + hr) + min ;

	return ((double) mins) ;
}

double min2day (double mins) {
	return (mins/60.0/24.0 + get_day(string("01/01/2000").c_str())) ;
}

int fix_path(const string& in, string& out) {
	// fprintf(stderr, "Converting path \'%s\'\n", in.c_str());
	// fflush(stderr);

	map<string, string> folders ;
	folders["W"] = "Work" ;
	folders["U"] = "UsersData" ;
	folders["X"] = "Temp" ;
	folders["P"] = "Products" ;
	folders["T"] = "Data" ;
	// fprintf(stderr, "Initialized network drive table\n");

#ifndef _WIN32
	// on Linux, handle first the Windows native format: \\\\nas1\\Work\..  
	if (in.length() >= 2 && in.substr(0, 2) =="\\\\") {
		// just switching '\' to '/'; works, but adjacent slashes should be unified
		out = in;
		char revSlash = '\\';
		char fwdSlash = '/';
		std::replace(out.begin(), out.end(), revSlash, fwdSlash);
		fprintf(stderr, "Converted path \'%s\' to \'%s\'\n", in.c_str(), out.c_str());
		fflush(stderr);

		return 0;
	}
#endif

	// Work only on "X:/...."  or  "/cygdrive/X/...." input strings
	if ((in.length() < 3 || in.substr(1,2) != ":/") && 
		(in.length() < 12 || in.substr(0, 10) != "/cygdrive/" || in.substr(11, 1) != "/")) {
		out = in ;
		return 0 ;
	}

	char driveLetter = (in.substr(1, 2) == ":/") ? in.substr(0, 1)[0] : (char)toupper(in.substr(10, 1)[0]);
	string drive = string(1, driveLetter);

	// Special handling of local drives (C: or D:)
	if (drive == "C" || drive == "D" || drive == "H") {
		out = in ;
		return 0 ;
	}

	int pathPos = (in.substr(1, 2) == ":/") ? 3 : 12;

	if (folders.find(drive) == folders.end()) {
		fprintf(stderr, "Unknown Folder Map %s. Trying to work with original path\n", drive.c_str()) ;
		out = in;
		return 0 ;
	}

#ifdef _WIN32
	out = "\\\\nas1\\" + folders[drive] ;
#else
	out = "/nas1/" + folders[drive] ;
#endif

	istringstream in_stream(in.substr(pathPos, in.length() - pathPos)) ;	
	string token ;

	while (getline(in_stream,token,'/'))
#ifdef _WIN32
		out += "\\" + token ;
#else
		out += "/" + token ;
#endif
	fprintf(stderr, "Converted path \'%s\' to \'%s\'\n", in.c_str(), out.c_str());
	fflush(stderr);

	return 0 ;

}
