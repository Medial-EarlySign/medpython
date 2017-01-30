//condor_runner : An object for managing condor submissions

#define _CRT_SECURE_NO_WARNINGS

#pragma once

#include "medial_utilities/medial_utilities/medial_utilities.h"

#include <stdio.h>
#include <stdlib.h>

#include <string>
#include <vector> 
#include <map>
#include <algorithm> 
#include <sstream> 
#include <iostream> 

#define MAX_BUFFER_SIZE 2048

using namespace std ; 

class command {
public:
	map<string,string> entries ;
	double w ;
} ;

struct compare_commands {
    bool operator()(const command& left, const command& right) {
                return (left.w > right.w) ;
    }
} ;


class condor_runner {
private:
	map<string,string> entries ;
	vector<command> commands ;
	string file_name ;
	int stage ; // 0 - Initial, 1 - Printed, 2 - Submitted

public:
	condor_runner() ;
	condor_runner(string file_name) ;

	void clear() ;

	void set_field(string key, string value) ;
	void add_command()  ;
	void add_command(command& cmd) ;
	void add_command(double w) ;

	int set_weight(double w) ;
	int set_weight(int i, double w) ;
	int set_command_field(string key, string value) ;
	int set_command_field(int i, string key, string value) ;

	void order_commands() ;

	void set_file_name(string file_name) ;
	int print(string file_name) ;
	int print() ;
	
	int print_after(string file_name,string header_name) ;
	int print_after(string header_name) ;

	int submit() ;
	int wait() ;
	int wait(string log_file) ;
};

