#include "condor_runner.h"

condor_runner::condor_runner() {
	this->commands.clear() ;
	this->entries.clear() ;
	this->file_name.clear() ;
	this->stage = 0 ;
}
	
condor_runner::condor_runner(string file_name) {
	this->commands.clear() ;
	this->entries.clear() ;
	this->file_name = file_name ;
	this->stage = 0 ;
}

void condor_runner::clear() {
	entries.clear() ;
	commands.clear() ;
	file_name.clear() ;
	stage = 0 ;
}

void condor_runner::set_field(string key, string value) {
	this->entries[key] = value ;
}

void condor_runner::add_command() {
	command new_command ;
	new_command.w = 0.0 ;
	this->commands.push_back(new_command) ;
}
	
void condor_runner::add_command(command& cmd) {
	this->commands.push_back(cmd) ;
}

void condor_runner::add_command(double w) {
	command new_command ;
	new_command.w = w ;
	this->commands.push_back(new_command) ;
}

int condor_runner::set_weight(double w) {

	if (this->commands.empty()) {
		fprintf(stderr,"Cannot set weight on a non-existing command\n") ;
		return -1 ;
	}

	this->commands.back().w = w ;
	return 0 ;
}

int condor_runner::set_weight(int i, double w) {

	if (this->commands.size() <= i) {
		fprintf(stderr,"Cannot set weight on a non-existing command\n") ;
		return -1 ;
	}

	this->commands[i].w = w ;
	return 0 ;
}

int condor_runner::set_command_field(string key, string value) {

	if (this->commands.empty()) {
		fprintf(stderr,"Cannot set field on a non-existing command\n") ;
		return -1 ;
	}

	(this->commands.back().entries)[key] = value ;
	return 0 ;
}

int condor_runner::set_command_field(int i, string key, string value) {
	
	if (this->commands.size() <= i) {
		fprintf(stderr,"Cannot set field on a non-existing command\n") ;
		return -1 ;
	}

	this->commands[i].entries[key] = value ;
	return 0 ;
}

void condor_runner::order_commands() {
	sort(this->commands.begin(),this->commands.end(),compare_commands()) ;
}

void condor_runner::set_file_name(string file_name) {
	this->file_name = file_name ;
}

int condor_runner::print(string file_name) {
	this->set_file_name(file_name) ;
	return this->print() ;
}

int condor_runner::print_after(string file_name, string header_name) {
	this->set_file_name(file_name) ;
	return this->print_after(header_name) ;
}

int condor_runner::print_after(string header_name) {

	if (this->file_name.empty()) {
		fprintf(stderr,"Cannot print without setting file-name\n") ;
		return -1 ;
	}

	FILE *ifp = safe_fopen(header_name.c_str(), "r", false) ;
	if (ifp==NULL) {
		fprintf(stderr,"Cannot open header file %s for reading\n",header_name.c_str()) ;
		return -1 ;
	}

	FILE *ofp = safe_fopen(this->file_name.c_str(), "w", false) ;
	if (ofp==NULL) {
		fprintf(stderr,"Cannot open submission file %s for writing\n",this->file_name.c_str()) ;
		return -1 ;
	}	

	char buffer[MAX_BUFFER_SIZE] ;
	while (!feof(ifp)) {
		fgets(buffer,MAX_BUFFER_SIZE,ifp) ;

		if (ferror(ifp)) {
			fprintf(stderr,"Cannot read from header file %s\n",header_name.c_str()) ;
			return -1 ;
		}
		if (feof(ifp))
			break ;

		fprintf(ofp,buffer) ;
	}

	fclose(ifp) ;
	fclose(ofp) ;

	this->stage = 1 ;
	return this->print() ;
}

int condor_runner::print() {

	if (this->file_name.empty()) {
		fprintf(stderr,"Cannot print without setting file-name\n") ;
		return -1 ;
	}

	FILE *fp ;

	if (this->stage == 2) {
		fprintf(stderr,"Cannot change sumbission file after submission\n") ;
		return -1 ;
	} else if (this->stage == 1)
		fp = safe_fopen(this->file_name.c_str(), "a", false) ;
	else
		fp = safe_fopen(this->file_name.c_str(), "w", false) ;

	if (fp == NULL) {
		fprintf(stderr,"Canont open submission file \'%s\'\n",this->file_name.c_str()) ;
		return -1 ;
	}

	for (map<string,string>::iterator it=this->entries.begin(); it!=this->entries.end(); it++)
		fprintf(fp,"%s= %s\n",it->first.c_str(),it->second.c_str()) ;

	for (int icmd=0; icmd < this->commands.size(); icmd++) {
		fprintf(fp,"\n") ;
		for (map<string,string>::iterator cit=this->commands[icmd].entries.begin(); cit!=this->commands[icmd].entries.end(); cit++)
			fprintf(fp,"%s= %s\n",cit->first.c_str(),cit->second.c_str()) ;
		fprintf(fp,"queue\n") ;
	}

	this->stage = 1 ;
	fclose(fp) ;

	return 0 ;
}
		
int condor_runner::submit() {

	if (this->stage != 1)  {
		fprintf(stderr,"Canont submit at stage %d\n",this->stage) ;
		return -1 ;
	}

	char command[1024] ;
	sprintf(command,"condor_submit %s",this->file_name.c_str()) ;

	int rc ;
	if ((rc = system(command)) != 0) {
		fprintf(stderr,"Submission failed (\'%s\', rc = %d)\n",command,rc) ;
		return -1 ;
	}

	this->stage = 2 ;
	return 0 ;
}

int condor_runner::wait() {

	if (this->stage != 2)  {
		fprintf(stderr,"Canont wait at stage %d\n",this->stage) ;
		return -1 ;
	}

	char command[1024] ;
	map<string,int> log_files ;

	string initial_dir = "" ;
	if (this->entries.find("initialdir") != this->entries.end())
		initial_dir = this->entries["initialdir"] ;

	if (this->entries.find("log") != this->entries.end()) {
		string log_file = initial_dir + "/" + this->entries["log"] ;
		log_files[log_file] = 1 ;
	}

	for (int icmd=0; icmd < this->commands.size(); icmd++) {
		if (this->commands[icmd].entries.find("log") != this->commands[icmd].entries.end()) {
			string log_file = initial_dir + "/" + this->commands[icmd].entries["log"] ;
			log_files[log_file] = 1 ;
		}
	}

	for (map<string,int>::iterator it=log_files.begin(); it != log_files.end(); it++) {
		sprintf(command,"condor_wait %s",it->first.c_str()) ;
		fprintf(stderr,"Waiting for %s\n",it->first.c_str()) ;

		if (system(command) != 0) {
			fprintf(stderr,"Waiting failed\n") ;
			return -1 ;
		}
	}

	return 0 ;
}

int condor_runner::wait(string log_file) {

	if (this->stage != 2)  {
		fprintf(stderr,"Canont wait at stage %d\n",this->stage) ;
		return -1 ;
	}

	char command[1024] ;
	sprintf(command,"condor_wait %s",log_file.c_str()) ;
	fprintf(stderr,"Waiting for %s\n",log_file.c_str()) ;

	if (system(command) != 0) {
		fprintf(stderr,"Submission failed\n") ;
		return -1 ;
	}

	return 0 ;
}


