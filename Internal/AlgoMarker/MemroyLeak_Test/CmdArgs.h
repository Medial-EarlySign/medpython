#ifndef __CMD_ARGS__
#define __CMD_ARGS__
#include <MedUtils/MedUtils/MedUtils.h>
#include <fstream>
#include <MedProcessTools/MedProcessTools/MedSamples.h>

#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL LOG_DEF_LEVEL

class ProgramArgs : public medial::io::ProgramArgs_base {
private:
	void post_process() {
	}
public:
	string amconfig;
	string amlib;
	string in_jsons;
	string jreq;
	int prediction_time;
	int pid_id;

	int test_stage;

	bool calc_by_type;

	ProgramArgs() {
		po::options_description options("Required options");
		options.add_options()
			("amconfig", po::value<string>(&amconfig)->required(), "ampath")
			("amlib", po::value<string>(&amlib)->default_value(""), "amlib")
			("in_jsons", po::value<string>(&in_jsons)->required(), "in_jsons")
			("jreq", po::value<string>(&jreq)->default_value(""), "jreq")
			("pid_id", po::value<int>(&pid_id)->default_value(1), "pid id")
			("prediction_time", po::value<int>(&prediction_time)->default_value(20210101), "prediction_time")
			("test_stage", po::value<int>(&test_stage)->default_value(0), "test stage")
			("calc_by_type", po::value<bool>(&calc_by_type)->default_value(false), "calc_by_type")

			;

		init(options);
	}
};


#endif
