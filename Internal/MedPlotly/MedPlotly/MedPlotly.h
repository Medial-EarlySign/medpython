#pragma once
//
// Class to handle creation of charts using the plotly js package
//
// Main functionalities:
// (1) Configurable page content
// (2) Single signal charts
// (3) Panel charts
// (4) drugs heatmap
// (5) Major demographics table (TBD as a table)
// (6) Registries View (TBD)
// (7) ReadCodes View
//
//



#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <boost/program_options.hpp>
#include <boost/algorithm/string/split.hpp>
#include <InfraMed/InfraMed/InfraMed.h>
#include <InfraMed/InfraMed/MedPidRepository.h>
#include <MedProcessTools/MedProcessTools/MedModel.h>
#include <MedTime/MedTime/MedTime.h>

#include <Logger/Logger/Logger.h>
#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL	LOG_DEF_LEVEL
using namespace std;
namespace po = boost::program_options;


//=========================================================================================================
namespace PlotlyColorDefaults {


	static const vector<string> bg_opaque_colors = { "rgba(162, 217, 206,0.333)" ,  "rgba(199, 210, 232, 0.333)" ,  "rgba(220, 199, 232,0.333)",
		"rgba(232, 199, 204,0.333)",  "rgba(200, 232, 199,0.333)",  "rgba(232, 217, 199,0.333)",  "rgba(215, 216, 192,0.333)" };

	static string get_bg_color(int i) { return PlotlyColorDefaults::bg_opaque_colors[i % (int)PlotlyColorDefaults::bg_opaque_colors.size()]; }

};

//=========================================================================================================
class SignalParams : public SerializableObject {
public:
	int null_zeros = 1;
	int log_scale = 1;
	int time_chan = 0;
	int val_chan = 0;

	int init(map<string, string>& _map);
};

//=========================================================================================================
class PanelInfo : public SerializableObject {
public:
	string name = "";			// to be used when adding it
	string title = "";			// title of panel chart
	vector<string> sigs;		// sigs to be drawn on panel chart.
	vector<string> drugs;		// list of drug groups to visualize
	vector<string> drug_colors;			// bg colors to use. (when not initialized, will be initialized from a defaulted list).
	int null_zeros = -1;		// -1: means take default
	int log_scale = -1;			// -1: means take default
	int width = -1;
	int height = -1;
	int block_mode = -1;

	PanelInfo() {};
	PanelInfo(string _name, string _title, const vector<string> &_sigs, const vector<string> &_drugs, const vector<string> &_drug_colors, int _null_zeros, int _log_scale) {
		name = _name; title = _title; sigs = _sigs; drugs = _drugs; drug_colors = _drug_colors; null_zeros = _null_zeros; log_scale = _log_scale;
	}

	int init(map<string, string>& _map);
};


//=========================================================================================================
class DrugsHeatMapParams : public SerializableObject {
public:
	vector<string> drugs;
	string color0;
	string color1;
	int granularity_months = 1;
	int min_date = 20000101;
	string drug_sig = "Drug";
	int width = 1200;
	int height = 200;


	int init(map<string, string>& _map);
};

//=========================================================================================================
// General Params:
// Either defaulted , or read from a config file
//=========================================================================================================
class MedPlotlyParams : public SerializableObject {
public:
	string page_title = "Patient Report";

	// signal specific params
	map<string, SignalParams> sig_params;

	// Drugs groups
	map<string, vector<string>> drugs_groups;

	// Panels
	map<string, PanelInfo> panels;

	// drugs heatmap
	DrugsHeatMapParams dhm_params;

	// view order
	vector<string> views;

	// javascript
	string js_dir = "";
	vector<string> js_files;

	// zeros handling
	int null_zeros_default = 1; // in default we treat zeros as outliers and don't draw them.
	int log_scale_default = 1;
	int width_default = 600;
	int height_default = 400;
	int block_mode_default = 1;

	// time unit , currently supporting Date and Minutes
	int rep_time_unit = MedTime::Date; // default is Date.

	MedModel model_rep_processors;
	bool load_dynamically = true;
	vector<int> phisical_read_sigs;

	int read_config(const string &config_fname);

};

//=========================================================================================================
// Used in order to be able to mark some interesting times on generated charts
class ChartTimeSign : public SerializableObject {
public:
	int time;
	string name;
	string color = ""; // if empty will use default color lists

	ChartTimeSign() {};
	ChartTimeSign(int t, const string &_name, const string _color) { time = t; name = _name; color = _color; }
	int init(map<string, string>& _map);

	// example of how to initialize a vector:
	//vector<ChartTimeSign> vec = { ChartTimeSign(1,"a","b") , ChartTimeSign(2,"d","e") };

};

class LocalViewsParams {
public:
	LocalViewsParams() {};
	LocalViewsParams(int _pid, int _from, int _to) { pid = _pid; from_date = _from; to_date = _to; }
	int pid = 0;
	int from_date = 0;
	int to_date = 0;


};

class PidDataRec {
public:
	PidRec rec_static;
	PidDynamicRec rec_dynamic;
	bool use_dynamic = false;

	void uget(int sid, UniversalSigVec &usv) {
		if (use_dynamic)
			rec_dynamic.uget(sid, rec_dynamic.get_n_versions() - 1, usv);
		else
			rec_static.uget(sid, usv);
	}
	void uget(const string &sid, UniversalSigVec &usv) {
		if (use_dynamic)
			rec_dynamic.uget(sid, rec_dynamic.get_n_versions() - 1, usv);
		else
			rec_static.uget(sid, usv);
	}
	void *get(const string &sid, int &len) {
		if (use_dynamic)
			return rec_dynamic.get(sid, rec_dynamic.get_n_versions() - 1, len);
		else
			return rec_static.get(sid, len);
	}
	int pid() {
		if (use_dynamic)
			return rec_dynamic.pid;
		else
			return rec_static.pid;
	}
	MedRepository *my_base_rep() {
		if (use_dynamic)
			return rec_dynamic.my_base_rep;
		else
			return rec_static.my_base_rep;
	}
};

//=========================================================================================================
// MedPatientPlotlyDate:
// - Built to handle date based repositories
// - allows configurable viewing of a given PidRec
//=========================================================================================================
class MedPatientPlotlyDate {

public:
	MedPlotlyParams params;

	int read_config(const string &config_name) { return params.read_config(config_name); }

	// gets a pid record and returns a string representing the html for its view page based on configuration
	// mode is either file (for an html file:/// version) or server (for use with a web server returning it)
	int get_rec_html(string &shtml, LocalViewsParams &lvp, PidDataRec &rec, const string &mode, const vector<ChartTimeSign> &sign_times, const vector<string> &view);
	int get_rec_html(string &shtml, LocalViewsParams &lvp, PidDataRec &rec, const string &mode, const vector<ChartTimeSign> &sign_times) {
		return get_rec_html(shtml, lvp, rec, mode, sign_times, params.views);
	}

	int init_rep_processors(MedPidRepository &rep, const string &rep_conf);

private:
	// builders for html
	int add_html_header(string &shtml, const string &mode);
	int add_basic_demographics(string &shtml, PidDataRec &rec, vector<ChartTimeSign> &times);
	int add_panel_chart(string &shtml, LocalViewsParams &lvp, PidDataRec &rec, const PanelInfo &pi, const vector<ChartTimeSign> &times);
	int add_drugs_heatmap(string &shtml, PidDataRec &rec);

	// THIN_RC report
	void add_thin_rc_chart(string &shtml, PidDataRec &rec, const vector<ChartTimeSign> &times);

	// categorical signal , add as table
	void add_categorical_table(string sig, string &shtml, PidDataRec &rec, const vector<ChartTimeSign> &times);

	// heatmap creation
	void get_drugs_heatmap(PidDataRec &rec, vector<int> &_xdates, vector<string> &_sets_names, vector<vector<float>> &_hmap, const vector<string> &drugs);


	// adding a specific dataset of a sepcific sig + channel
	void add_xy_js(string &shtml, vector<int> &dates, vector<float> &vals, int null_zeros_flag, string prefix);
	void add_xy_js(string &shtml, UniversalSigVec &usv, int time_chan, int chan, int null_zeros_flag, string prefix);

	void add_dataset_js(string &html, UniversalSigVec &usv, int time_chan, int chan, int null_zeros_flag, string prefix, string sname, int yaxis, string sig);
	void add_bg_dataset_js(string &shtml, vector<int> &dates, vector<float> &vals, int null_zeros_flag, string color, string prefix, string sname, int yaxis, string name);

	// helpers
	string  date_to_string(int date);
	string	time_to_string(int time, int time_unit);
	string	time_to_string(int time) { return time_to_string(time, params.rep_time_unit); }
	void get_usv_min_max(UniversalSigVec &usv, float &vmin, float &vmax);

};
