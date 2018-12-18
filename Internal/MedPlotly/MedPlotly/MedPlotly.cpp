//
// MedPlotly Implementation
//

#include "MedPlotly.h"
#include <MedUtils/MedUtils/MedGenUtils.h>
#include <MedUtils/MedUtils/MedGlobalRNG.h>

//------------------------------------------------------------------------------------------------
int SignalParams::init(map<string, string>& _map)
{
	for (auto entry : _map) {
		string field = entry.first;
		if (field == "null_zeros") { null_zeros = stoi(entry.second); }
		else if (field == "log_scale") { log_scale = stoi(entry.second); }
		else if (field == "time_chan") { time_chan = stoi(entry.second); }
		else if (field == "val_chan") { val_chan = stoi(entry.second); }

	}

	return 0;
}

//------------------------------------------------------------------------------------------------
int PanelInfo::init(map<string, string>& _map)
{
	for (auto entry : _map) {
		string field = entry.first;
		if (field == "name") name = entry.second;
		else if (field == "title") title = entry.second;
		else if (field == "sigs") split(sigs, entry.second, boost::is_any_of(","));
		else if (field == "drugs") split(drugs, entry.second, boost::is_any_of(","));
		else if (field == "drug_colors") split(drug_colors, entry.second, boost::is_any_of(","));
		else if (field == "null_zeros") null_zeros = stoi(entry.second);
		else if (field == "log_scale") log_scale = stoi(entry.second);
		else if (field == "width") width = stoi(entry.second);
		else if (field == "height") log_scale = stoi(entry.second);
		else if (field == "block_mode") block_mode = stoi(entry.second);
	}

	if (title == "") title = name;

	return 0;
}


//------------------------------------------------------------------------------------------------
int DrugsHeatMapParams::init(map<string, string>& _map)
{
	for (auto entry : _map) {
		string field = entry.first;
		if (field == "granularity") { granularity_months = stoi(entry.second); }
		else if (field == "min_date") { min_date = stoi(entry.second); }
		else if (field == "color0") { color0 = entry.second; }
		else if (field == "color1") { color1 = entry.second; }
		else if (field == "drugs") { drugs.clear(); split(drugs, entry.second, boost::is_any_of(",")); }
	}

	return 0;
}

//------------------------------------------------------------------------------------------------
int ChartTimeSign::init(map<string, string>& _map)
{
	for (auto entry : _map) {
		string field = entry.first;
		if (field == "time") { time = stoi(entry.second); }
		else if (field == "name") { name = stoi(entry.second); }
		else if (field == "color") { color = entry.second; }
	}

	return 0;
}


//------------------------------------------------------------------------------------------------
int MedPlotlyParams::read_config(const string &fname)
{
	ifstream inf(fname);

	if (!inf) {
		MERR("MedSignals: read: Can't open file %s\n", fname.c_str());
		return -1;
	}

	MLOG("MedPlotlyParams :: reading config file %s\n", fname.c_str());
	string curr_line;
	while (getline(inf, curr_line)) {
		if ((curr_line.size() > 1) && (curr_line[0] != '#')) {
			vector<string> fields;
			split(fields, curr_line, boost::is_any_of("\t\n\r"));

			//MLOG("##>> parsing line: %s\n", curr_line.c_str());

			if (fields.size() >= 2) {
				if (fields[0] == "NULL_ZEROS") null_zeros_default = stoi(fields[1]);
				else if (fields[0] == "LOG_SCALE") log_scale_default = stoi(fields[1]);
				else if (fields[0] == "WIDTH") width_default = stoi(fields[1]);
				else if (fields[0] == "HEIGHT") height_default = stoi(fields[1]);
				else if (fields[0] == "BLOCK_MODE") block_mode_default = stoi(fields[1]);
				else if (fields[0] == "SIG") {
					vector<string> flist;
					split(flist, fields[1], boost::is_any_of(","));
					MLOG("%s : %d elements\n", fields[1].c_str(), flist.size());
					SignalParams sp;
					sp.init_from_string(fields[2]);
					for (auto f : flist) {
						sig_params[f] = sp;
						MLOG("sig %s time_chan %d loaded\n", f.c_str(), sig_params[f].time_chan);
					}
				}
				else if (fields[0] == "DRUG_GROUP") {
					vector<string> f;
					split(f, fields[2], boost::is_any_of(","));
					drugs_groups[fields[1]] = f;
					// we also keep a default order of those in dhm_params
					dhm_params.drugs.push_back(fields[1]);
				}
				else if (fields[0] == "PANEL") {
					PanelInfo pi;
					pi.init_from_string(fields[1]);
					panels[pi.name] = pi;
				}
				else if (fields[0] == "DRUGS_HEATMAP") {
					dhm_params.init_from_string(fields[1]);
				}
				else if (fields[0] == "VIEW") {
					vector<string> f;
					split(f, fields[1], boost::is_any_of(","));
					views.insert(views.end(), f.begin(), f.end());
				}
				else if (fields[0] == "JSDIR") js_dir = fields[1];
				else if (fields[0] == "JSFILES") {
					vector<string> f;
					split(f, fields[1], boost::is_any_of(","));
					js_files.insert(views.end(), f.begin(), f.end());
				}
				else if (fields[0] == "TIME_UNIT") {
					rep_time_unit = med_time_converter.string_to_type(fields[1]);
				}
			}
		}
	}

	return 0;
}


//------------------------------------------------------------------------------------------------
int MedPatientPlotlyDate::add_html_header(string &shtml, const string &mode)
{
	string fprefix = "";
	if (mode == "file")
		fprefix = "file:///" + params.js_dir + "/";
	shtml += "<!DOCTYPE HTML>\n<html>\n<head>\n";
	shtml += "\t<link rel=\"stylesheet\" style=\"text/css\" href=\"" + fprefix + "w3/w3.css \" >\n";
	for (auto &f : params.js_files)
		shtml += "\t<script type=\"text/javascript\" src=\"" + fprefix + f + "\"></script>\n";
	shtml += "\t<style> body { margin: 20px;} </style>\n";
	shtml += "</head>\n";
	return 0;
}

//------------------------------------------------------------------------------------------------
int MedPatientPlotlyDate::add_basic_demographics(string &shtml, PidRec &rec, vector<ChartTimeSign> &times)
{
	float age;
	int len;
	int death = 0;
	UniversalSigVec usv;

	MLOG("add_basic_demographics 1\n");
	int byear = (int)(((SVal *)rec.get("BYEAR", len))[0].val);
	int gender = (int)(((SVal *)rec.get("GENDER", len))[0].val);
	rec.uget("DEATH", usv);
	if (usv.len > 0) {
		if (usv.n_time_channels() > 0)
			death = (int)usv.Time(0);
		else
			death = (int)usv.Val(0);
	}

	MLOG("add_basic_demographics 2\n");
	shtml += "<h1> Patient Report </h1>\n";

	shtml += "<h3> pid " + to_string(rec.pid) + " , ";
	if (gender == 1)
		shtml += "Male , ";
	else
		shtml += "Female ,";

	shtml += "Birth Year : " + to_string(byear);
	if (death > 0) {
		string ds = time_to_string(death);
		shtml += " , Death " + ds;
		if (params.rep_time_unit == MedTime::Minutes) death += 1400; // ISSUE bypass death is in days... we move it to end of day
		ChartTimeSign cts(death, "Death", "'red'");
		times.push_back(cts);
	}
	shtml += "</h3>\n";

	MLOG("add_basic_demographics 3\n");
	if (times.size() > 0) {
		shtml += "<h3> Anchor Dates : </h3>\n";
		for (auto &t : times) {
			shtml += "<h3> ";
			if (t.name == "Death")
				age = get_age(med_time_converter.convert_times(params.rep_time_unit, MedTime::Date, t.time), byear);
			else
				age = get_age(t.time, byear);

			stringstream s;
			s << fixed << setprecision(2) << age;
			shtml += "age " + s.str() + " , ";
			if (t.name == "Death")
				shtml += " date: " + time_to_string(t.time);
			else
				shtml += " date: " + date_to_string(t.time);
			if (t.name != "") shtml += " [" + t.name + "] ";
			shtml +="</h3>\n";
		}
	}

	MLOG("add_basic_demographics 4\n");

	return 0;
}

//----------------------------------------------------------------------------------------------------------------------
void MedPatientPlotlyDate::get_drugs_heatmap(PidRec &rec, vector<int> &_xdates, vector<string> &_sets_names, vector<vector<float>> &_hmap, const vector<string> &drugs)
{
	int month_granularity = params.dhm_params.granularity_months; // 1 , 2, 3, 4, 6, 12 (should divide 12)
	string drug_sig = params.dhm_params.drug_sig;
	int hmap_min_date = params.dhm_params.min_date;
	int h_lines = (int)drugs.size();

	if (rec.my_base_rep->sigs.sid(drug_sig) < 0) return; // working with a rep with no Drug signal...

	// read drug data
	UniversalSigVec usv;
	rec.uget(drug_sig, usv);

	vector<int> xdates;
	vector<vector<float>> hmap;

	if (usv.len == 0) return;

	// calculate min/max cells dates, xdays and xdates
	int min_date = usv.Time(0, 0);
	if (min_date < hmap_min_date) min_date = hmap_min_date;
	int max_date = usv.Time(usv.len-1, 0);
	min_date = (min_date/10000)*10000+101;
	max_date = (max_date/10000)*10000+(12-month_granularity)*100 + 1;
	vector<int> xdays;
	xdates.clear();
	int t = min_date;
	while (t<=max_date) {
		int days = med_time_converter.convert_date(MedTime::Days, t);
		xdates.push_back(t);
		xdays.push_back(days);

		t += (month_granularity)*100;
		if (t % 10000 > 1231) t = (t/10000 + 1)*10000 + 101;
	}
	max_date = t;
	xdays.push_back(med_time_converter.convert_date(MedTime::Days, t));

	hmap.clear();
	hmap.resize(h_lines);

	vector<int> drug_days;
	for (int i=0; i<usv.len; i++)
		drug_days.push_back(med_time_converter.convert_date(MedTime::Days, usv.Time(i, 0)));

	int section_id = rec.my_base_rep->dict.section_id(drug_sig);
	vector<int> sets_sum(h_lines, 0);
	int first_nonz = 99999999, last_nonz = 0;
	for (int s=0; s<h_lines; s++) {
		vector<unsigned char> lut;
		string drug_group = drugs[s];
		rec.my_base_rep->dict.prep_sets_indexed_lookup_table(section_id, params.drugs_groups[drug_group], lut);
		hmap[s].resize(xdates.size(), (float)0);
		int last_day_counted = xdays[0]-1;
		int curr_cell = 0;
		for (int i=0; i<usv.len; i++) {
			if (lut[(int)usv.Val(i, 0)]) {
				//MLOG("Taking the drug at i=%d %d,%d\n", i, usv.Time(i,0), (int)usv.Val(i, 1));
				sets_sum[s] += (int)usv.Val(i, 1);
				if (usv.Time(i, 0) < first_nonz && usv.Time(i, 0) >= min_date) first_nonz = usv.Time(i, 0);
				if (usv.Time(i, 0) > last_nonz && usv.Time(i, 0) <= max_date) last_nonz = usv.Time(i, 0);
				int to_day = drug_days[i] + (int)usv.Val(i, 1);
				if (to_day > last_day_counted) {
					int from_day = max(last_day_counted+1, drug_days[i]);
					// get the curr_cell from_day is contained in
					while ((curr_cell < xdays.size()-1) && (xdays[curr_cell+1] <= from_day)) curr_cell++;
					if (curr_cell >= xdays.size()-1) break; // we finished this part
					int left = to_day - from_day;
					if (to_day < xdays[0]) left = 0;
					//MLOG("from %d , to %d , curr %d , left %d\n", from_day, to_day, curr_cell, left);
					while (left > 0) {
						if (from_day + left < xdays[curr_cell+1]) {
							// all contained in current
							//MLOG("add1: before : from %d , to %d , curr %d , left %d\n", from_day, to_day, curr_cell, left);
							hmap[s][curr_cell] += (float)(left+1);
							last_day_counted = to_day;
							left = 0;
							//MLOG("add1: after : from %d , to %d , curr %d , left %d\n", from_day, to_day, curr_cell, left);
						}
						else {
							// count to the end and advance to next 
							//MLOG("add2: before : from %d , to %d , curr %d , left %d add %d \n", from_day, to_day, curr_cell, left, xdays[curr_cell+1] - from_day);
							hmap[s][curr_cell] += (float)(xdays[curr_cell+1] - from_day);
							curr_cell++;
							if (curr_cell >= xdays.size()-1) break; // we finished this part
							from_day = xdays[curr_cell];
							left = to_day - from_day;
							last_day_counted = from_day - 1;
							//MLOG("add2: after : from %d , to %d , curr %d , left %d\n", from_day, to_day, curr_cell, left);
						}
					}
				}
			}
		}

		// hmap[s] needs to be normalized to coverage
		for (int i=0; i<hmap[s].size(); i++)
			hmap[s][i] /= (float)(xdays[i+1] - xdays[i]);
	}

	// Shrinkage
	// We leave at most 2 years of no drugs at start, and 1 at the end.
	// Also we get rid of drugs with 0 usage
	first_nonz = first_nonz - 20000;
	last_nonz = last_nonz + 10000;
	int start_cell = -1, end_cell = -1;
	for (int i=0; i<xdates.size(); i++) {
		if (first_nonz >= xdates[i]) start_cell++;
		if (last_nonz >= xdates[i]) end_cell++;
	}
	if (start_cell < 0) start_cell = 0;
	_xdates.clear();
	//MLOG("min_date %d max_date %d first_nonz %d last_nonz %d start_cell %d end_cell %d xdates.size %d\n", min_date, max_date, first_nonz, last_nonz, start_cell, end_cell, xdates.size());
	for (int i=start_cell; i<=end_cell; i++) _xdates.push_back(xdates[i]);
	_hmap.clear();
	_sets_names.clear();
	for (int s=0; s<sets_sum.size(); s++) {
		if (sets_sum[s] > 0) {
			_sets_names.push_back(drugs[s]);
			_hmap.push_back({});
			for (int i=start_cell; i<=end_cell; i++) _hmap.back().push_back(hmap[s][i]);
		}
	}



	// debug

#if 0
	MLOG("xdates: ");
	for (int i=0; i<xdates.size(); i++) MLOG("[%d] %d:%d ", i, xdates[i], xdays[i]);
	MLOG("\n");
	for (int s=0; s<sets.size(); s++) {
		MLOG("%s : ", sets[s][0].c_str());
		for (int i=0; i<hmap[s].size(); i++) MLOG("[%d] %f ", i, hmap[s][i]);
		MLOG("\n");
	}
#endif

}

//----------------------------------------------------------------------------------------
string MedPatientPlotlyDate::date_to_string(int date)
{
	int y = date/10000;
	int m = (date % 10000)/100;
	int d = date % 100;

	string s = "'" + to_string(y) + "-" + to_string(m) + "-" + to_string(d) + "'";
	return s;
}

//----------------------------------------------------------------------------------------
string MedPatientPlotlyDate::time_to_string(int time, int time_unit)
{
	if (time_unit < 0) {
		if (time < 21000000)
			time_unit = MedTime::Date;
		else
			time_unit = MedTime::Minutes;
	}

	if (time_unit == MedTime::Date)
		return date_to_string(time);

	// left is the minutes case , we want to get to "YYYY-MM-DD hh:mm:ss" ss is always 0
	string s = med_time_converter.convert_times_S(MedTime::Minutes, MedTime::DateTimeString, time);

	string out_s = "'" + s.substr(0, 4) + "-" + s.substr(4, 2) + "-" + s.substr(6, 2) + " " + s.substr(9, 2) + ":" + s.substr(12, 2) + ":00'";
	//MLOG("stime is %s out_s %s\n", s.c_str(), out_s.c_str());

	return out_s;
}

//----------------------------------------------------------------------------------------
void MedPatientPlotlyDate::get_usv_min_max(UniversalSigVec &usv, float &vmin, float &vmax)
{
	vmin = (float)1e10;
	vmax = -vmin;
	for (int chan=0; chan<usv.n_val_channels(); chan++)
		for (int i=0; i<usv.len; i++) {
			float v = usv.Val(i, chan);
			if (v < vmin) vmin = v;
			if (v > vmax) vmax = v;
		}
}

//----------------------------------------------------------------------------------------
void MedPatientPlotlyDate::add_xy_js(string &shtml, UniversalSigVec &usv, int time_chan, int chan, int null_zeros_flag, string prefix)
{
	// dates
	shtml += prefix+"x: [";
	for (int i=0; i<usv.len; i++) {
		shtml += time_to_string(usv.Time(i, time_chan));
		if (i < usv.len - 1)	shtml += ",";
	}
	shtml += "],\n";

	// vals
	shtml += prefix + "y: [";
	for (int i=0; i<usv.len; i++) {
		float v = usv.Val(i, chan);
		if (null_zeros_flag == 0 || v > 0)
			shtml += to_string(v);
		else
			shtml += "null";
		if (i < usv.len - 1)	shtml += ",";
	}
	shtml += "],\n";

}

//----------------------------------------------------------------------------------------
void MedPatientPlotlyDate::add_xy_js(string &shtml, vector<int> &dates, vector<float> &vals, int null_zeros_flag, string prefix)
{
	// dates
	shtml += prefix+"x: [";
	for (int i=0; i<dates.size(); i++) {
		shtml += time_to_string(dates[i]);
		if (i < dates.size() - 1)	shtml += ",";
	}
	shtml += "],\n";

	// vals
	shtml += prefix + "y: [";
	for (int i=0; i<vals.size(); i++) {
		float v = vals[i];
		if (null_zeros_flag == 0 || v > 0)
			shtml += to_string(v);
		else
			shtml += "null";
		if (i < vals.size() - 1)	shtml += ",";
	}
	shtml += "],\n";

}

//----------------------------------------------------------------------------------------
void MedPatientPlotlyDate::add_dataset_js(string &shtml, UniversalSigVec &usv, int time_chan, int chan, int null_zeros_flag, string prefix, string sname, int yaxis, string sig)
{

	shtml += prefix + "var " + sname + " = {\n";
	add_xy_js(shtml, usv, time_chan, chan, null_zeros_flag, prefix + "\t");

	// types/general defs
	shtml += prefix + "\ttype: 'scatter',\n";
	shtml += prefix + "\tmode: 'lines+markers',\n";
	shtml += prefix + "\tline: {shape: 'spline', width: 2, smoothing: 0.75},\n";
	shtml += prefix + "\tyaxis: 'y" + to_string(yaxis) + "',\n";
	if (null_zeros_flag) shtml += prefix + "\tconnectgaps: true,\n";
	shtml += prefix + "\tname: '" + sig +"'\n";

	// close it
	shtml += prefix + "};\n";
}


//----------------------------------------------------------------------------------------
void MedPatientPlotlyDate::add_bg_dataset_js(string &shtml, vector<int> &dates, vector<float> &vals, int null_zeros_flag, string color, string prefix, string sname, int yaxis, string name)
{

	shtml += prefix + "var " + sname + " = {\n";
	add_xy_js(shtml, dates, vals, null_zeros_flag, prefix + "\t");
	/*
	name: 'Statins',
	yaxis: 'y4',
	type: 'scatter',
	mode: 'none',
	fill: 'tozeroy',
	fillcolor: 'rgba(162, 217, 206,0.333)',
	line: {shape: 'hv'}
	*/
	// types/general defs
	shtml += prefix + "\ttype: 'scatter',\n";
	shtml += prefix + "\tmode: 'none',\n";
	shtml += prefix + "\thoverinfo: 'none',\n";
	shtml += prefix + "\tfill: 'tozeroy',\n";
	shtml += prefix + "\tfillcolor: '"+color+"',\n";
	shtml += prefix + "\tline: {shape: 'hv'},\n";
	shtml += prefix + "\tyaxis: 'y" + to_string(yaxis) + "',\n";
	if (null_zeros_flag) shtml += prefix + "\tconnectgaps: true,\n";
	shtml += prefix + "\tname: '" + name +"'\n";

	// close it
	shtml += prefix + "};\n";
}


//----------------------------------------------------------------------------------------
int MedPatientPlotlyDate::add_panel_chart(string &shtml, LocalViewsParams &lvp, PidRec &rec, const PanelInfo &pi, const vector<ChartTimeSign> &times)
{
	int pid = rec.pid;
	int def_time_chan = 0;
	int pwidth = (pi.width < 0) ? params.width_default : pi.width;
	int pheight = (pi.height < 0) ? params.height_default : pi.height;
	int block_mode = (pi.block_mode < 0) ? params.block_mode_default : pi.block_mode;
	int null_zeros = (pi.null_zeros < 0) ? params.null_zeros_default : pi.null_zeros;
	int log_scale = (pi.log_scale < 0) ? params.log_scale_default : pi.log_scale;

	vector<string> titles = pi.sigs;

	// div_name
	string div_name = "div";
	for (auto &s : pi.sigs) div_name += "_" + s;
	div_name += to_string(rand_N(10000));


	// computing datasets

	UniversalSigVec usv;
	int cnt = 0;
	int tot_len = 0;
	string shtml_sets;
	int n_yaxis = (int)pi.sigs.size();
	vector<float> vmin(pi.sigs.size()), vmax(pi.sigs.size());
	for (int i=0; i<pi.sigs.size(); i++) {
		rec.uget(pi.sigs[i], usv);
		int time_chan = def_time_chan;
		if (params.sig_params.find(pi.sigs[i]) != params.sig_params.end())
			time_chan = params.sig_params[pi.sigs[i]].time_chan;

		tot_len += usv.len;
		for (int chan=0; chan<usv.n_val_channels(); chan++) {
			add_dataset_js(shtml_sets, usv, time_chan, chan, null_zeros, "\t\t", "set" + to_string((++cnt)), i+1, pi.sigs[i]);
		}
		get_usv_min_max(usv, vmin[i], vmax[i]);
	}

	if (tot_len == 0) return 0;

	if (pi.drugs.size() > 0) {
		for (int i=0; i<pi.drugs.size(); i++) {
			vector<string> dname;
			vector<int> xdates;
			vector<vector<float>> hmap;
			get_drugs_heatmap(rec, xdates, dname, hmap, { pi.drugs[i] });
			if (xdates.size()>0 && hmap.size()>0) {
				string color = PlotlyColorDefaults::bg_opaque_colors[i % PlotlyColorDefaults::bg_opaque_colors.size()];
				if (i < pi.drug_colors.size()) color = pi.drug_colors[i];
				add_bg_dataset_js(shtml_sets, xdates, hmap[0], 0, color, "\t\t", "set" + to_string(++cnt), n_yaxis+1, pi.drugs[i]);
				n_yaxis++;
				vmin.push_back(0);
				vmax.push_back(0);
				titles.push_back(pi.drugs[i]);
			}
		}

	}


	// set height , width, and block_mode
	shtml += "\t<div id=\"" + div_name + "\" style=\"width:" + to_string(pwidth) + "px;height:" + to_string(pheight) + "px;";
	if (block_mode) shtml += "display: inline-block;";
	shtml+= "\"></div>\n";
	shtml += "\t<script>\n";

	shtml += shtml_sets;

	// prep layout

	shtml += "\t\tvar layout = {\n";
	shtml += "\t\t\ttitle: '" + pi.title + "',\n";
	//float psize = (float)1.0 - n_yaxis*(float)0.02;
	float psize = (float)0.98;
	// deal with multiple yaxis
	// deal with multiple yaxis
	for (int i=0; i<n_yaxis; i++) {
		if (i == 0)
			shtml += "\t\t\tyaxis" + to_string(i+1) +": {title: '" + titles[i] + "', showline: false";

		if (i > 0) {
			shtml += "\t\t\tyaxis" + to_string(i+1) +": {showline: false";
			//shtml += ", overlaying: 'y', side: 'right', position: " + to_string(psize+0.02*i) + ", tick: '', showticklabels: false";
			shtml += ", overlaying: 'y', side: 'right', position: " + to_string(psize) + ", tick: '', showticklabels: false";
		}
		if (log_scale && vmin[i] > 0 && usv.n_val_channels() < 2) shtml += ",type: 'log', autorange: true";
		shtml += "},\n";
	}
	// xaxis setup
	string from_t, to_t;
	if (lvp.from_date > 0 && lvp.to_date >= lvp.from_date) {
		from_t = date_to_string(lvp.from_date);
		to_t = date_to_string(lvp.to_date);
		if (params.rep_time_unit == MedTime::Minutes) {
			from_t.pop_back();
			from_t += " 00:00'";
			to_t.pop_back();
			to_t += " 23:59'";			
		}
	}

	shtml += "\t\t\txaxis: { omain: [0," + to_string(psize) +"], ";
	if (from_t != "") shtml += "range: [" + from_t + "," + to_t + "], ";
	if (params.rep_time_unit == MedTime::Date) 
		shtml += "hoverformat: '%Y/%m/%d'},\n";
	else
		shtml += "hoverformat: '%Y/%m/%d %H:%M'},\n";
	if (times.size() > 0) {
		shtml += "\t\t\tshapes: [";
		for (auto &t : times) {
			string ts = date_to_string(t.time);
			if (t.name == "Death") ts = time_to_string(t.time);
			shtml += "{type: 'line', x0: " + ts + ", y0: " + to_string(vmin[0]);
			shtml += ", x1: " + ts + ", y1: " + to_string(vmax[0]);
			string color = t.color; //"'black'";
			shtml += ", line: { color: " + color + "} },";
		}
		shtml += "]\n";
	}
	shtml += "\t\t};\n";

	// prep data variable
	shtml += "\t\tvar data = [";
	for (int i=0; i<cnt; i++) {
		string set_name = "set" + to_string(i+1);
		shtml += set_name;
		if (i < cnt-1) shtml += ",";
	}
	shtml += "];\n";

	// actual plot
	shtml += "\t\tPlotly.plot('" + div_name + "', data, layout);\n";

	shtml += "\t</script>\n";

	return 0;
}

//-------------------------------------------------------------------------------------------------------------------------------
void MedPatientPlotlyDate::add_thin_rc_chart(string &shtml, PidRec &rec, const vector<ChartTimeSign> &times)
{
	int pid = rec.pid;
	int time_chan = 0;
	int pwidth = 1200; //vm["pwidth"].as<int>();
	int pheight = 600; //vm["pheight"].as<int>();
	int block_mode = 0; //vm["block_mode"].as<int>();
	int rc_acc_days = 60;
	unordered_map<string, double> code2wgt ={ { "0", 0.1 },{ "1", 0.1 },{ "2", 0.1 },{ "3", 0.1 },{ "4", 0.1 },{ "5", 0.5 },{ "6", 0.1 },{ "7", 1.0 },{ "8", 0.1 },{ "9", 0.1 },
	{ "A", 2.0 },{ "B", 5.0 },{ "C", 2.5 },{ "D", 1.5 },{ "E", 1.5 },{ "F", 2.0 },{ "G", 3.0 },{ "H", 3.0 },{ "J", 2.0 },{ "K", 1.0 },
	{ "L", 2.5 },{ "M", 0.8 },{ "N", 1.0 },{ "P", 1.0 },{ "Q", 2.0 },{ "R", 0.9 },{ "S", 1.0 },{ "T", 1.0 },{ "U", 2.0 },{ "Z", 0.1 }
	};
	double p_decay = 0.5, q_decay = 0.2;
	double alpha;

	alpha = -log(q_decay/p_decay)/(double)rc_acc_days;

	if (rec.my_base_rep->sigs.sid("RC") <= 0)
		return; // NO RC signal, nothing to do, maybe not THIN??

	// plan:
	// calculate some accumulator in a time window and add a point of (x=time,y=accumulator,text=drugs in day x)

	vector<string> xdates;
	vector<float> yvals;
	vector<int> days;
	vector<float> days_cnt;
	vector<string> hovertext;

	UniversalSigVec usv;
	rec.uget("RC", usv);

	if (usv.len == 0) return; // nothing to do - 0 RC records.

	int curr_date = 0;
	int curr_days = 0;
	int section_id = rec.my_base_rep->dict.section_id("RC");
	for (int i=0; i<usv.len; i++) {
		int i_date = usv.Time(i, 0);
		int i_val = (int)usv.Val(i, 0);
		int d = med_time_converter.convert_date(MedTime::Days, i_date);
		if (d != curr_date) {
			days.push_back(d); days_cnt.push_back(0);
			xdates.push_back(date_to_string(i_date));
			float y = 0;
			int j = (int)days_cnt.size() - 1;
			while (j >= 0) {
				int ddays = d - days[j];
				if (ddays > rc_acc_days) break;
				double factor = p_decay * exp(-alpha*(double)ddays);
				y += (float)factor*days_cnt[j];
				j--;
			}
			yvals.push_back(y);
			hovertext.push_back("");
		}
		curr_date = d;
		//days_cnt.back()++;
		//yvals.back()++;

		// recover curr text
		string curr_text = "";
		if (rec.my_base_rep->dict.dict(section_id)->Id2Names.find(i_val) != rec.my_base_rep->dict.dict(section_id)->Id2Names.end())
			for (int j = 0; j < rec.my_base_rep->dict.dict(section_id)->Id2Names[i_val].size(); j++) {
				string sname = rec.my_base_rep->dict.dict(section_id)->Id2Names[i_val][j];
				string scode = sname.substr(0, 1);

				curr_text += "|" + sname;

				if (j == 0) { // assumes the first name is the READCODE !!!!!!!!
					float wgt = (float)code2wgt[scode];
					days_cnt.back() += wgt;
					yvals.back() += wgt;
					//MLOG("date %d sname= %s :: scode = %s :: weight = %f :: acc %f\n", i_date, sname.c_str(), scode.c_str(), code2wgt[scode], days_cnt.back());
				}
			}
		replace(curr_text.begin(), curr_text.end(), '\"', '@');
		replace(curr_text.begin(), curr_text.end(), '\'', '@');

		if (hovertext.back() != "") hovertext.back() += "<br>";
		hovertext.back() += curr_text;
	}

	// prep x , y, text arrays
	string ax, ay, atext;
	float ymax = 0, ymin = 9999999;
	for (int i=0; i<xdates.size(); i++) {
		ax += xdates[i];
		ay += to_string(yvals[i]);
		atext += "\"" + hovertext[i] + "\"";
		if (yvals[i] > ymax) ymax = yvals[i];
		if (yvals[i] < ymin) ymin = yvals[i];
		//atext += "\"TTT" + to_string(i) + "\"";
		if (i < xdates.size()-1) {
			ax += ",";
			ay += ",";
			atext += ",";
		}
		//MLOG("i %d :: days %d,%d :: xdate %s :: yval %d :: %s\n", i, days[i], days_cnt[i], xdates[i].c_str(), yvals[i], hovertext[i].c_str());
	}

	// write RC div
	// div_name
	string div_name = "div_RC_";
	div_name += to_string(rand_N(10000));

	//<div id="chart" style="width:1200px;height:500px;"></div>
	shtml += "\t<div id=\"" + div_name + "\" style=\"width:" + to_string(pwidth) + "px;height:" + to_string(pheight) + "px;";
	if (block_mode) shtml += "display: inline-block;";
	shtml+= "\"></div>\n";
	shtml += "\t<script>\n";

	shtml += "\t\tvar set1 = {\n";
	shtml += "\t\t\tx: [" + ax + "],\n";
	shtml += "\t\t\ty: [" + ay + "],\n";
	shtml += "\t\t\ttext: [" + atext + "],\n";
	shtml += "\t\t\tyaxis: 'y1',\n";
	shtml += "\t\t\ttype: 'scatter',\n";
	shtml += "\t\t\tname: 'ReadCodes',\n";
	shtml += "\t\t\tmode: 'lines+markers',\n";
	shtml += "\t\t\tline: {shape: 'spline', width: 2, smoothing: 0.75},\n";
	shtml += "\t\t\thoverinfo: 'x+text'\n";
	shtml += "\t\t};\n";

	shtml += "\t\tvar layout ={\n";
	shtml += "\t\t\ttitle: 'ReadCodes',\n";
	shtml += "\t\t\tyaxis: {autorange: true},\n";
	if (times.size() > 0) {
		shtml += "\t\t\tshapes: [";
		for (auto &t : times) {
			shtml += "{type: 'line', x0: " + date_to_string(t.time) + ", y0: " + to_string(ymin);
			shtml += ", x1: " + date_to_string(t.time) + ", y1: " + to_string(ymax);
			string color = t.color; //"'black'";
			shtml += ", line: { color: " + color + "} },";
		}
		shtml += "]\n";

	}

	//if (time > 0) shtml += "\t\t\tshapes: [{type: 'line', x0: " + date_to_string(time) + ", y0: " + to_string(ymin) + ", x1: " + date_to_string(time) + ", y1: " + to_string(ymax) + "}]\n";
	shtml += "\t\t};\n";

	shtml += "\t\tPlotly.plot('" + div_name+"', [set1], layout);\n";

	shtml += "\t</script>\n";

}

//-------------------------------------------------------------------------------------------------------------------------------
int MedPatientPlotlyDate::add_drugs_heatmap(string &shtml, PidRec &rec)
{
	vector<int> xdates;
	vector<vector<float>> hmap;

	vector<string> sets_names;

	//MLOG("##>> heatmap for: "); for (auto d : params.dhm_params.drugs) MLOG("%s ", d.c_str()); MLOG("\n");

	get_drugs_heatmap(rec, xdates, sets_names, hmap, params.dhm_params.drugs);

	//MLOG("##>> got only : "); for (auto d : sets_names) MLOG("%s ", d.c_str()); MLOG("\n");

	int pid = rec.pid;
	int time_chan = 0;
	int dwidth = params.dhm_params.width;
	int dheight = params.dhm_params.height + 30*(int)sets_names.size();
	int block_mode = 0;

	string div_name = "div_drug_heatmap" + to_string(rand_N(10000));

	shtml += "\t<div id=\"" + div_name + "\" style=\"width:" + to_string(dwidth) + "px;height:" + to_string(dheight) + "px;";
	if (block_mode) shtml += "display: inline-block;";
	shtml+= "\"></div>\n";
	shtml += "\t<script>\n";

	// xValues is xdates
	shtml += "\t\tvar xValues = [";
	for (int i=0; i<xdates.size(); i++) {
		//MLOG("====> xValues %d : %d\n", i, xdates[i]);
		shtml += date_to_string(xdates[i]); if (i<xdates.size()-1) shtml += ",";
	}
	shtml += "];\n";

	// yValues is the sets names
	shtml += "\t\tvar yValues = [";
	for (int i=0; i<sets_names.size(); i++) {
		shtml += "'" + sets_names[i] + "'"; if (i<sets_names.size()-1) shtml += ",";
	}
	shtml += "];\n";

	// zValues is tje actual heatmap
	shtml += "\t\tvar zValues = [";
	for (int i=0; i<sets_names.size(); i++) {
		shtml += "[";
		for (int j=0; j<xdates.size(); j++) {
			shtml += to_string(hmap[i][j]);
			if (j<xdates.size()-1) shtml += ",";
		}
		shtml += "]";
		if (i<sets_names.size()-1) shtml += ",";
	}
	shtml += "];\n";

	// color scales
	//shtml += "\t\tvar colorscaleValue = [ [0, '#3D9970'], [1, '#001f3f']];\n";
	shtml += "\t\tvar colorscaleValue = [ [0, '#001f3f'], [1, '#2ecc71']];\n";

	// data
	shtml += "\t\tvar data = [{ x: xValues,	y: yValues,	z: zValues,	type: 'heatmap', colorscale: colorscaleValue, showscale: false}];\n";

	// layout
	shtml += "\t\tvar layout = { title: 'Drugs HeatMap', margin: {l: 200}};\n";

	// actual plot
	shtml += "\t\tPlotly.plot('" + div_name + "', data, layout);\n";

	shtml += "\t</script>\n";

	return 0;
}

//----------------------------------------------------------------------------------------------------------------------------------------
int MedPatientPlotlyDate::get_rec_html(string &shtml, LocalViewsParams &lvp, PidRec &rec, const string &mode, const vector<ChartTimeSign> &sign_times, const vector<string> &view)
{
	shtml = "";

	add_html_header(shtml, mode);

	vector<ChartTimeSign> local_sign_times = sign_times;

	shtml += "<body>\n";

	add_basic_demographics(shtml, rec, local_sign_times);

	MLOG("After demographics\n");
	for (auto v : view) {
		MLOG("Working on view %s\n", v.c_str());
		if (v == "demographic") {
			add_basic_demographics(shtml, rec, local_sign_times);
		}
		else if (params.panels.find(v) != params.panels.end()) {
			// add a panel
			add_panel_chart(shtml, lvp, rec, params.panels[v], local_sign_times);
		}
		else if (v == "RC") {
			add_thin_rc_chart(shtml, rec, local_sign_times);
		}
		else if (rec.my_base_rep->sigs.sid(v) > 0) {
			// add a signal (as a simple panel)
			int null_zeros = -1;
			int log_scale = -1;
			if (params.sig_params.find(v) != params.sig_params.end()) {
				null_zeros = params.sig_params[v].null_zeros;
				log_scale = params.sig_params[v].log_scale;
			}
			PanelInfo pi(v, v, { v }, {}, {}, null_zeros, log_scale);
			params.panels[v] = pi;
			add_panel_chart(shtml, lvp, rec, params.panels[v], local_sign_times);
		}
		else if (v == "drugs_heatmap") {
			add_drugs_heatmap(shtml, rec);
		}

	}
	
	// add_RCs_to_js(rep, vm, shtml);
	MLOG("Finished preparing\n");
	shtml += "</body>\n";
	shtml += "</html>\n";

	return 0;
}
