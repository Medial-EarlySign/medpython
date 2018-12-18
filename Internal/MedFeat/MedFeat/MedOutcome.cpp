//
// MedOutcome
//

#define _SCL_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

#include "MedOutcome.h"
#include "Logger/Logger/Logger.h"
#include "MedUtils/MedUtils/MedUtils.h"

#include <fstream>
#include <algorithm>
#include <thread>

#define LOCAL_SECTION LOG_MEDSTAT
#define LOCAL_LEVEL	LOG_DEF_LEVEL
extern MedLogger global_logger;


// OLD CODE NOT NEEDED

#if 0

//-------------------------------------------------------------------------------------------
int MedOutcome::read_from_file(const string &fname)
{
	ifstream inf(fname);

	MLOG("MedOutcome: reading %s\n", fname.c_str());
	if (!inf) {
		MERR("MedOutcome: can't open file %s for read\n",fname.c_str());
		return -1;
	}

	name = "NONE";
	desc = "";
	type = "single";
	n_categ = -1;
	out.clear();
	ignore_pids.clear();
	pid2use.clear();

	string curr_line;
	while (getline(inf,curr_line)) {
		//MLOG("--> %s\n",curr_line.c_str());
		if ((curr_line.size() > 1) && (curr_line[0] != '#')) {

			if (curr_line[curr_line.size()-1] == '\r')
				curr_line.erase(curr_line.size()-1) ;

			vector<string> fields;
			split(fields, curr_line, boost::is_any_of("\t"));

			if (fields.size() >= 2) {

				if (fields[0] == "NAME") name = fields[1];
				if (fields[0] == "DESC") desc += fields[1]+"\n";
				if (fields[0] == "TYPE") type = fields[1];
				if (fields[0] == "NCATEG") n_categ = stoi(fields[1]);
				if (fields[0] == "EVENT" && fields.size() >= 4) {
					//MLOG("-->### %s %s (%d)\n",fields[0].c_str(),fields[1].c_str(),out.size());
					OutcomeRecord orec;
					orec.pid = stoi(fields[1]);
					orec.time = stoll(fields[2]);
					if (orec.time < 30000000)
						orec.date = (int)orec.time;
					else
						orec.date = 0;
					orec.categ = stof(fields[3]);
					if (fields.size()>4)
						orec.length = stoi(fields[4]);
					else
						orec.length = 0;
					if (fields.size()>5) orec.raw_time = stoll(fields[5]); else orec.raw_time = orec.time;
					if (orec.raw_time < 30000000)
						orec.raw_date = (int)orec.raw_time;
					else
						orec.raw_date = 99999999;
					out.push_back(orec);
					//MLOG("(%d) %d 1\n", pid2use.size(), orec.pid);
					if (pid2use.find(orec.pid) == pid2use.end())
						pid2use[orec.pid] = 1;
				}
				if (fields[0] == "IGNORE" && fields.size() == 2)  {
					int pid = stoi(fields[1]);
					ignore_pids.push_back(pid);
					pid2use[pid] = -1;
				}

			}
		}
	}

	inf.close();
	return 0;
}

//-------------------------------------------------------------------------------------------
int MedOutcome::write_to_file(const string &fname)
{
	ofstream of;

	of.open(fname,ios::out);

	if (!of) {
		MERR("MedOutcome: can't open file %s for write\n",fname.c_str());
		return -1;
	}

	of << "NAME\t" << name.c_str() << "\n";
	of << "DESC\t" << desc.c_str() << "\n";
	of << "TYPE\t" << type.c_str() << "\n";
	for (int i=0; i<out.size(); i++) {
		unsigned long long datetime = out[i].time;
		unsigned long long raw_datetime = out[i].raw_time;
		if (datetime == 0) datetime = out[i].date;
		if (raw_datetime == 0) raw_datetime = out[i].raw_date;
		of << "EVENT\t" << out[i].pid << "\t" << datetime << "\t" << out[i].categ << "\t" << out[i].length << "\t" << raw_datetime << "\n";
	}

	for (int i=0; i<ignore_pids.size(); i++) {
		of << "IGNORE\t" << ignore_pids[i] << "\n";
	}

	of.close();
	return 0;
}

//-------------------------------------------------------------------------------------------
void MedOutcome::subsample(float prob)
{
	vector<OutcomeRecord> sub_out;

	for (int i=0; i<out.size(); i++)
		if (rand_1() < prob)
			sub_out.push_back(out[i]);

	out = sub_out;
}

//-------------------------------------------------------------------------------------------
void OutcomeFilter::init_from_string(const string &s)
{
	vector<string> fields;
	split(fields, s, boost::is_any_of(",="));
	sigs.clear();

	for (int i=0; i<fields.size(); i++) {

		if (fields[i] == "filter_flag") filter_flag = stoi(fields[++i]);
		else if (fields[i] == "filter_just_demographics") filter_just_demographics = stoi(fields[++i]);
		else if (fields[i] == "min_date") min_date = stoi(fields[++i]);
		else if (fields[i] == "max_date") max_date = stoi(fields[++i]);
		else if (fields[i] == "train_to_take") train_to_take = stoi(fields[++i]);
		else if (fields[i] == "gender") gender = stoi(fields[++i]);
		else if (fields[i] == "min_age") min_age = stof(fields[++i]);
		else if (fields[i] == "max_age") max_age = stof(fields[++i]);
		else if (fields[i] == "sigs") sigs.push_back(fields[++i]);
		else if (fields[i] == "sigs_win_min") sigs_win_min = stoi(fields[++i]);
		else if (fields[i] == "sigs_win_max") sigs_win_max = stoi(fields[++i]);
		else if (fields[i] == "min_sigs_time_points_in_win") min_sigs_time_points_in_win = stoi(fields[++i]);
		else if (fields[i] == "stick_to_sigs") stick_to_sigs = stoi(fields[++i]);
		else if (fields[i] == "stick_to_sigs") stick_to_sigs = stoi(fields[++i]);
		else if (fields[i] == "n_points_to_take_in_win") n_points_to_take_in_win = stoi(fields[++i]);
		else if (fields[i] == "clock_sigs") clock_sigs.push_back(fields[++i]);
		else if (fields[i] == "match_flag") match_flag = stoi(fields[++i]);
		else if (fields[i] == "match_age_flag") match_age_flag = stoi(fields[++i]);
		else if (fields[i] == "age_bin_width") age_bin_width = stof(fields[++i]);
		else if (fields[i] == "match_gender_flag") match_gender_flag = stoi(fields[++i]);
		else if (fields[i] == "match_const_flag") match_const_flag = stoi(fields[++i]);
		else if (fields[i] == "match_do_shuffle_flag") match_do_shuffle_flag = stoi(fields[++i]);
		else if (fields[i] == "match_dates_flag") match_dates_flag = stoi(fields[++i]);
		else if (fields[i] == "dates_bin_width") dates_bin_width = stof(fields[++i]);
		else if (fields[i] == "match_location_flag") match_location_flag = stoi(fields[++i]);
		else if (fields[i] == "match_signals_flag") match_signals_flag = stoi(fields[++i]);
		else if (fields[i] == "match_sig") {sigs_to_match.push_back(fields[++i]); sig_bin_width.push_back(stof(fields[++i]));}
		else if (fields[i] == "match_min_size") match_min_size = stoi(fields[++i]);
		else if (fields[i] == "match_max_ratio") match_max_ratio = stoi(fields[++i]);
		else if (fields[i] == "match_event_case_price_ratio") match_event_case_price_ratio = stof(fields[++i]);		
		else if (fields[i] == "match_max_samples") match_max_samples = stoi(fields[++i]);
		else MWARN("OutcomeFilter: unknown token [%s] \n", fields[i].c_str());
	}
}



//====================================================================================================
// MedFeatList
//====================================================================================================
int FeatParams::init(const string &in_name)
{
	full_name = in_name;

	vector<string> fields;
	split(fields, in_name, boost::is_any_of("$"));

	sig_name = fields[0];
	if (fields.size() > 0)	sig_name = fields[0]; else sig_name = "";
	if (fields.size() > 1)	feat_name = fields[1]; else feat_name = "";

	// Age - special case
	if (sig_name == "Age" || sig_name == "AGE") {
		if (feat_name == "Param") {
			f_type = F_Age_Param;
		}
		else {
			f_type = F_Age_Calculated;
			return 0;
		}
	}

	// Constant signals (SVals) that are related only to the id (demogrpahics, etc)
	if (feat_name == "" || feat_name == "Demographic") f_type = F_Const;
	if (feat_name.compare(0, 4, "Last")==0) f_type = F_Last;
	if (feat_name.compare(0, 5, "TLast")==0) f_type = F_TLast;
	if (feat_name.compare(0, 3, "Avg")==0) f_type = F_Avg;
	if (feat_name.compare(0, 3, "Max")==0) f_type = F_Max;
	if (feat_name.compare(0, 5, "Delta")==0) f_type = F_Delta;
	if (feat_name.compare(0, 5, "Slope")==0) f_type = F_Slope;
	if (feat_name.compare(0, 4, "NumV")==0) f_type = F_NumV;

	// sid and sig_type will be initialized someplace else, as we need a specific repository for that
	if (f_type == F_Last) { ind = stoi(feat_name.substr(4)); }
	if (f_type == F_TLast) { ind = stoi(feat_name.substr(5)); }

	vector<FeatureType> from_to_params ={ F_Avg, F_Max, F_Delta, F_Slope, F_NumV };
	if (find(from_to_params.begin(), from_to_params.end(), f_type) != from_to_params.end()) {
		vector<string> fields2;
		split(fields2, feat_name, boost::is_any_of("_"));

		if (fields2.size() > 1) {
			from_time = stoll(fields2[1]);
			from_day = (int)from_time;
		}

		if (fields2.size() > 2) {
			to_time = stoll(fields2[2]);
			to_day = (int)to_time;
		}

	}

	return 0;
}

int FeatParams::init_sid(MedRepository &rep)
{
	if (sig_name != "") {
		sid = rep.sigs.sid(sig_name);
		sig_type = rep.sigs.Sid2Info[sid].type;
	}

	return 0;
}

void MedFeatList::push_feature(string feat) {
	MLOG("%d=%s, ", features.size(), feat.c_str());
	features.push_back(feat);
}
//................................................................................................
int MedFeatList::read_from_file(const string &fname)
{
	clear();

	ifstream inf(fname);

	if (!inf) {
		MERR("MedFeatList: can't open file %s for read\n",fname.c_str());
		return -1;
	}
	MLOG("reading features from: %s\n", fname.c_str());
	features.clear();
	map<string,int> signals_in;
	string curr_line;
	while (getline(inf,curr_line)) {
		if ((curr_line.size() > 1) && (curr_line[0] != '#')) {
			if (curr_line[curr_line.size()-1] == '\r')
				curr_line.erase(curr_line.size()-1) ;

			vector<string> fields;
			boost::split(fields, curr_line, boost::is_any_of(" \t"), boost::algorithm::token_compress_on);
			if (fields.size() == 1)
				push_feature(fields[0]); // regular <sig>$<feature type> kind of feature

			if (fields.size() >= 2) {
				if (fields[0] == "DESC") 
					desc = fields[1];
				else {
					// we do a cartesian product of signals against feature_types
					vector<string> sigs;
					vector<string> feat_types;
					boost::split(sigs, fields[0], boost::is_any_of(","));
					boost::split(feat_types, fields[1], boost::is_any_of(","));
					assert(sigs.size() > 0);
					assert(feat_types.size() > 0);
					for (int i=0; i<sigs.size(); i++) {
						// maintaining signals vectos
						string my_sig = "";
						if (feat_types[0] == "Special") {
							// NOTE: when adding a special signal, you must update this code = MedFeatList creation AND in collect_features
							if (sigs[i].compare(0, 5, "Hosp_") == 0) my_sig = "Hospitalization";
							else if (sigs[i].compare(0, 4, "Dep_") == 0) my_sig = "?";
							else if (sigs[i].compare(0, 5, "Drug_") == 0) my_sig = "Drug";
							else if (sigs[i].compare(0, 4, "Fit_") == 0) my_sig = "?";
							else if (sigs[i].compare(0, 7, "Anemia_") == 0) my_sig = "?";
							else if (sigs[i].compare(0, 17, "SMOKING_ENRICHED_") == 0) my_sig = "SMOKING_ENRICHED";
							else if (sigs[i].compare(0, 11, "TimeToEvent") == 0) my_sig = "##TimeToEvent";
							else {
								MERR("Unknown special signal: %s", sigs[i].c_str());
								throw exception();
							}
						}
						else {
							if (sigs[i].compare(0, 2, "BP") == 0) my_sig = "BP";
							else if (sigs[i].compare(0, 7, "ALCOHOL") == 0) my_sig = "ALCOHOL";
							else my_sig = sigs[i];
						}
						if (my_sig.compare(0,2,"##") != 0 && signals_in.find(my_sig) == signals_in.end()) {
							signals.push_back(my_sig);
							signals_in[my_sig] = 1;
						}
						MLOG("feature %s: ", sigs[i].c_str());
						// cartesian product with feature types
						for (int j = 0; j < feat_types.size(); j++) {
							push_feature(sigs[i] + "$" + feat_types[j]);							
						}
						MLOG("\n");
					}
				}
			}

		}
	}
	feats.resize(features.size());
	for (int i=0; i<features.size(); i++) {
		//MLOG("i %d - %s\n", i, features[i].c_str());
		feats[i].init(features[i]);
	}
	inf.close();

	return 0;
}
		
//................................................................................................
int MedFeatList::write_to_file(const string &fname)
{
	ofstream of;

	of.open(fname,ios::out);

	if (!of) {
		MERR("MedFeatList: can't open file %s for write\n",fname.c_str());
		return -1;
	}

	of << "DESC\t" << desc << "\n";
	for (int i=0; i<features.size(); i++) {
		of << features[i] << "\n";
	}

	of.close();

	return 0;
}

//====================================================================================================
//-------------------------------------------------------------------------------------------
// filtering a single outcome record according to the rules in filter
int filter_outcome_record(MedRepository &rep, OutcomeFilter &filter, OutcomeRecord &orec, vector<OutcomeRecord> &recs_out, bool verbose = false)
{
	if (verbose)
		MLOG("considering filter for pid [%d]\n", orec.pid);
	recs_out.clear();

	if (filter.filter_flag == 0) { 
		if (verbose) MLOG("filter_flag == 0, not filtering!\n");
		recs_out.push_back(orec); return 1; 
	}

	if (filter.min_date>=0 && orec.date>0 && (orec.date<filter.min_date || orec.date>filter.max_date)) return 0;

	int len;
	int gender = (int)((SVal *)rep.get(orec.pid,"GENDER",len))[0].val;
	if (len == 0) return 0;

	if (verbose) MLOG("gender = [%d]\n", gender);
	gender = (1 << (gender-1));	
	if ((filter.gender&gender) == 0) {
		if (verbose) MLOG("filtered by gender\n");
		return 0;
	}

	int train = (int)((SVal *)rep.get(orec.pid,"TRAIN",len))[0].val;
	train = (1 << (train-1));
	if (len > 0 && (filter.train_to_take&train)==0)
		return 0;

	int byear = (int)((SVal *)rep.get(orec.pid,"BYEAR",len))[0].val;
	if (orec.date > 0) {
		float age = get_age(orec.date,byear);
		if (verbose) MLOG("age = [%f]\n", age);
		if (age < filter.min_age || age > filter.max_age) {
			if (verbose) MLOG("filtered by age\n");
			return 0;
		}
	}

	if (filter.history_time_before > 0) {
		vector<int> dates;
		rep.get_dates_with_signal(orec.pid,filter.clock_sigs,dates);
		if (dates.size()>0) {
			if ((date_to_days(orec.date) - date_to_days(dates[0])) < filter.history_time_before)
				return 0;
		}
	}

	if (filter.filter_just_demographics) {
		if (verbose) MLOG("record passed as filter_just_demographics == 1\n"); 
		recs_out.push_back(orec); return 1; 
	} else if (verbose) MLOG("filter_just_demographics == 0, continuing to filter by signals...\n");


	if (filter.sigs.size() > 0 && orec.date > 0) {

		vector<int> dates;
		rep.get_dates_with_signal(orec.pid,filter.sigs,dates);
		if (dates.size() == 0) return 0;

		if (filter.min_sigs_time_points_in_win > 0) {

			// we now have to count how many signal points there are in the relevant window
			int reference = date_to_days(orec.date);
			int r_from = reference - filter.sigs_win_max;
			int r_to = reference - filter.sigs_win_min;
			if (r_to < r_from) swap(r_to,r_from);

			int n_in = 0;
			int dates_first = -1;
			int dates_last = -1;
			int closest_id = -1, closest_dist = 10000000;
			vector<pair<int, int>> dist_i;
			for (int i=0; i<dates.size(); i++) {
				int d = date_to_days(dates[i]);
				if (d>r_from && d<=r_to) {
					n_in++;
					pair<int, int> p_i;
					p_i.first = abs(reference - d);
					p_i.second = i;
					dist_i.push_back(p_i);
					if (dates_first < 0) dates_first = i;
					dates_last = i;
					if(abs(reference - d) < closest_dist) {
						closest_dist = abs(reference-d);
						closest_id = i;
					}
				}
			}
			
			if (n_in < filter.min_sigs_time_points_in_win)
				return 0;

			// we now have to choose the actual points, up to n_points_to_take_in_win
			if (filter.stick_to_sigs) {

				sort(dist_i.begin(), dist_i.end());
				for (int i = 0; i < min(filter.n_points_to_take_in_win, (int)dist_i.size()); i++) {
					recs_out.push_back(orec);
					recs_out[i].date = dates[dist_i[i].second];
				}
				orec.date = dates[closest_id];
			} else {
				// need to choose random points
				// TBD
			}
		}

	}
	if (verbose) MLOG("[%d] records passed\n", recs_out.size());
	return 1;
}


//-------------------------------------------------------------------------------------------
// filtering a whole outcome, returning the new filtered outcome in mout_f
void get_filtered_outcome(MedRepository &rep, MedOutcome &mout, OutcomeFilter &filter, MedOutcome &mout_f)
{
	int i;

	mout_f.name = mout.name;
	mout_f.desc = mout.desc;
	mout_f.type = mout.type;
	mout_f.ignore_pids = mout.ignore_pids;
	mout_f.out.clear();
	mout_f.pid2use.clear();

	for (i=0; i<mout_f.ignore_pids.size(); i++) 
		mout_f.pid2use[mout_f.ignore_pids[i]] = -1;

	for (i=0; i<mout.out.size(); i++) {

		OutcomeRecord orec = mout.out[i];

		vector<OutcomeRecord> recs;
		if (filter_outcome_record(rep, filter, orec, recs, i < 10)) {
			//if (recs.size() > 1) MLOG("got %d recs (categ %d) (size up to now %d)\n", recs.size(), (int)orec.categ, mout_f.out.size());
			for (int j = 0; j < recs.size(); j++) {
				mout_f.out.push_back(recs[j]);
//				mout_f.out.push_back(orec);
			}
			//if (recs.size() > 1) MLOG("got %d recs (categ %d) (size up to now(2) %d)\n", recs.size(), (int)orec.categ, mout_f.out.size());
			mout_f.pid2use[orec.pid] = (int)orec.categ;
		} else {
			// ??? no real need to ignore these...we simply filter them
			//mout_f.ignore_pids.push_back(orec.pid);
			//mout_f.pid2use[orec.pid] = -1;
		}

	}
}

//-------------------------------------------------------------------------------------------
// gets a single type outcome, and transforms it to a binary outcome by randomizing a measure point (a point in time with at least one signal from sigs)
void simple_single_to_binary_outcome(MedRepository &rep, MedOutcome &mout_s, MedOutcome &mout_b, vector<string> &sigs)
{
	int i;
	vector<int> dates;

	mout_b.clear();
	mout_b.name = mout_s.name;
	mout_b.desc = mout_b.desc;
	mout_b.type = "binary";
	mout_b.ignore_pids.clear();
	mout_b.out = mout_s.out; // initializing binary outcome with all the events
	mout_b.pid2use.clear();

	OutcomeRecord orec;
	// adding all controls
	for (i=0; i<rep.index.pids.size(); i++) {
		int pid = rep.index.pids[i];

		if (mout_s.pid2use.find(pid) == mout_s.pid2use.end()) { // not in ignore list and not an event -> a control
			dates.clear();
			rep.get_dates_with_signal(pid,sigs,dates);
			
			if (dates.size() > 0) {

				int r = rand_N((int)dates.size());
				orec.pid = pid;
				orec.date = dates[r];
				orec.raw_date = dates[r];
				orec.categ = 0;
				orec.length = MAX_OUTCOME_LENGTH;
				mout_b.out.push_back(orec);
			}
					
		} 
	}

}

//-------------------------------------------------------------------------------------------
void get_outcome_record_signature(MedRepository &rep, OutcomeFilter &filter, OutcomeRecord &orec, string &signature)
{
	int len;

	signature = "";

	if (filter.match_const_flag) {
		signature += "Const : ";
	}

	if (filter.match_age_flag) {
		int byear = (int)((SVal *)rep.get(orec.pid,"BYEAR",len))[0].val;
		float age = get_age(orec.date,byear);
		int age_bin = (int)(age/filter.age_bin_width);
		signature += "Age " + to_string(age_bin*int(filter.age_bin_width)) + " : ";
	}

	if (filter.match_gender_flag) {
		int gender = (int)((SVal *)rep.get(orec.pid,"GENDER",len))[0].val;		
		signature += "Gender " + to_string(gender) + " : ";
	}

	if (filter.match_dates_flag) {
		int days = date_to_days(orec.date);
		float d_bin = (float)days/filter.dates_bin_width;
		signature += "Date " + to_string((int)d_bin) + " : ";
	}

	if (filter.match_location_flag) {
//		int branch = (int)((SVal *)rep.get(orec.pid,"BRANCH",len))[0].val;
		int branch = (int)((SVal *)rep.get(orec.pid,"DISTRICT",len))[0].val;
		signature += "Branch " + to_string(branch) + " : ";
	}

	if (filter.match_signals_flag) {
		for (int i=0; i<filter.sigs_to_match.size(); i++) {
			SDateVal *sdv = (SDateVal *)rep.get_date(orec.pid,filter.sigs_to_match[i],orec.date,"==");
			signature += filter.sigs_to_match[i] + " ";
			if (sdv != NULL) {
				signature += to_string((int)(sdv[0].val/filter.sig_bin_width[i]));
			} else {
				signature += "--";
			}
			signature += " : ";
		}
	}

}

//-------------------------------------------------------------------------------------------
// we "simulate" what we did with the event cases, and add n_points_ in a window in a time
// frame before (or after) each control point.

//void add_controls_in_pre_window

//-------------------------------------------------------------------------------------------
void matched_single_to_binary_outcome(MedRepository &rep, OutcomeFilter filter, MedOutcome &mout_s, MedOutcome &mout_b)
{
	int i;

	// general plan:
	// 1) filter the current single events
	// 2) create a signature for each event and count them in a map.
	// 3) randomize samples for control within filtering limit, and then signature and map them as well.
	// 4) go over map values of events and calculate minimal factor in controls
	// 5) decide on a sampling factor (with some adjustments for small numbers) and choose the actual random set of controls.

	MedOutcome mout_f;
	MLOG("Before first event filtering: %d events\n",mout_s.out.size());
	get_filtered_outcome(rep,mout_s,filter,mout_b);
	MLOG("After first event filtering: %d left (ignore %d)\n",mout_b.out.size(),mout_b.ignore_pids.size());
	if (filter.match_flag == 0)
		return;


	MedOutcome mout_c;
	simple_single_to_binary_outcome(rep,mout_b,mout_c,filter.sigs);
	MLOG("After adding initial random controls: %d elements (ignore %d) (n_points %d)\n",mout_c.out.size(),mout_c.ignore_pids.size(),filter.n_points_to_take_in_win);
	filter.filter_just_demographics = 1;
	get_filtered_outcome(rep,mout_c,filter,mout_f);

	MLOG("Filtered outcome : %d elements\n", mout_f.out.size());
	map<string,pair<int,int>> cnts;
	map<string,vector<int>> control_ids;
	string signature;
	pair<int,int> p0;
	vector<int> empty;
	p0.first = 0; p0.second = 0;
	for (i=0; i<mout_f.out.size(); i++) {
		get_outcome_record_signature(rep,filter,mout_f.out[i],signature);
		if (cnts.find(signature) == cnts.end()) {
			cnts[signature] = p0;
			control_ids[signature] = empty;
		}
		if (mout_f.out[i].categ > 0) {
			cnts[signature].first++;
		}
		else {
			cnts[signature].second++;
			control_ids[signature].push_back(i);
		}
	}

	map<string,pair<int,int>>::iterator it;
	float min_ratio = 1e9;
	int n_ev = 0;
	for (it=cnts.begin(); it!=cnts.end(); it++) {
		int ev_cnt = it->second.first;
		int ctrl_cnt = it->second.second;
		float ratio = (float)(1+ctrl_cnt)/(float)(1+ev_cnt);
		if (ev_cnt > 0 && ctrl_cnt > 0) {
			if (ratio < min_ratio && ev_cnt>=filter.match_min_size) min_ratio = ratio;
			n_ev += ev_cnt;
		}

		MLOG("events/control count: %s :: %d <-> %d :: %8.3f\n",it->first.c_str(), ev_cnt, ctrl_cnt,ratio);
	}
	
	int factor = (int)min_ratio;
	if (factor < 1) factor = 1;
	if (factor > filter.match_max_ratio) {
		MLOG("updating factor {%8.3f} to match_max_ratio %d\n", factor, filter.match_max_ratio);
		factor = filter.match_max_ratio;
	}
	MLOG("min match ratio is %8.3f (%d) (%d events , %d expected controls)\n",min_ratio,factor,n_ev,n_ev*factor);

	// now actually choosing a subset in a second pass and building the binary set
	mout_b.type = "binary";
	mout_b.ignore_pids.clear(); // no need for ignores in a binary outcome
	mout_b.pid2use.clear();	    // same...
	// mout_b.out already contains the filtered events, we need to add a sampling of the controls
	for (it=cnts.begin(); it!=cnts.end(); it++) {
		signature = it->first;
		random_shuffle(control_ids[signature].begin(), control_ids[signature].end());
		int ntake = it->second.first * factor;
		if (ntake == 0) ntake = 1 + factor/2;
		if (ntake > (int)control_ids[signature].size()) ntake = (int)control_ids[signature].size();
		MLOG("%s : ntake is %d size is %d %d\n",signature.c_str(),ntake,control_ids[signature].size(),mout_b.out.size());
		for (i=0; i<ntake; i++) {
			mout_b.out.push_back(mout_f.out[control_ids[signature][i]]);
		}
		
	}

	MLOG("Added %d controls to events, with a factor of %d\n",(int)mout_b.out.size() - (int)mout_s.out.size(),factor);
}

void match_binary_outcome(MedRepository &rep, const string &match_params, MedOutcome &mout_in, MedOutcome &mout_out)
{
	OutcomeFilter filter;
	filter.init_from_string(match_params);
	match_binary_outcome(rep, filter, mout_in, mout_out);
}


void try_pairing_ratio(int& opt_cnt1, int& opt_cnt2, double& opt_r, map<string, pair<int, int>> cnts, double w, double r, int match_max_samples) {
	// Find Optimal ratio - cnt1 and cnt2 are the overall number of samples we're giving up on 
	int cnt1 = 0, cnt2 = 0;
	int left_samples = 0;
	for (auto it = cnts.begin(); it != cnts.end(); it++) {
		int ev_cnt = it->second.first;
		int ctrl_cnt = it->second.second;
		if (ev_cnt == 0.0)
			cnt1 += ctrl_cnt;
		else if (ctrl_cnt == 0.0)
			cnt2 += ev_cnt;
		else {
			double iratio = (float)(ctrl_cnt) / (float)(ev_cnt);
			if (iratio < r)
				// more events than we want, have to give up on some:
				cnt2 += (ev_cnt - (int)(ctrl_cnt / r + 0.5));
			else
				cnt1 += (ctrl_cnt - (int)(ev_cnt * r + 0.5));
		}
		left_samples += ev_cnt + ctrl_cnt;
	}
	left_samples -= cnt1;
	left_samples -= cnt2;
	double current_price = cnt1 + w*cnt2;
	double opt_price = opt_cnt1 + w*opt_cnt2;

	MLOG("considered ratio %8.3f price %8.3f lose %d cases and %d controls, left with %d samples\n", r, current_price, cnt2, cnt1, left_samples);
	if ((left_samples < match_max_samples) && (opt_r == 0 || current_price < opt_price)) {
		opt_r = r;
		MLOG("updated opt_r to %8.3f\n", opt_r);
		opt_cnt1 = cnt1;
		opt_cnt2 = cnt2;
	}
}
// search for the optimal ratio between control/event samples
// the price of giving up 1 control is 1.0, the price of giving up 1 event is w 
double get_pairing_ratio(map<string, pair<int, int>> cnts, double w, int n_steps, int match_max_samples) {

	// Get ratio for w->0
	double min_ratio = -1, max_ratio = -1;
	int total_ev = 0, total_ctrl = 0;
	for (auto it = cnts.begin(); it != cnts.end(); it++) {
		int ev_cnt = it->second.first;
		int ctrl_cnt = it->second.second;		
		total_ev += ev_cnt;
		total_ctrl += ctrl_cnt;
		if (ev_cnt > 0 && ctrl_cnt > 0) {
			float iratio = (float)(ctrl_cnt) / (float)(ev_cnt);
			if (min_ratio == -1 || iratio < min_ratio)
				min_ratio = iratio;
			if (max_ratio == -1 || iratio > max_ratio)
				max_ratio = iratio;
		}
	}

	MLOG("min ratio %8.3f max ratio %8.3f match_event_case_price_ratio %8.3f\n", min_ratio, max_ratio, w);
	double r_for_max_samples = max(double(match_max_samples) / (total_ev + 0.5) - 1, 1.0);
	MLOG("total_ctrl %d total_ev %d match_max_samples %d r_for_max_samples %8.3f \n", total_ctrl, total_ev, match_max_samples, r_for_max_samples);

	int opt_cnt1 = -1, opt_cnt2 = -1;
	double opt_r = 0;
	if (r_for_max_samples < max_ratio)
		try_pairing_ratio(opt_cnt1, opt_cnt2, opt_r, cnts, w, r_for_max_samples, match_max_samples);
	double r_step = (max_ratio - min_ratio) / n_steps;
	for (int ir = 0; ir<n_steps; ir++) {
		double r = min_ratio + ir*r_step;
		try_pairing_ratio(opt_cnt1, opt_cnt2, opt_r, cnts, w, r, match_max_samples);
	}
	MLOG("opt ratio %8.3f price %8.3f lose %d cases and %d controls\n", opt_r, opt_cnt1 + w*opt_cnt2, 
		opt_cnt2, opt_cnt1);

	return opt_r;
}


void match_binary_outcome(MedRepository &rep, OutcomeFilter &filter, MedOutcome &mout_in_unfiltered, MedOutcome &mout_out)
{
	MedOutcome mout_in;
	get_filtered_outcome(rep, mout_in_unfiltered, filter, mout_in);
	// now actually choosing a subset in a second pass and building the binary set
	mout_out = mout_in;
	mout_out.ignore_pids.clear(); // no need for ignores in a binary outcome
	mout_out.pid2use.clear();	    // same...
	mout_out.out.clear();

	MLOG("In outcome : before filtering [%d], after filtering [%d] elements\n", mout_in_unfiltered.out.size(), mout_in.out.size());
	map<string, pair<int, int>> cnts;
	map<string, vector<int>> control_ids;
	map<string, vector<int>> event_ids;
	string signature;
	pair<int, int> p0(0,0);
	vector<int> empty;
	for (int i=0; i<mout_in.out.size(); i++) {
		if (filter.min_date < 0 || (mout_in.out[i].date >= filter.min_date && mout_in.out[i].date <= filter.max_date)) {
			get_outcome_record_signature(rep, filter, mout_in.out[i], signature);
			if (cnts.find(signature) == cnts.end()) {
				cnts[signature] = p0;
				control_ids[signature] = empty;
				event_ids[signature] = empty;
			}
			if (mout_in.out[i].categ > 0) {
				cnts[signature].first++;
				event_ids[signature].push_back(i);
			}
			else {
				cnts[signature].second++;
				control_ids[signature].push_back(i);
			}
		}
	}

	double opt_ratio = get_pairing_ratio(cnts, filter.match_event_case_price_ratio, 20, filter.match_max_samples);
	float factor = (float)opt_ratio;
	if (factor < 1) factor = 1;
	if (factor > filter.match_max_ratio) {
		MLOG("updating factor {%8.3f} to match_max_ratio %d\n", factor, filter.match_max_ratio);
		factor = (float)filter.match_max_ratio;
	}
	MLOG("opt ratio is %8.3f (effective factor = %8.3f)\n", opt_ratio, factor);

	
	// sample controls and events into mout
	int n_events = 0, n_ctrl = 0;
	if (!filter.match_do_shuffle_flag)
		MLOG("NOTE: match_do_shuffle_flag=0, taking controls and cases in the order they were given in the input file\n");
	for (auto it=cnts.begin(); it!=cnts.end(); it++) {
		signature = it->first;
		int ev_cnt = it->second.first;
		int ctrl_cnt = it->second.second;
		if (ev_cnt == 0 || ctrl_cnt == 0) {
			MLOG("%s : ntake_ctrl: %d/%d ntake_case: %d/%d total size is %d\n", signature.c_str(), 0, control_ids[signature].size(),
				0, event_ids[signature].size(), mout_out.out.size());			
			continue;
		}
		int ntake_ctrl = (int)((float)ev_cnt * factor);
		//if (ntake == 0) ntake = 1 + (int)factor/2;
		if (ntake_ctrl > (int)control_ids[signature].size()) ntake_ctrl = (int)control_ids[signature].size();

		int ntake_ev = (int)((float)ctrl_cnt / factor);
		//if (ntake == 0) ntake = 1 + (int)factor/2;
		if (ntake_ev > (int)event_ids[signature].size()) ntake_ev = (int)event_ids[signature].size();

		if (filter.match_do_shuffle_flag) {
			random_shuffle(control_ids[signature].begin(), control_ids[signature].end());
			random_shuffle(event_ids[signature].begin(), event_ids[signature].end());
		}
		
		for (int i=0; i<ntake_ctrl; i++)
			mout_out.out.push_back(mout_in.out[control_ids[signature][i]]);
		for (int i = 0; i<ntake_ev; i++)
			mout_out.out.push_back(mout_in.out[event_ids[signature][i]]);
		n_ctrl += ntake_ctrl;
		n_events += ntake_ev;
		MLOG("%s : ntake_ctrl: %d/%d ntake_case: %d/%d total size is %d\n", signature.c_str(), ntake_ctrl, control_ids[signature].size(),
			ntake_ev, event_ids[signature].size(), mout_out.out.size());
	}

	MLOG("Added %d controls, %d cases, with a factor of %f\n", n_ctrl, n_events, n_events, factor);
	mout_out.sort();
}

// for each entry in mout_in, adds ALL the samples in the range of 0-days (with the same categ)
// this is done to enrich samples to contain samples in all periods before the event
void expand_binary_outcome_prefix(MedRepository &rep,int days, const string &sig, MedOutcome &mout_in, MedOutcome &mout_out)
{
	int i,j,len;

	mout_out = mout_in;
	mout_out.out.clear();

	for (i=0; i<mout_in.out.size(); i++) {
		OutcomeRecord outr = mout_in.out[i];
		int ref_days = date_to_days(outr.date);
		mout_out.out.push_back(outr);
		SDateVal *sdv = rep.get_before_date(outr.pid,sig, outr.date, len);
		for (j=len-1; j>=0; j--) {
			if ((ref_days - date_to_days(sdv[j].date)) <= days) {
				outr.date = sdv[j].date;
				mout_out.out.push_back(outr);
			} else {
				break;
			}
		}
	}
}


//....................................................................................................................................................................
//....................................................................................................................................................................
//....................................................................................................................................................................
struct collect_features_thread_info {
	int thread_id;
	int from_out, to_out;
	MedOutcome *mout;
	MedRepository *rep;
	MedFeaturesData *mf;
	vector<string> *sigs;
	vector<string> *snames;
	vector<string> *tnames;
	int state;
	float missing;
};

#endf

//....................................................................................................................................................................
float get_Hosp_days_feature(string &feat_name, MedRepository *rep, int pid, int date)
{
	// format for this feature
	// Hosp_<from win>_<to_win> :: counting the number of hospitalization days within the given window

	vector<string> fields;
	boost::split(fields, feat_name, boost::is_any_of("_"));
	int from_days = stoi(fields[1]);
	int to_days = stoi(fields[2]);
	int len;
	SDateVal *sdv = (SDateVal *)rep->get_before_date(pid, "Hospitalization", date, len);

	//if (len == 0) return 0;
	int i = len - 1;
	int days = date_to_days(date);
	float count = 0;
	while (i >= 0) {
		int idays = date_to_days(sdv[i].date);
		if (days - idays > to_days) break;
		if (days - idays >= from_days)
			count += sdv[i].val;
		i--;
	}
//	MLOG("pid %d : date %d : %s : %f\n", pid, date, feat_name.c_str(), count);
	return count;

}

float get_Smoking_feature(string &feat_name, MedRepository *rep, int pid, int date, float missing)
{
	int current_smoker, ex_smoker;
	int years_since_quitting, smoking_years;
	float pack_years;

	int len;
	bool never_smoked = true;

	SValShort4 *smx_status = (SValShort4 *)rep->get(pid, "SMOKING_ENRICHED", len);
	if (len > 0)
		never_smoked = (smx_status[0].val1 == -1);
	assert(len <= 1);

	if (len == 0) { // No Data
		current_smoker = ex_smoker = (int)missing;
		years_since_quitting = smoking_years = (int)missing;
		pack_years = missing;
	}
	else if (never_smoked) { // Non Smoker
		current_smoker = ex_smoker = 0;
		years_since_quitting = 100;
		smoking_years = 0;
		pack_years = 0.0;
	}
	else { // (Ex)Smoker
		int start_year = smx_status[0].val1;
		int end_year = smx_status[0].val2;
		int target_year = date / 10000;
		if (target_year < end_year) {
			// still in smoking period
			smoking_years = target_year - start_year;
			years_since_quitting = 0;
			current_smoker = 1;
		}
		else {
			// maybe done smoking
			current_smoker = smx_status[0].val4;
			smoking_years = end_year - start_year; // we are merciful
			if (!current_smoker)
				years_since_quitting = target_year - end_year;
			else
				years_since_quitting = 0;
		}
		pack_years = ((float)smx_status[0].val3 / 20) * smoking_years;
		ex_smoker = 1 - current_smoker;
	}

	if (feat_name.compare(17, 100, "Current_Smoker") == 0) return (float) current_smoker;
	else if (feat_name.compare(17, 100, "Ex_Smoker") == 0) return (float) ex_smoker;
	else if (feat_name.compare(17, 100, "Smoking_Years") == 0) return (float) smoking_years;
	else if (feat_name.compare(17, 100, "Smok_Years_Since_quitting") == 0) return (float) years_since_quitting;
	else if (feat_name.compare(17, 100, "Smok_Pack_Years") == 0) return (float) pack_years;
	else {
		MERR("Dont know how to handle: %s", feat_name.c_str());
		throw exception();
	}
}

//....................................................................................................................................................................
float get_Anemia_feature(string &feat_name, MedRepository *rep, int pid, int date)
{
	// format for this feature
	// Anemia_Type1 : 0 - no anemia, 1 - Women: Hg < 11 RBC>=4 , Men: Hg<12.5 RBC>=4.5 , 2 Women Hg<11,RBC<4, Men: Hg<12.5,RBC<4.5 (first approximation before a deeper study/model)

	vector<string> fields;
	boost::split(fields, feat_name, boost::is_any_of("_"));

	int len_hg, len_rbc, len_gen,gender;
	SDateVal *sdv_hg = (SDateVal *)rep->get_before_date(pid, "Hemoglobin", date, len_hg);
	SDateVal *sdv_rbc = (SDateVal *)rep->get_before_date(pid, "RBC", date, len_rbc);
	SVal *sv_gender = (SVal *)rep->get(pid, "GENDER", len_gen);
	gender = (int)sv_gender[0].val;

	if (len_hg == 0 || len_rbc == 0) return 0;

	if (gender == 1) {
		if (sdv_hg[len_hg - 1].val < (float)12.5) {
			if (sdv_rbc[len_rbc - 1].val >= (float)4.5)
				return 1;
			else
				return 2;
		}
		return 0;
	}
	else {
		if (sdv_hg[len_hg - 1].val < (float)11) {
			if (sdv_rbc[len_rbc - 1].val >= (float)4.0)
				return 1;
			else
				return 2;
		}
		return 0;
	}


	return 0;

}
//....................................................................................................................................................................
float get_Dep_days_feature(string &feat_name, MedRepository *rep, int pid, int date)
{
	// format for this feature
	// Dep_<from win>_<to_win>_<Department_name> :: counting the number of hospitalization days within the given window

	vector<string> fields;
	boost::split(fields, feat_name, boost::is_any_of("_"));
	int from_days = stoi(fields[1]);
	int to_days = stoi(fields[2]);
	string dep = fields[3];
	for (int j = 4; j < fields.size(); j++)
		dep += "_" + fields[j];

	int len;
	SDateVal *sdvh = (SDateVal *)rep->get_before_date(pid, "Hospitalization", date, len);
	SDateVal *sdv = (SDateVal *)rep->get_before_date(pid, "Hospitalization_Dep", date, len);

	//if (len == 0) return 0;
	int i = len - 1;
	int days = date_to_days(date);
	float count = 0;
	while (i >= 0) {
		int idays = date_to_days(sdv[i].date);
		if (days - idays > to_days) break;
		if (days - idays >= from_days && rep->dict.is_in_set((int)sdv[i].val,dep))
			count += sdvh[i].val;
		i--;
	}
	//MLOG("pid %d : date %d : %s : %f\n", pid, date, feat_name.c_str(), count);
	return count;

}

//....................................................................................................................................................................
float get_Fit_feature(string &feat_name, MedRepository *rep, int pid, int date)
{
	// format for this feature
	// Fit_Last_<k> or Fit_Max
	float empty = (float)0.001;
	vector<string> fields;
	boost::split(fields, feat_name, boost::is_any_of("_"));
	string type = fields[1];
	int bef = 0;
	if (type == "Last") bef = stoi(fields[2]);

	int len;
	SDateVal *sdv = (SDateVal *)rep->get_before_date(pid, "FECAL_TEST", date, len);

	if (len == 0) return empty;

	if (bef > 0) {
		if (len >= bef) {
			int npos = ((int)sdv[len - bef].val) / 100;
			int ntest = ((int)sdv[len - bef].val) % 100;
			if (ntest > 0)
				return  ((float)npos / (float)ntest);
		}
		return empty;
	}

	// Max case
	if (type == "Max") {
		float max_fit = -1;
		for (int i = 0; i < len; i++) {
			int npos = ((int)sdv[len - bef].val) / 100;
			int ntest = ((int)sdv[len - bef].val) % 100;
			if (ntest > 0) {
				float fit_val = (float)npos / (float)ntest;
				if (fit_val > max_fit)
					max_fit = fit_val;
			}
		}
		return max_fit;
	}

	return empty;

}
//....................................................................................................................................................................
float get_Drug_days_feature(string &feat_name, MedRepository *rep, int pid, int date, map<string, vector<char>>& lut_cache)
{
	// format for this feature
	// Dep_<from win>_<to_win>_<Department_name> :: counting the number of hospitalization days within the given window

	vector<string> fields;
	boost::split(fields, feat_name, boost::is_any_of("_"));
	int min_days_before = stoi(fields[1]);
	int max_days_before = stoi(fields[2]);

	// separate the drug set name, which is everything after fields[2] (and might contain underscores...)
	std::vector<std::string> before_and_after;
	iter_split(before_and_after, feat_name, boost::algorithm::first_finder(fields[2] + "_"));

	if (lut_cache.find(feat_name) == lut_cache.end()) {
		// better to cache a lut (lookup table) on the first pid and reuse it for all pids
		string set_name = before_and_after[1];
		int section_id = rep->dict.section_id("Drug");
		int my_set = rep->dict.dict(section_id)->Name2Id[set_name];
		if (my_set == -1) {
			cerr << "could not find definition for " << set_name << "\n";
			throw exception();
		}
		vector<char> lut;
		rep->dict.prep_sets_lookup_table(section_id, { set_name }, lut);
		lut_cache[feat_name] = lut;
		//fprintf(stderr, "prepared lut for set_id: [%d]\n", my_set);
	}

	int len;
	SDateShort2 *sdv2 = (SDateShort2 *)rep->get(pid, "Drug", len);

//	if (len == 0) return 0;
	int i = len - 1;
	int days = date_to_days(date);
	float count = 0;
	while (i >= 0) {
		if (sdv2[i].date <= date) {
			int idays = date_to_days(sdv2[i].date);			
			if (days - idays >= min_days_before && days - idays <= max_days_before && lut_cache[feat_name][(int)sdv2[i].val1]) {
				count = 1.0;
				break;
			}
		}
		i--;
	}
	//if (count > 0)
	//	MLOG("pid %d : date %d : %s : %s : %f\n", pid, date, feat_name.c_str(), drug.c_str(), count);
	return count;

}

//....................................................................................................................................................................
// threaded - as this part might take a long time
void collect_features_for_outcome_thread(void *p)
{
	collect_features_thread_info *tp = (collect_features_thread_info *)p;

	int i,k,len;
	
	// precalculating signals names and types in order to improve performance (split is not efficient)
	vector<string> snames;
	vector<string> tnames;
	vector<float *> dsigs;
	
	map<string, vector<char>> lut_cache;
//	snames.resize(tp->sigs->size());
//	tnames.resize(tp->sigs->size());
	dsigs.resize(tp->sigs->size());
	
	for (i=0; i<tp->sigs->size(); i++) {
		//vector<string> fields;
		//split(fields, (*tp->sigs)[i], boost::is_any_of("$@<>"));
		//snames[i] = fields[0];
		//tnames[i] = "";
		//if (fields.size()>=2)
		//	tnames[i] = fields[1];
		dsigs[i] = VEC_DATA(tp->mf->data[(*tp->sigs)[i]]);
	}


	// going over outcome and accumulating features
	//MLOG(">>>>>> Starting thread %d from: %d to: %d\n",tp->thread_id,tp->from_out,tp->to_out);
	for (k=tp->from_out; k<=tp->to_out; k++) {
		int pid = tp->mout->out[k].pid;
		if (!(tp->rep->contains_pid(pid))) {
			MERR(">>>> PID %d not in repository\n",pid);
			continue;
		}
		int byear = (int)(((SVal *)tp->rep->get(pid,"BYEAR",len))[0].val);
		float age = get_age(tp->mout->out[k].date,byear);
		int pdate = tp->mout->out[k].date;
		int raw_date = tp->mout->out[k].raw_date;

		// label
		tp->mf->label[k] = tp->mout->out[k].categ;

		//MLOG("k=%d pid %d byear %d age %f pdate %d label %f\n",k,pid,byear,age,pdate, tp->mf->label[k]);

		for (i=0; i<tp->sigs->size(); i++) {

			string sig = (*tp->sigs)[i];
			//float *sigd = VEC_DATA((tp->mf->data[sig]));

//			if (sig.find("$") == string::npos) {
			if ((*tp->tnames)[i] == "") {
				// demographics
				if (sig == "AGE") {
					tp->mf->data["AGE"][k] = age;
				}
				else if (sig == "AGE_INT") {
					tp->mf->data["AGE_INT"][k] = (float)(int)age;
				}
				else if (sig == "Date") {
					tp->mf->data["Date"][k] = (float)pdate;
				}
				else {
					// demographics in which name is exactly as it is in the repository
					float val = ((SVal *)tp->rep->get(pid,sig,len))[0].val;
					dsigs[i][k] = val;
					//tp->mf->data[sig][k] = val;
				}

			} else {


				// date val signals
				string sname, stype;

				// split signal to name and type
				//vector<string> fields;
				//split(fields, sig, boost::is_any_of("$@<>"));
				//sname = fields[0];
				//stype = fields[1];
				sname = (*tp->snames)[i];
				stype = (*tp->tnames)[i];

				if (stype == "Special") {
					// NOTE: if you add a special feature update it here = collect_features AND in MedFeatList creation code
					// In this case we create the feature via a call to a special routine
					if (sname.compare(0, 5, "Hosp_") == 0) dsigs[i][k] = get_Hosp_days_feature(sname, tp->rep, pid, pdate);
					else if (sname.compare(0, 4, "Dep_") == 0) dsigs[i][k] = get_Dep_days_feature(sname, tp->rep, pid, pdate);
					else if (sname.compare(0, 5, "Drug_") == 0) dsigs[i][k] = get_Drug_days_feature(sname, tp->rep, pid, pdate, lut_cache);
					else if (sname.compare(0, 4, "Fit_") == 0) dsigs[i][k] = get_Fit_feature(sname, tp->rep, pid, pdate);
					else if (sname.compare(0, 7, "Anemia_") == 0) dsigs[i][k] = get_Anemia_feature(sname, tp->rep, pid, pdate);
					else if (sname.compare(0, 17, "SMOKING_ENRICHED_") == 0) dsigs[i][k] = get_Smoking_feature(sname, tp->rep, pid, pdate, tp->missing);
					else if (sname.compare(0, 11, "TimeToEvent") == 0) dsigs[i][k] = (float)(date_to_days(raw_date) - date_to_days(pdate));
				}

				int pdate_add = 1;
				/*
				if (sname == "Hospitalization" || sname == "FECAL_TEST" || sname == "Colonscopy") {
					//pdate_add = -10000; // TEST!!! : to make sure we take these questionable items MUCH BEFORE cancer dates
					int m = (pdate % 10000) / 100;
					if (m > 6)
						pdate_add = -600;
					else
						pdate_add = 600 - 10000;
					
				}
				*/
				vector<SDateVal> buffer;
				SDateVal *sdv;				
				if (sname.substr(0, 2) == "BP" || sname.substr(0, 7) == "ALCOHOL") {
					string signal_sname;
					if (sname.substr(0, 2) == "BP")
						signal_sname = "BP";
					else 
						signal_sname = "ALCOHOL";
					int orig_len;
					SDateShort2 *sdv2 = (SDateShort2 *)tp->rep->get(pid, signal_sname, orig_len);
					char which_field = sname.back();
					len = 0;
					for (size_t i = 0; i < orig_len; i++)
					{
						if (sdv2[i].date > pdate + pdate_add)
							break;
						SDateVal newval;
						newval.date = sdv2[i].date;
						if (which_field == '1')
							newval.val = sdv2[i].val1;
						else if (which_field == '2')
							newval.val = sdv2[i].val2;
						else {
							throw runtime_error("field name: " + sname + " must end with '1' or '2'");
						}
						buffer.push_back(newval);
						len++;
					}
					sdv = &buffer[0];
				}
				else {
					sdv = (SDateVal *)tp->rep->get_before_date(pid, sname, pdate + pdate_add, len);
				}

				if (stype == "NVals") {dsigs[i][k] = (float)len; continue;}

				// Last1, Last1T , Min, Max
				if (len > 0) {

					// Last1
					if (stype == "Last1") {dsigs[i][k] = sdv[len-1].val; continue;}

					// Last1T
					if (stype == "Last1T") {dsigs[i][k] = age - get_age(sdv[len-1].date,byear); continue;}

					// Min, Max
					if (stype == "Max") {dsigs[i][k] = sdv_get_max(sdv,len); continue;}
					if (stype == "Min") {dsigs[i][k] = sdv_get_min(sdv,len); continue;}

					// Avg,Std
					if (stype == "Avg") {dsigs[i][k] = sdv_get_avg(sdv,len); continue;}
					if (stype == "Std") {dsigs[i][k] = sdv_get_std(sdv,len); continue;}

					// Before2Y
					if (stype == "Before2Y" || stype == "Delta2Y") {
						int pos = len-1;
						while (pos >= 0) {
							if (get_day_approximate(sdv[pos].date) < get_day_approximate(pdate-720)) {
								if (stype == "Before2Y") dsigs[i][k] = sdv[pos].val;
								if (stype == "Delta2Y") dsigs[i][k] = sdv[len-1].val - sdv[pos].val;
								pos = -1;
							}
							pos--;
						}
						continue;
					}


				}

				// Last1delta , Last2, Last2T

				if (len > 1) {
					// LinVal, LinDelta
					if (stype == "LinVal") {dsigs[i][k] = sdv_get_linear_val(sdv,len-1,pdate); continue;}
					if (stype == "LinDelta") {dsigs[i][k] = sdv_get_linear_delta(sdv,len-1,pdate); continue;}

					// Last1delta
					if (stype == "Last1delta") {dsigs[i][k] = sdv[len-1].val - sdv[len-2].val; continue;}

					// Last2
					if (stype == "Last2") {dsigs[i][k] = sdv[len-2].val; continue;}

					// Last2T
					if (stype == "Last2T") {dsigs[i][k] = age - get_age(sdv[len-2].date,byear); continue;}

				}

				// Slope4Y, Avg4Y
				SDateVal *sdv_s = sdv;
				int jj=0;
				int slen = len;
				while (jj<len && sdv[jj].date < (pdate - 2*1460)) {jj++; sdv_s++; slen--; /*MLOG("j=%d len=%d sate=%d ldate=%d\n",j,len,sdv[j].date,ldate);*/}
				if (slen > 1) {
					//MLOG("2Y: pid %d len %d date0 %d ldate %d slen %d slope %f avg %f\n",pid,len,sdv[0].date,ldate,slen,sdv_get_slope(sdv_s,slen),sdv_get_avg(sdv_s,slen));
					if (stype == "Slope4Y") {dsigs[i][k] = sdv_get_slope(sdv_s,slen); continue;}
					if (stype == "Avg4Y") {dsigs[i][k] = sdv_get_avg(sdv_s,slen); continue;}
				}
				// Slope2Y, Avg2Y
				while (jj<len && sdv[jj].date < (pdate - 1460)) {jj++; sdv_s++; slen--;}
				if (slen > 1) {
					if (stype == "Slope2Y") {dsigs[i][k] = sdv_get_slope(sdv_s,slen); continue;}
					if (stype == "Avg2Y") {dsigs[i][k] = sdv_get_avg(sdv_s,slen); continue;}
				}

			}

		}


	}
	//MLOG(">>>>>> Ending thread %d\n",tp->thread_id);
	tp->state = 1;
}

//....................................................................................................................................................................
// Features should be given in the following way:
// demographic signals: just the name (as appears in the repository)
// date_val signals: given in the form of
int collect_features_for_outcome(MedRepository &rep, MedOutcome &mout, vector<string> &sigs, float missing, MedFeaturesData &mf)
{
	int i;

	mf.clear();
	mf.label_name = mout.name;

	// init mf
	mf.signals = sigs;
	mf.label.resize(mout.out.size());
	for (i=0; i<sigs.size(); i++) {
		mf.data[sigs[i]].resize(mout.out.size());
		fill(mf.data[sigs[i]].begin(),mf.data[sigs[i]].end(),missing);
	}

	mf.row_ids.resize(mout.out.size());
	for (i=0; i<mout.out.size(); i++) {
		mf.row_ids[i].id = mout.out[i].pid;
		mf.row_ids[i].date = mout.out[i].date;
	}

	vector<string> snames(sigs.size());
	vector<string> tnames(sigs.size());

	for (i=0; i<sigs.size(); i++) {
		vector<string> fields;
		split(fields, sigs[i], boost::is_any_of("$@<>"));
		snames[i] = fields[0];
		tnames[i] = "";
		if (fields.size()>=2)
			tnames[i] = fields[1];
	}

	//MLOG("collect_features: After defining names, overall %d features\n",mf.signals.size());

	// thread infos
	int nthreads = 24;
	vector<collect_features_thread_info> tp(nthreads);
	int tlen = (int)mout.out.size()/nthreads;
	for (i=0; i<nthreads; i++) {
		tp[i].thread_id = i;
		tp[i].from_out = i*tlen;
		tp[i].to_out = i*tlen + tlen - 1;
		tp[i].mout = &mout;
		tp[i].rep = &rep;
		tp[i].mf = &mf;
		tp[i].sigs = &sigs;
		tp[i].snames = &snames;
		tp[i].tnames = &tnames;
		tp[i].state = 0;
		tp[i].missing = missing;
	}
	tp[nthreads-1].to_out = (int)mout.out.size()-1;

	// sending threads
	vector<thread> th_handle(nthreads);
	for (int i=0; i<nthreads; i++) {
		th_handle[i] = thread(collect_features_for_outcome_thread, (void *)&tp[i]);
	}

//	MLOG("Waiting for all threads\n");

	int n_state = 0;
	while (n_state < nthreads) {
		this_thread::sleep_for(chrono::milliseconds(10));
		n_state = 0;
		for (int i=0; i<nthreads; i++)
			n_state += tp[i].state;
	}
//	MLOG("Waiting for all threads join()\n");

	for (int i=0; i<nthreads; i++)
		th_handle[i].join();
		

	return 0;
}

//....................................................................................................................
// mode - train/test ... isplit - nsplits: means all (for both train or test)
int get_split_from_outcome(MedOutcome &mout, MedSplit &spl, int isplit, const string &mode, MedOutcome &mout_s)
{
	mout_s.clear();
	mout_s.name = mout.name;
	mout_s.desc = mout.desc + " --> Split " + to_string(isplit);
	mout_s.type = mout.type;
	mout_s.ignore_pids = mout.ignore_pids;
	//MLOG("size in mout: %d , mout_s: %d\n", mout.pid2use.size(), mout_s.pid2use.size());
	//for (auto & pidval : mout.pid2use) {
	//	MLOG("%d %d\n", pidval.first, pidval.second);
	//	//mout_s.pid2use[pidval.first] = pidval.second;
	//}
	//mout_s.pid2use.insert(mout.pid2use.begin(), mout.pid2use.end());

	if (isplit > spl.nsplits) {
		MERR("Wrong split number %d (only %d splits)\n",isplit,spl.nsplits);
		return -1;
	}

	for (int i=0; i<mout.out.size(); i++) {

		//MLOG("i %d\n", i);
		if (spl.pid2split.find(mout.out[i].pid) != spl.pid2split.end()) {

			if ( (isplit == spl.nsplits) ||
				((mode == "all") && (isplit <= spl.nsplits)) ||
				((mode == "train") && (isplit < spl.nsplits) && (spl.pid2split[mout.out[i].pid] != isplit)) ||
				((mode == "train_debug") && (isplit < spl.nsplits) && (spl.pid2split[mout.out[i].pid] > 0) && (spl.pid2split[mout.out[i].pid] != isplit)) ||
				((mode == "test") && (isplit < spl.nsplits) && (spl.pid2split[mout.out[i].pid] == isplit)))

				mout_s.out.push_back(mout.out[i]);

		}
	}

	return 0;
}

//=============================================================================================================================
// Prediction Results handling
//=============================================================================================================================

// formating a line ready to be printed
void PredictedRes::format(string &s)
{
	char fmts[200];

	bool is_int_outcome = (outcome == (float)((int)outcome));
	if (is_int_outcome)
		sprintf(fmts, "%8d\t%10lld\t%8.5f\t%10lld\t%3d\t%3d\t%3d", pid, date, pred, outcome_date, (int)outcome, split, nsplits);
	else
		sprintf(fmts, "%8d\t%10lld\t%8.5f\t%10lld\t%8.5f\t%3d\t%3d", pid, date, pred, outcome_date, outcome, split, nsplits);

	s = string(fmts);

}

// parsing a line containing a record
int PredictedRes::parse(string &line)
{
	int rc = sscanf(line.c_str(), "%d\t%lld\t%f\t%lld\t%f\t%d\t%d", &pid, &date, &pred, &outcome_date, &outcome, &split, &nsplits);

	if (rc != 7)
		return -1;
	return 0;
}

// writing PredictionList to file (as text)
int PredictionList::write_to_file(const string &fname)
{
	sort();	// better to sort first, otherwise bootstrap might get confused
	ofstream of;

	of.open(fname, ios::out);

	if (!of) {
		MERR("PredictionList: can't open file %s for write\n", fname.c_str());
		return -1;
	}

	for (auto &pr : preds) {
		string s;
		pr.format(s);
		of << s << '\n';
	}

	of.close();

	return 0;
}

// reading text PredictionList file
int PredictionList::read_from_file(const string &fname)
{
	ifstream inf(fname);

	if (!inf) {
		MERR("PredictionList: can't open file %s for read\n", fname.c_str());
		return -1;
	}

	preds.clear();
	string curr_line;
	int skipped = 0;
	while (getline(inf, curr_line)) {
		if ((curr_line.size() > 1) && (curr_line[0] != '#')) {
			PredictedRes pr;
			if (pr.parse(curr_line) < 0) {
				MERR("PredictionList: read: wrong format in line : %s\n", curr_line.c_str());
				if (skipped++ > 10)
					return -1;
			} else preds.push_back(pr);

		}
	}

	inf.close();
	return 0;
}

// create a list containing only the given splits
int PredictionList::get_splits(vector<int> splits, PredictionList &pl)
{
	pl.clear();

	vector<int> lut(1024, 0); // assuming less than 1024 splits....
	for (auto sp : splits) lut[sp] = 1;
	for (auto &pr : preds)
		if (pr.split>=0 && lut[pr.split])
			pl.preds.push_back(pr);
	return 0;
}

// get preds and y ready to be used in performance
void PredictionList::get_preds_and_y(vector<float> &predictions, vector<float> &y)
{
	predictions.clear();
	y.clear();
	for (auto &pr: preds) {
		predictions.push_back(pr.pred);
		y.push_back(pr.outcome);
	}
}

// same for a specific split
void PredictionList::get_split_preds_and_y(int split, vector<float> &predictions, vector<float> &y)
{
	predictions.clear();
	y.clear();
	for (auto &pr: preds) {
		if (pr.split == split) {
			predictions.push_back(pr.pred);
			y.push_back(pr.outcome);
		}
	}
}
#endif
