#include "MedEmbed.h"
#include <Logger/Logger/Logger.h>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>                                                                                                                                                
#include <MedUtils/MedUtils/MedGenUtils.h>
#include <MedIO/MedIO/MedIO.h>

#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL	LOG_DEF_LEVEL
using namespace std;

//=====================================================================================
// EmbeddingSig
//=====================================================================================
//-------------------------------------------------------------------------------------
EmbeddedCodeType EmbeddingSig::type_name_to_code(string name)
{
	if (name == "categ" || name == "categorial") return ECTYPE_CATEGORIAL;
	if (name == "cont" || name == "contiuous") return ECTYPE_CONTINUOUS;
	if (name == "age" || name == "AGE" || name == "Age") return ECTYPE_AGE;
	if (name == "dummy" || name == "DUMMY") return ECTYPE_DUMMY;
	if (name == "model" || name == "Model") return ECTYPE_MODEL;
	return ECTYPE_UNDEFINED;
}

//-------------------------------------------------------------------------------------
int EmbeddingSig::init(map<string, string>& _map)
{
	for (auto entry : _map) {
		//MLOG("es init(): %s -> %s\n", entry.first.c_str(), entry.second.c_str());
		string field = entry.first;
		if (field == "sig") sig = entry.second;
		else if (field == "add_hierarchy") add_hierarchy = stoi(entry.second);
		else if (field == "do_shrink") do_shrink = stoi(entry.second);
		else if (field == "do_counts") do_counts = stoi(entry.second);
		else if (field == "time_chan") time_chan = stoi(entry.second);
		else if (field == "val_chan") val_chan = stoi(entry.second);
		else if (field == "win_from") win_from = stoi(entry.second);
		else if (field == "win_to") win_to = stoi(entry.second);
		else if (field == "type") type = type_name_to_code(entry.second);
		else if (field == "sig_time_unit" || field == "time_unit") sig_time_unit = med_time_converter.string_to_type(entry.second);
		else if (field == "win_time_unit") win_time_unit = med_time_converter.string_to_type(entry.second);
		else if (field == "model_file") {
			model_file = entry.second;
			model = new MedModel;
			MLOG("reading model file %s\n", model_file.c_str());
			if (model->read_from_file(model_file) < 0)
				MTHROW_AND_ERR("ERROR: Could not read model %s\n", model_file.c_str());
		}

		else if (field == "ranges") {

			// example ranges=10,100:200,210,220:300,600 will create ranges of (10,100),(200,210),(210,220),(300,600)

			vector<string> f1;
			boost::split(f1, entry.second, boost::is_any_of(":"));

			for (auto &f : f1) {
				vector<string> f2;
				boost::split(f2, f, boost::is_any_of(","));
				if (f2.size()>=2) {
					for (int j=0; j<f2.size()-1; j++) {
						vector<float> r;
						r.push_back(stof(f2[j]));
						r.push_back(stof(f2[j+1]));
						ranges.push_back(r);
					}
				}
			}
		}
		else if (field == "categories") {
			// example categories=list:drug_codes : will create categories of all the sets that are in the file drug_codes
			// here we see them as a comma separated list 
			vector<string> f;
			boost::split(f, entry.second, boost::is_any_of(","));
			categories_to_embed.insert(f.begin(), f.end());
		}

	}

	return 0;
}

//-------------------------------------------------------------------------------------
// given a categorial val : return the sets it is contained in that are in the list of requested categories (= "in range")
int EmbeddingSig::get_categ_orig(int val, vector<int> &members)
{
	if (sig_members2sets_in_range.find(val) == sig_members2sets_in_range.end()) return 0;
	for (auto j : sig_members2sets_in_range[val]) members.push_back(categ_convert[j]);
	return 0;
}

//-------------------------------------------------------------------------------------
// given a categorial val : return all the codes it adds before shrinkage
int EmbeddingSig::get_categ_codes(int val, vector<int> &codes, int use_shrink)
{
	if (use_shrink) get_categ_shrunk_codes(val, codes);
	vector<int> members;
	get_categ_orig(val, members);
	for (auto i : members)	if (Orig2Code.find(i) != Orig2Code.end()) codes.push_back(Orig2Code[i]);
	return 0;
}


//-------------------------------------------------------------------------------------
// given a categorial val : return all the codes it adds after shrinkage
int EmbeddingSig::get_categ_shrunk_codes(int val, vector<int> &codes)
{
	vector<int> members;
	get_categ_orig(val, members);
	for (auto i : members)	if (Orig2ShrunkCode.find(i) != Orig2ShrunkCode.end())	codes.push_back(Orig2ShrunkCode[i]);
	return 0;
}

//-------------------------------------------------------------------------------------
// given a val : get the serial range number it is contained in.
int EmbeddingSig::get_continuous_orig(float val)
{
	for (int i=0; i<ranges.size(); i++) {
		if (val >= ranges[i][0] && val < ranges[i][1])
			return i;
	}

	return -1; // not in range
}

//-------------------------------------------------------------------------------------
// given a val : get the orig (pre shrinking) code for its serial range
int EmbeddingSig::get_continuous_codes(float val, int use_shrink)
{
	if (use_shrink) return get_continuous_shrunk_codes(val);
	int j = get_continuous_orig(val);
	if (j < 0 || (Orig2Code.find(j) == Orig2Code.end())) return -1;
	return Orig2Code[j];
}

//-------------------------------------------------------------------------------------
int EmbeddingSig::get_continuous_shrunk_codes(float val)
{
	int j = get_continuous_orig(val);
	if (j < 0 || (Orig2ShrunkCode.find(j) == Orig2ShrunkCode.end())) return -1;
	return Orig2ShrunkCode[j];
}

//-------------------------------------------------------------------------------------
string EmbeddingSig::print_to_string(int verbosity)
{
	stringstream buffer;

	buffer << "#==================================\n";
	buffer << "# ES description:\n";
	buffer << "#==================================\n";

	buffer << "Sig: " << sig << " type: " << type << "\n";
	buffer << "add_hierarchy: " << add_hierarchy << " do_shrink: " << do_shrink << " channels: " << time_chan << "(t) " << val_chan << "(v)\n";

	buffer << "categories: size " << categories_to_embed.size() << "\n";
	// ranges
	buffer << "Ranges : (" << ranges.size() << ") :\n";
	for (int j=0; j<ranges.size(); j++)
		buffer << "[" << j << "] : " << ranges[j][0] << " - " << ranges[j][1] << "\n";

	// Orig2Code, Orig2Name , Orig2ShrunkCode
	buffer << "Codes: (" << Orig2Code.size() << ") :\n";
	for (auto &e : Orig2Code) {
		int i = e.first;
		buffer << "[" << i << "] code " << Orig2Code[i]; 
		if (Orig2ShrunkCode.find(i) != Orig2ShrunkCode.end()) buffer << " shrunk " << Orig2ShrunkCode[i];
		else buffer << " shrunk ---";
		if (Orig2Name.find(i) != Orig2Name.end()) buffer << " name " << Orig2Name[i];
		else buffer << " name __NO_NAME__";
		buffer << "\n";

	}

	return buffer.str();
}

//---------------------------------------------------------------------------------------------------------------------------
int EmbeddingSig::init_dummy()
{
	MLOG("ES:init_dummy : dummy variable\n");
	Orig2Code[0] = 0;
	Orig2Name[0] = "Dummy_variable_always_1";
	Orig2ShrunkCode[0] = 0; // 0 is always guaranteed to pass through shrinkage
	return 0;
}

//---------------------------------------------------------------------------------------------------------------------------
int EmbeddingSig::init_continous(int &curr_code)
{
	if (type != ECTYPE_CATEGORIAL && type != ECTYPE_AGE) return 0;
	MLOG("ES:init_continous : age/continous sig %s\n", sig.c_str());
	int c = 0;
	for (auto &r : ranges) {
		Orig2Code[c] = curr_code++;
		Orig2Name[c++] = "sig: " + sig + ".t" + to_string(time_chan) + ".v" + to_string(val_chan) + ".win:" + to_string(win_from) + "_" + to_string(win_to) + " range: " + to_string(r[0]) + " - " + to_string(r[1]);
	}
	return 0;
}

//---------------------------------------------------------------------------------------------------------------------------
int EmbeddingSig::init_categorial_tables(MedDictionarySections &dict)
{
	if (type != ECTYPE_CATEGORIAL) return 0;

	int section_id = dict.section_id(sig);

	// first we initialize Name2Id if it is not initialized from serialization
	if (Name2Id.size() == 0) {
		// we only need the names for the needed categories
		for (auto &c : categories_to_embed)
			Name2Id[c] = dict.dicts[section_id].Name2Id[c];
	}

	vector<int> categ_convert(dict.dicts[section_id].Id2Name.rbegin()->first + 1, -1);

	for (auto &c : categories_to_embed) {
		categ_convert[dict.dicts[section_id].Name2Id[c]] = Name2Id[c];
	}

	// init sig_members2sets according to the right add_hierarchy option
	if (add_hierarchy == 0) {
		// in this case each category will only affect its own
		for (auto &e : dict.dicts[section_id].Id2Name)
			sig_members2sets[e.first] = { e.first };
	}
	else {
		vector<int> members;
		dict.dicts[section_id].get_members_to_all_sets(members, sig_members2sets);
	}

	// get sig_members2sets_in_range
	for (auto &e : sig_members2sets) {
		vector<int> in_range;
		for (auto i : e.second)	if (categ_convert[i] >= 0) in_range.push_back(i);
		if (in_range.size() > 0)
			sig_members2sets_in_range[e.first] = in_range;
	}

	// reminder, at this point sig_members2sets and sig_members2sets_in_range are in the numbering set of the apply repository
	// the Orig2X tables are in the numbering of the one used in building the matrix for the first time.
	// Hence to use we have to go through the sig_members2X map, and then convert using categ_convert to use Orig2X.

	return 0;
}

//---------------------------------------------------------------------------------------------------------------------------
int EmbeddingSig::init_categorial(MedDictionarySections &dict, int &curr_code)
{
	if (type != ECTYPE_CATEGORIAL) return 0;
	
	MLOG("ES:init_categorial : model categorial sig %s\n", sig.c_str());

	int section_id = dict.section_id(sig);
	// get Orig2Code, Orig2name
	for (auto &c : categories_to_embed) {
		if (Name2Id.find(c) != Name2Id.end()) {
			int j = Name2Id[c];
			Orig2Code[j] = curr_code++;
			Orig2Name[j] = "sig: " + sig + ".t" + to_string(time_chan) + ".v" + to_string(val_chan) + ".win:" + to_string(win_from) + "_" + to_string(win_to);
			Orig2Name[j] += ".Orig_" + to_string(j);
			for (auto &s : dict.dicts[section_id].Id2Names[j])
				Orig2Name[j] += "|" + s;
		}
	}

	return 0;
}

//-------------------------------------------------------------------------------------------------------------
int EmbeddingSig::get_feat_for_model(MedPidRepository &rep, vector<pair<int, int>> &pids_times)
{
	feat.clear();
	MedSamples samples;

	MLOG("============> 1 pids_times %d\n", pids_times.size());
	for (auto &pt : pids_times)
		samples.insertRec(pt.first, pt.second);
	//samples.time_unit = MedTime::Minutes;
	//global_default_time_unit = MedTime::Minutes;
	//global_default_windows_time_unit = MedTime::Minutes;
	MLOG("============> 2\n");
	samples.normalize();

	MLOG("============> 3 :: Samples %d\n", samples.idSamples.size());

	model->features.clear();
	model->apply(rep, samples);

	MLOG("============> 4\n");

	feat = model->features;

	MLOG("============> 5\n");

	int k = 0;
	for (auto &s : feat.samples) {
		pid_time2idx[pair<int, int>(s.id, s.time)] = k++;
	}

	MLOG("============> 6 \n");

	model->features.clear();
	MLOG("get_feat_for_model : pid_times %d , feat %d x %d , pid_time2idx %d\n", pids_times.size(), feat.data.size(), feat.data.begin()->second.size(), pid_time2idx.size());

	return 0;

}

//----------------------------------------------------------------------------------------------------------------------------------------
int EmbeddingSig::add_sig_to_lines(UniversalSigVec &usv, int pid, int time, int use_shrink, map<int, map<int, float>> &out_lines)
{
	if (out_lines.find(time) == out_lines.end())
		out_lines[time] = map<int, float>();

	if (type == ECTYPE_DUMMY) {	out_lines[time][0] = 1.0f;	return 0; }

	vector<int> codes;
	int from_time = med_time_converter.diff_times(time, sig_time_unit, win_to, win_time_unit, sig_time_unit);
	int to_time = med_time_converter.diff_times(time, sig_time_unit, win_from, win_time_unit, sig_time_unit);
	//MLOG("add_sig_to_lines() pid %d sig %s time %d from_time %d to_time %d usv len: %d\n", pid, es.sig.c_str(), time, from_time, to_time, usv.len);

	if (type == ECTYPE_AGE) {
		// in this case usv must be "BYEAR"
		float age = (float)med_time_converter.get_age(time, sig_time_unit, (int)usv.Val(0));
		codes.push_back(get_continuous_codes(age, use_shrink));
	}

	else if (type == ECTYPE_CATEGORIAL) {

		if (usv.n_time_channels() > 0) {
			for (int j = 0; j<usv.len; j++) {
				int i_time = usv.Time(j, time_chan);
				if (i_time > from_time && i_time <= to_time)
					get_categ_codes((int)usv.Val(j, val_chan), codes, use_shrink); // ??? Consider translating values from current dictionary to originally train dictionary
				
			}
		}
		else {
			// this is for non time signals like GENDER
			for (int j = 0; j<usv.len; j++) 
				get_categ_codes((int)usv.Val(j, val_chan), codes, use_shrink);
		}

	}

	else if (type == ECTYPE_CONTINUOUS) {

		for (int j = 0; j<usv.len; j++) {
			int i_time = usv.Time(j, time_chan);
			if (i_time > from_time && i_time <= to_time)
				codes.push_back(get_continuous_codes(usv.Val(j, val_chan), use_shrink));
		}
	}

	for (auto &c : codes)
		if (c >= 0) {
			if (out_lines[time].find(c) == out_lines[time].end())
				out_lines[time][c] = 0; // initialization of entry
			if (do_counts)
				out_lines[time][c]++;
			else
				out_lines[time][c] = 1.0f;
		}

	return 0;
}


//=====================================================================================
// EmbedMatsCreator
//=====================================================================================



//-------------------------------------------------------------------------------------
// Embedding params init from string
int EmbedMatCreator::init(map<string, string>& _map)
{
	for (auto entry : _map) {
		string field = entry.first;
		if (field == "rep_time_unit") { rep_time_unit = med_time_converter.string_to_type(entry.second); }
		else if (field == "win_time_unit") { rep_time_unit = med_time_converter.string_to_type(entry.second); }
		else if (field == "sigs") {

			// example sigs={sig=Drug;type=categorial;ranges=100000,250000;add_hierarchy=1|sig=Age;type=age;ranges=0,5,18,30,40,50,60,70,80,1000;do_shrink=0}

			vector<string> f;
			boost::split(f, entry.second, boost::is_any_of("|"));
			for (auto &s : f) {
				EmbeddingSig es;
				es.init_from_string(s);
				embed_sigs.push_back(es);
			}
		}
	}

	return 0;
}

//-------------------------------------------------------------------------------------
// prepare :
// (1) prepare the list of signals to load
// (2) read repository on all needed signals and pids
// (3) prepare categ2sets list for every signal
// (4) For each categorial signal calculate:
//     (i) a map from a value to all sets containing it (if needed)
//     (ii) a range of values to map to.
//     (iii) a map from original value to matrix value (before shrinking)
// (5) For each non categorial signal calculate:
//     (i) the number of possible values
//     (ii) values to map to
//     (iii) map for ranges
// (3) add Age signal to codes if needed
//-------------------------------------------------------------------------------------
int EmbedMatCreator::prepare(MedPidRepository &rep)
{
	// initializing sigs_to_load
	sigs_to_load.clear();

	for (auto &es : embed_sigs) {
		if (es.type == ECTYPE_CATEGORIAL || es.type == ECTYPE_CONTINUOUS)
			sigs_to_load.push_back(es.sig);
		if (es.type == ECTYPE_AGE)
			sigs_to_load.push_back("BYEAR"); // special case, assuming BYEAR exists if Age is asked for
		if (es.type == ECTYPE_MODEL) {
			MLOG("Model type prepare() (%s) \n", es.sig.c_str());
			MLOG("model ptr %x\n", es.model);
			vector<string> sigs;
			es.model->get_required_signal_names(sigs);
			MLOG("model sigs (%d) : ", sigs.size());
			for (auto &s : sigs) MLOG("%s,", s.c_str());
			MLOG("\n");

			sigs_to_load.insert(sigs_to_load.end(), sigs.begin(), sigs.end());

			// ??? why is this here?...
			es.model->init_for_apply_rec(rep);
			es.model->features.print_csv();
			rep.sigs.get_sids(sigs, es.model_sids);
			MLOG("sids size %d\n", es.model_sids.size());
		}
	}

	if ((start_sid = rep.sigs.sid("STARTDATE")) > 0) { sigs_to_load.push_back("STARTDATE"); }
	if ((end_sid = rep.sigs.sid("ENDDATE")) > 0) { sigs_to_load.push_back("ENDDATE"); }
	if ((death_sid = rep.sigs.sid("DEATH")) > 0) { sigs_to_load.push_back("DEATH"); }

	// preparing coding space (pre shrinking) , and relevant maps

	curr_code = 1; // 0 kept for dummy
	for (auto &es : embed_sigs) {
		
		// prepare sig_members2sets , sig_members2sets_in_range, Orig2Code , Orig2Name and increase curr_code
		prep_memebers_to_sets(rep, es);

	}

	return 0;

}


//-------------------------------------------------------------------------------------
void EmbedMatCreator::prep_memebers_to_sets(MedPidRepository &rep, EmbeddingSig &es)
{

	es.sig_members2sets.clear();
	es.sig_members2sets_in_range.clear();
	es.Orig2Name.clear();
	es.Orig2Code.clear();
	es.Orig2ShrunkCode.clear();
	es.Name2Id.clear();

	if (es.type == ECTYPE_DUMMY) es.init_dummy(); 

	if (es.type == ECTYPE_MODEL) {

		MLOG("prep_memebers_to_sets : model variable\n");
		es.do_shrink = 0; // we never shrink model features

		// first get an empty features file with all the given names
		es.model->features.clear();
		MedSamples samples;
		es.model->verbosity = 0;
		es.model->apply(rep, samples, MED_MDL_APPLY_FTR_GENERATORS, MED_MDL_APPLY_FTR_PROCESSORS);
		int c = 0;
		for (auto &feat : es.model->features.data) {
			es.Orig2Code[c] = curr_code++;
			es.Orig2Name[c] = "From Model: " + feat.first;
			es.Orig2ShrunkCode[c] = es.Orig2Code[c];
			c++;
		}
		es.model->features.clear();

		return;
	}

	if (es.type == ECTYPE_AGE || es.type == ECTYPE_CONTINUOUS) es.init_continous(curr_code);

	if (es.type == ECTYPE_CATEGORIAL) es.init_categorial(rep.dict, curr_code);

}

//---------------------------------------------------------------------------------------------------------------------------------------------------
int EmbedMatCreator::add_model_feats_to_lines(EmbeddingSig &es, PidDynamicRec &pdr, vector<int> &times, int use_shrink, map<int, map<int, float>> &out_lines)
{
	for (auto time : times)
		if (out_lines.find(time) == out_lines.end())
			out_lines[time] = map<int, float>();

	MedIdSamples mis;
	mis.id = pdr.pid;
	for (auto t : times) {
		MedSample s;
		s.id = pdr.pid;
		s.time = t;
		mis.samples.push_back(s);
	}
	MedFeatures _feat;
	es.model->apply_rec(pdr, mis, _feat, true);
	//_feat.print_csv();
	assert(_feat.data.size() == es.Orig2Name.size());

	int c = 0;
	for (auto &f : _feat.data) {

		int code = c;
		if (use_shrink) code = es.Orig2ShrunkCode[c];

		for (int i=0; i<_feat.samples.size(); i++) {
			int time = _feat.samples[i].time;
			if (out_lines[time].find(code) == out_lines[time].end())
				out_lines[time][code] = 0; // initialization of entry
			out_lines[time][code] = f.second[i];
		}

		c++;
	}

	return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------
int EmbedMatCreator::add_model_feats_to_lines(EmbeddingSig &es, int pid, vector<int> &times, int use_shrink, map<int, map<int, float>> &out_lines)
{
	for (auto time : times)
		if (out_lines.find(time) == out_lines.end())
			out_lines[time] = map<int, float>();

	for (auto time : times) {
		pair<int, int> p(pid, time);
		if (es.pid_time2idx.find(p) == es.pid_time2idx.end()) {
			MLOG("ERROR: could not find mat for %d %d\n", pid, time);
			continue;
		}

		int i = es.pid_time2idx[p];

		//MLOG("add feats : pid %d time %d i %d (mat %d x %d) \n", pid, time, i, es.feat.data.size(), es.feat.data.begin()->second.size());

		int c = 0;
		for (auto &f : es.feat.data) {
			int code = c;
			if (use_shrink) code = es.Orig2ShrunkCode[c];
			int time = es.feat.samples[i].time;
			if (out_lines[time].find(code) == out_lines[time].end())
				out_lines[time][code] = 0; // initialization of entry
			out_lines[time][code] = f.second[i];
			c++;
		}
	}

	return 0;
}


//--------------------------------------------------------------------------------------------------------------------
int EmbedMatCreator::add_pid_lines(PidDynamicRec &pdr, MedSparseMat &smat, vector<int> &times, int use_shrink)
{
	// if n_versions is 0 : we use just version 0 (original)
	// otherwise n_versions must match the number of times we have, and we use each version.
	if (pdr.get_n_versions() != 0 && pdr.get_n_versions() != times.size()) {
		MTHROW_AND_ERR("ERROR: Can't create lines with pdr: n_versions is %d , while ntimes is %d\n", pdr.get_n_versions(), (int)times.size());
		return -1;
	}


	// prepare out_lines using embed_sigs mechanisms

	map<int, map<int, float>> out_lines;
	UniversalSigVec usv;
	for (auto &es : embed_sigs) {

		if (es.type == ECTYPE_MODEL) {
			add_model_feats_to_lines(es, pdr, times, use_shrink, out_lines);
		}
		else {
			for (int t=0; t<times.size(); t++) {
				int ver = t;
				if (pdr.get_n_versions() == 0) ver = 0;

				if (es.type != ECTYPE_DUMMY) pdr.uget(es.sid, ver, usv);

				es.add_sig_to_lines(usv, pdr.pid, times[t], use_shrink, out_lines);
			}
		}

	}


	// push to smat - must be critical to ensure thread safeness
#pragma omp critical
	{
		SparseMatRowMetaData meta;
		meta.pid = pdr.pid;
		for (int t=0; t<times.size(); t++) {
			meta.time = times[t];
			smat.insert_line(meta, out_lines[times[t]]);
		}
	}

	return 0;
}

//-----------------------------------------------------------------------------
void EmbedMatCreator::init_sids(MedPidRepository &rep)
{
	start_sid = rep.sigs.sid("STARTDATE");
	end_sid = rep.sigs.sid("ENDDATE");
	death_sid = rep.sigs.sid("DEATH");
	for (auto &es : embed_sigs)
		es.sid = rep.sigs.sid(es.sig);
}


//--------------------------------------------------------------------------------------------
// calculating the sub columns we need to produce by throwing away the columns that have too
// few or too much elements.
int EmbedMatCreator::get_shrinked_dictionary(MedSparseMat &smat, float min_p, float max_p)
{
	int max_code = 0;
	for (auto &es : embed_sigs) {
		if (es.Orig2Code.rbegin()->second > max_code)
			max_code = es.Orig2Code.rbegin()->second; // ??? who says max code is in the last one since we moved to go over categories.
	}

	vector<int> code2count(max_code+1, 0);

	for (auto &line : smat.lines) {
		for (auto & e : line)
			code2count[e.first]++;
	}
	int nlines = (int)smat.lines.size();
	if (nlines == 0) nlines = 1;

	int curr_code = 1;
	int n_before = 0, n_after = 0;
	for (auto &es : embed_sigs) {

		if (es.type == ECTYPE_DUMMY) continue;

		es.Orig2ShrunkCode.clear();
		if (es.do_shrink) {
			// passing only those codes that appear in the min_p:max_p range
			for (auto &e : es.Orig2Code) {
				int j = e.second;
				float p = (float)code2count[j]/(float)nlines;
				if (p>=min_p && p<=max_p) {
					es.Orig2ShrunkCode[e.first] = curr_code++;
				}
			}

		}
		else {
			// passing all codes to next stage
			for (auto &e : es.Orig2Code)
				es.Orig2ShrunkCode[e.first] = curr_code++;
		}

		MLOG("Shrinking %s : %d -> %d\n", es.sig.c_str(), es.Orig2Code.size(), es.Orig2ShrunkCode.size());
		n_before += (int)es.Orig2Code.size();
		n_after += (int)es.Orig2ShrunkCode.size();

	}

	MLOG("Total shrinkage: %d -> %d\n", n_before, n_after);

	return 0;
}

//--------------------------------------------------------------------------------------------
int EmbedMatCreator::shrink_mat(MedSparseMat &smat, MedSparseMat &shrunk_smat)
{
	// first we need to calculate a reverse from codes to shrunk codes
	map<int, int> code2shrunk;
	for (auto &es: embed_sigs) {
		for (auto &e : es.Orig2ShrunkCode)
			code2shrunk[es.Orig2Code[e.first]] = e.second;
	}

	shrunk_smat.clear();
	for (int i=0; i<smat.meta.size(); i++) {
		vector<pair<int, float>> new_line;
		for (auto &e : smat.lines[i]) {
			if (code2shrunk.find(e.first) != code2shrunk.end())
				new_line.push_back(pair<int, float>(code2shrunk[e.first], (float)e.second));
		}
		if (new_line.size() > 0) {
			shrunk_smat.meta.push_back(smat.meta[i]);
			shrunk_smat.lines.push_back(new_line);
		}
	}
	
	return 0;
}


//--------------------------------------------------------------------------------------------
int EmbedMatCreator::write_dict_to_file(string fname, int only_shrink)
{
	ofstream dict_f;

	dict_f.open(fname);

	if (!dict_f.is_open()) {
		MERR("Can't open file %s for writing\n", fname.c_str());
		return -1;
	}

	if (only_shrink) {
		dict_f << "shrunk_val\torig_val\tname\n";
		for (auto &es : embed_sigs) {
			for (auto &e : es.Orig2ShrunkCode) {
				dict_f << e.second << "\t" << e.first << "\t" << es.Orig2Name[e.first] << "\n";
			}
		}
	}
	else {
		dict_f << "code\tshrunk_val\torig_val\tname\n";
		for (auto &es : embed_sigs) {
			for (auto &e : es.Orig2Code) {
				int shrunk = -1;
				if (es.Orig2ShrunkCode.find(e.first) != es.Orig2ShrunkCode.end()) shrunk = es.Orig2ShrunkCode[e.first];
				dict_f << e.second << "\t" << shrunk << "\t" << e.first << "\t" << es.Orig2Name[e.first] << "\n";
			}
		}
	}

	return 0;
}

//--------------------------------------------------------------------------------------------
int EmbedMatCreator::get_sparse_mat(MedPidRepository &rep, MedSamples &samples, int use_outcome_time, int use_shrink, MedSparseMat &smat)
{
	vector<MedSample> samples_vec;

	samples.export_to_sample_vec(samples_vec);

	vector<pair<int, int>> pids_times;
	for (auto &s : samples_vec)
		if (use_outcome_time)
			pids_times.push_back(pair<int, int>(s.id, s.outcomeTime));
		else
			pids_times.push_back(pair<int, int>(s.id, s.time));

	return get_sparse_mat(rep, pids_times, use_shrink, smat);
}

//--------------------------------------------------------------------------------------------
int EmbedMatCreator::get_sparse_mat(MedPidRepository &rep, vector<pair<int, int>> &pids_times, int use_shrink, MedSparseMat &smat)
{
	MLOG("Generating sparse mat for %d lines\n", pids_times.size());
	
	smat.clear();
	smat.lines.resize(pids_times.size());
	smat.meta.resize(pids_times.size());

	int n = 0;
	//MedRepository *p_rep = (MedRepository *)&rep;

	for (auto &es : embed_sigs) {
		if (es.type == ECTYPE_MODEL) {
			MLOG("==============================> Generating model part of matrix\n");
			es.get_feat_for_model(rep, pids_times);
		}
	}

#pragma omp parallel for
	for (int i=0; i<pids_times.size(); i++) {

		int pid = pids_times[i].first;
		int time = pids_times[i].second;

		map<int, map<int, float>> out_lines;
		UniversalSigVec usv;
		for (auto &es : embed_sigs) {
			if (es.type == ECTYPE_MODEL) {
				//PidDynamicRec drec;
				vector<int> times = { time };
				//MLOG("============= %d  pid %d \n", es.model_sids.size(), pid);
				//drec.init_from_rep(&rep, pid, es.model_sids, times.size());
				//MLOG("Generated drec for pid %d time %d sids %d size %d\n", pid, time, es.model_sids.size(), drec.data_len);
				//add_model_feats_to_lines(es, drec, times, use_shrink, out_lines);
				add_model_feats_to_lines(es, pid, times, use_shrink, out_lines);
			}
			else {
				if (es.type != ECTYPE_DUMMY)	rep.uget(pid, es.sid, usv);
				es.add_sig_to_lines(usv, pid, time, use_shrink, out_lines);
			}
		}

		vector<pair<int, float>> line;
		smat.convert_map_to_line(out_lines[time], line);
		SparseMatRowMetaData meta;
		meta.time = time;
		meta.pid = pid;
		smat.insert_line(meta, i, line);

#pragma omp critical
		{
			n++;
			if (n % 10000 == 0) MLOG("Generated %d lines\n", n);
		}

	}

	return 0;
}


//--------------------------------------------------------------------------------------------
string EmbedMatCreator::print_to_string(int verbosity)
{
	string sout;

	sout += "sigs_to_load: ";
	for (auto &s : sigs_to_load)
		sout += s + ",";
	sout += "\n";
	sout += "embed_sigs: " + to_string(embed_sigs.size());
	for (auto &es : embed_sigs) {
		sout += es.print_to_string(verbosity);
	}

	return sout;
}



//-----------------------------------------------------------------------------
int EmbedTrainCreator::init(map<string, string>& _map)
{
	for (auto entry : _map) {
		string field = entry.first; 
		if (field == "rep_time_unit") { rep_time_unit = med_time_converter.string_to_type(entry.second); }
		if (field == "use_same_dictionaries") { use_same_dictionaries = stoi(entry.second); }
		else if (field == "win_time_unit") { win_time_unit = med_time_converter.string_to_type(entry.second); }
		else if (field == "byear_time_unit") { byear_time_unit = med_time_converter.string_to_type(entry.second); }
		else if (field == "min_time") { min_time = stoi(entry.second); }
		else if (field == "max_time") { max_time = stoi(entry.second); }
		else if (field == "min_age") { min_age = stoi(entry.second); }
		else if (field == "max_age") { max_age = stoi(entry.second); }
		else if (field == "npoints_per_pid") { npoints_per_pid = stoi(entry.second); }
		else if (field == "min_p") { min_p = stof(entry.second); }
		else if (field == "max_p") { max_p = stof(entry.second); }
		else if (field == "p_train") { p_train = stof(entry.second); }
		else if (field == "x_params") { x_params = entry.second; }
		else if (field == "y_params") { y_params = entry.second; }
		else if (field == "time_dist_range") {
			vector<string> f;
			boost::split(f, entry.second, boost::is_any_of(","));
			if (f.size() >= 2) {
				for (auto &e : f) time_dist_range.push_back(stoi(e));
			}
		}
		else if (field == "time_dist_points") {
			vector<string> f;
			boost::split(f, entry.second, boost::is_any_of(","));
			for (auto &e : f) time_dist_points.push_back(stoi(e));
		}

	}

	return 0;
}

// helpers: read/write a file of <pid> <xtime> <ytime> records
//---------------------------------------------------------------------------------------
int EmbedTrainCreator::read_xy_records(string xy_fname, vector<EmbedXYRecord> &xy)
{
	vector<vector<string>> res;
	if (read_text_file_cols(xy_fname, " \t", res) < 0) return -1;
	if (res.size() >= 0) {
		for (int i=0; i<res.size(); i++) {			
			if (res[i].size() >= 3) {
				EmbedXYRecord rec;
				rec.pid = stoi(res[i][0]);
				rec.x_time = stoi(res[i][1]);
				rec.y_time = stoi(res[i][2]);
				xy.push_back(rec);
			}
		}
	}
	return 0;
}

//----------------------------------------------------------------------------------------
int EmbedTrainCreator::write_xy_records(string xy_fname, vector<EmbedXYRecord> &xy)
{
	ofstream xy_f;

	xy_f.open(xy_fname);

	if (!xy_f.is_open()) {
		MERR("Can't open file %s for writing\n", xy_fname.c_str());
		return -1;
	}

	for (auto &rec : xy) {
		xy_f << rec.pid << " " << rec.x_time << " " << rec.y_time << "\n";
	}

	return 0;
}

//----------------------------------------------------------------------------------------
int EmbedTrainCreator::generate_from_xy_file(string xy_fname, string rep_fname, string out_prefix)
{
	MLOG("Generating an xy training sparse matrices for train/test : xy %s , rep %s\n", xy_fname.c_str(), rep_fname.c_str());

	// init repository
	MedPidRepository rep;

	if (rep.init(rep_fname) < 0) return -1;
	//if (rep.read_all(rep_fname) < 0) return -1;
	// initialize mat creators
	EmbedMatCreator x_emc, y_emc;
	if (y_params == "") y_params = x_params;
	x_emc.init_from_string(x_params);
	y_emc.init_from_string(y_params);
	x_emc.prepare(rep);
	y_emc.prepare(rep);
	x_emc.init_sids(rep);
	y_emc.init_sids(rep);


	x_emc.rep_time_unit = rep_time_unit;
	x_emc.win_time_unit = win_time_unit;
	y_emc.rep_time_unit = rep_time_unit;
	y_emc.win_time_unit = win_time_unit;

	// init lists
	vector<EmbedXYRecord> xy_recs;
	if (read_xy_records(xy_fname, xy_recs) < 0) return -1;

	vector<pair<int, int>> x_list_train, y_list_train, x_list_test, y_list_test;
	vector<int> pids;
	unordered_map<int, int> pids2group;
	for (auto &r : xy_recs) 
		if (pids2group.find(r.pid) == pids2group.end()) {
			if (rand_1() <= p_train)
				pids2group[r.pid] = 1;
			else
				pids2group[r.pid] = 0;
			pids.push_back(r.pid);
		}

	for (auto &r : xy_recs) {
		if (pids2group[r.pid]) {
			x_list_train.push_back(pair<int, int>(r.pid, r.x_time));
			y_list_train.push_back(pair<int, int>(r.pid, r.y_time));
		}
		else {
			x_list_test.push_back(pair<int, int>(r.pid, r.x_time));
			y_list_test.push_back(pair<int, int>(r.pid, r.y_time));
		}
	}

	MLOG("pids %d , p_train %f, x_list train %d, test %d, y_list train %d, test %d\n",
		pids.size(), p_train, x_list_train.size(), x_list_test.size(), y_list_train.size(), y_list_test.size());

	// load repository
	vector<string> sigs;
	sigs = x_emc.sigs_to_load;
	sigs.insert(sigs.end(), y_emc.sigs_to_load.begin(), y_emc.sigs_to_load.end());
	MLOG("Reading into rep : %d pids x %d sigs\n", pids.size(), sigs.size());
	if (rep.load(sigs, pids) < 0) MTHROW_AND_ERR("Failed reading pids x sigs from repository\n"); 

	//MLOG("x_emc:\n%s", x_emc.print_to_string(1).c_str());

	// handling the shrinkage first
	MedSparseMat x_train, y_train, x_test, y_test;

	// first calculating shrinkage on x_train group, which we need anyway
	MLOG("get x train mat without shrinkage\n");
	x_emc.get_sparse_mat(rep, x_list_train, 0, x_train);
	MLOG("Calculate x shrinked dictionary for probs %f-%f\n", min_p, max_p);
	x_emc.get_shrinked_dictionary(x_train, min_p, max_p);

	if (use_same_dictionaries) {
		MLOG("Use x shrinking in y as well\n");
		if (x_emc.embed_sigs.size() != y_emc.embed_sigs.size())
			MTHROW_AND_ERR("Can't use same dictionaries on x,y : different es sizes\n");
		for (int i=0; i<x_emc.embed_sigs.size(); i++) {
			if ((x_emc.embed_sigs[i].sig != y_emc.embed_sigs[i].sig) ||
				(x_emc.embed_sigs[i].type != y_emc.embed_sigs[i].type) ||
				(x_emc.embed_sigs[i].ranges.size() != y_emc.embed_sigs[i].ranges.size()) ||
				(x_emc.embed_sigs[i].Orig2Code.size() != y_emc.embed_sigs[i].Orig2Code.size()) ||
				(x_emc.embed_sigs[i].Orig2Name.size() != y_emc.embed_sigs[i].Orig2Name.size()))
				MTHROW_AND_ERR("x,y EmbedMatCreators do not match, can't use same dictionaries (i=%d)\n", i);
			y_emc.embed_sigs[i].Orig2ShrunkCode = x_emc.embed_sigs[i].Orig2ShrunkCode;
		}
	}
	else {
		MLOG("get y train mat without shrinkage\n");
		y_emc.get_sparse_mat(rep, y_list_train, 0, y_train);
		MLOG("Calculate x shrinked dictionary\n");
		y_emc.get_shrinked_dictionary(y_train, min_p, max_p);
	}


	// generating matrices
	x_train.clear();
	y_train.clear();

	MLOG("Generating x_train matrix (%d)\n", x_list_train.size());
	x_emc.get_sparse_mat(rep, x_list_train, 1, x_train);
	MLOG("writing x_train matrix (%d)\n", x_list_train.size());
	x_train.write_to_files(out_prefix+"_xtrain.mat", out_prefix+"_xtrain.meta", "");

	if (x_list_test.size() > 0) {
		MLOG("Generating x_test matrix (%d)\n", x_list_test.size());
		x_emc.get_sparse_mat(rep, x_list_test, 1, x_test);
		MLOG("writing x_test matrix (%d)\n", x_list_test.size());
		x_test.write_to_files(out_prefix+"_xtest.mat", out_prefix+"_xtest.meta", "");
	}

	MLOG("Generating y_train matrix (%d)\n", y_list_train.size());
	y_emc.get_sparse_mat(rep, y_list_train, 1, y_train);
	MLOG("writing y_train matrix (%d)\n", y_list_train.size());
	y_train.write_to_files(out_prefix+"_ytrain.mat", out_prefix+"_ytrain.meta", "");

	if (y_list_test.size() > 0) {
		MLOG("Generating y_test matrix (%d)\n", y_list_test.size());
		y_emc.get_sparse_mat(rep, y_list_test, 1, y_test);
		MLOG("writing y_test matrix (%d)\n", y_list_test.size());
		y_test.write_to_files(out_prefix+"_ytest.mat", out_prefix+"_ytest.meta", "");
	}

	MLOG("Writing schemes and dictionaries");
	x_emc.write_to_file(out_prefix + "_x.scheme");
	x_emc.write_dict_to_file(out_prefix + "_x.dict", 0);
	x_emc.write_dict_to_file(out_prefix + "_x_shrunk.dict", 1);
	y_emc.write_to_file(out_prefix + "_y.scheme");
	y_emc.write_dict_to_file(out_prefix + "_y.dict", 0);
	y_emc.write_dict_to_file(out_prefix + "_y_shrunk.dict", 1);

	return 0;
}

//-------------------------------------------------------------------------------------------------
int EmbedTrainCreator::generate_xy_list(string xy_fname, string pids_fname, string rep_fname)
{
	vector<int> pids;
	vector<string> spids;
	if (read_text_file_col(pids_fname, "#", " \t,;", 0, spids) < 0) MTHROW_AND_ERR("Error reading pids file %s\n", pids_fname.c_str());
	for (auto &s : spids) pids.push_back(stoi(s));
	MLOG("read %d pids from file %s\n", pids.size(), pids_fname.c_str());

	MedPidRepository rep;
	if (rep.init(rep_fname) < 0) MTHROW_AND_ERR("Error initializing repository %s\n", rep_fname.c_str());

	vector<string> sigs;
	int byear_sid, start_sid, end_sid, death_sid;

	if ((byear_sid = rep.sigs.sid("BYEAR")) > 0) sigs.push_back("BYEAR");
	if ((start_sid = rep.sigs.sid("STARTDATE")) > 0) sigs.push_back("STARTDATE");
	if ((end_sid = rep.sigs.sid("ENDDATE")) > 0) sigs.push_back("ENDDATE");
	if ((death_sid = rep.sigs.sid("DEATH") > 0)) sigs.push_back("DEATH");

	MLOG("sids byear %d start %d end %d death %d\n", byear_sid, start_sid, end_sid, death_sid);
	if (rep.load(sigs, pids) < 0) return -1;

	vector<EmbedXYRecord> xy_times;

	for (auto pid : pids) {

		int from_time = min_time;
		int to_time = max_time;

		if (start_sid > 0) {
			int start_t = medial::repository::get_value(rep, pid, start_sid);
			if (start_t > from_time) from_time = start_t;
		}

		if (end_sid > 0) {
			int end_t = medial::repository::get_value(rep, pid, end_sid);
			if (end_t < to_time) to_time = end_t;
		}

		
		if (byear_sid > 0) {
			// take age into account
			int byear = medial::repository::get_value(rep, pid, byear_sid);
			if (byear_time_unit == MedTime::Years) byear -= 1900;
			int min_age_t = med_time_converter.convert_times(byear_time_unit, rep_time_unit, byear+min_age);
			int max_age_t = med_time_converter.convert_times(byear_time_unit, rep_time_unit, byear+max_age);
			//MLOG("byear=%d min_age %d max_age %d min_age_t %d max_age_t %d\n", byear, min_age, max_age, min_age_t, max_age_t);

			if (from_time < min_age_t) from_time = min_age_t;
			if (to_time > max_age_t) to_time = max_age_t;
		}
		
		//MLOG("pid %d from %d to %d\n", pid, from_time, to_time);

		// at this stage we have to choose time points between from_time and to_time

		if (to_time > from_time) {

			for (int i=0; i<npoints_per_pid; i++) {

				// choose a random time :
				int x_rtime, y_rtime;

				// date case
				if (rep_time_unit == MedTime::Date) {
					int from_days = med_time_converter.convert_date(MedTime::Days, from_time);
					int to_days = med_time_converter.convert_date(MedTime::Days, to_time);
					int rdays = from_days + rand_N(to_days - from_days + 1);
					x_rtime = med_time_converter.convert_days(MedTime::Date, rdays);
				}
				else {
					// general case
					x_rtime = from_time + rand_N(to_time - from_time + 1);					
				}

				// choose a random time_dist and get y_rtime
				int yd = 0;
				if ((time_dist_range.size() > 0) && (time_dist_points.size() > 0)) {
					if (rand_1() < 0.5f)
						yd = time_dist_range[0] + rand_N(time_dist_range[1] - time_dist_range[0] + 1);
					else
						yd = time_dist_points[rand_N((int)time_dist_points.size())];
				}
				else if (time_dist_range.size() > 0) {
					yd = time_dist_range[0] + rand_N(time_dist_range[1] - time_dist_range[0] + 1);
				} else if (time_dist_points.size() > 0)
					yd = time_dist_points[rand_N((int)time_dist_points.size())];

				if (rep_time_unit == MedTime::Date) {
					y_rtime = med_time_converter.add_subtruct_days(x_rtime, yd);
				}
				else
					y_rtime = x_rtime + yd;

				EmbedXYRecord r;
				r.pid = pid;
				r.x_time = x_rtime;
				r.y_time = y_rtime;
				xy_times.push_back(r);
			}
		}
	}

	MLOG("Writing list file %s with %d records\n", xy_fname.c_str(), xy_times.size());
	write_xy_records(xy_fname, xy_times);
	return 0;
}
