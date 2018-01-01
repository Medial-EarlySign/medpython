#include "MedSamples.h"
#include "MedProcessTools/MedProcessTools/MedFeatures.h"
#include "Logger/Logger/Logger.h"

#define LOCAL_SECTION MED_SAMPLES_CV
#define LOCAL_LEVEL	LOG_DEF_LEVEL

//=======================================================================================
// MedSample
//=======================================================================================
// Get sample from tab-delimited string, where pos indicate the position of each field (fields are id,date,outcome,outcome_date,split,pred)
//.......................................................................................
int MedSample::parse_from_string(string &s, map <string, int> & pos) {
	if (pos.size() == 0)
		return parse_from_string(s);
	vector<string> fields; 
	boost::split(fields, s, boost::is_any_of("\t\n\r"));
	if (fields.size() == 0)
		return -1;
	try {
		if (pos["id"] != -1)
			id = (int)stod(fields[pos["id"]]);
		if (pos["date"] != -1)
			time = (int)stod(fields[pos["date"]]);
		if (pos["outcome"] != -1)
			outcome = stof(fields[pos["outcome"]]);
		if (pos["outcome_date"] != -1)
			outcomeTime = (int)stod(fields[pos["outcome_date"]]);
		if (pos["split"] != -1 && fields.size() > pos["split"])
			split = stoi(fields[pos["split"]]);
		if (pos["pred"] != -1 && fields.size() > pos["pred"])
			prediction.push_back(stof(fields[pos["pred"]]));
		return 0;
	}
	catch (std::invalid_argument e) {
		MLOG("could not parse [%s]\n", s.c_str());
		throw e;
	}
}

// Get sample from tab-delimited string, in old or new format (<split> and <prediction> optional, <predictions> can be several numbers (tab delimited))
//.......................................................................................
int MedSample::parse_from_string(string &s)
{
	vector<string> fields;
	boost::split(fields, s, boost::is_any_of("\t"));

	// old format is starting with EVENT
	// new format is starting with SAMPLE
	prediction.clear();

	// old format:
	// EVENT <id> <time> <outcome> <outcomeLen(dummy here)> <outcomeTime> <split> <predictions>
	// <split> and <prediction> optional
	// <predictions> can be several numbers (tab delimited)
	if (fields[0] == "EVENT") {
		if (fields.size() < 6) return -1;
		id = stoi(fields[1]);
		time = stoi(fields[2]);
		outcome = stof(fields[3]);
		int dummy_length = stoi(fields[4]);
		outcomeTime = stoi(fields[5]);
		if (fields.size() >= 7)
			split = stoi(fields[6]);
		if (fields.size() >= 8) {
			for (int i=7; i<fields.size(); i++)
				prediction.push_back(stof(fields[i]));
		}
		return 0;
	}

	// new format:
	// SAMPLE <id> <time> <outcome> <outcomeTime> <split> <predictions>
	// <split> and <prediction> optional
	// <predictions> can be several numbers (tab delimited)
	if (fields[0] == "SAMPLE") {
		if (fields.size() < 5) return -1;
		id = stoi(fields[1]);
		time = stoi(fields[2]);
		outcome = stof(fields[3]);
		outcomeTime = stoi(fields[4]);
		if (fields.size() >= 6)
			split = stoi(fields[5]);
		if (fields.size() >= 7) {
			for (int i=6; i<fields.size(); i++)
				prediction.push_back(stof(fields[i]));
		}
		return 0;
	}

	return -1;

}

// Write to string in new format
//.......................................................................................
int MedSample::write_to_string(string &s)
{
	s = "";
	s += "SAMPLE\t" + to_string(id) + "\t" + to_string(time) + "\t" + to_string(outcome) + "\t" + to_string(outcomeTime);
	s += "\t" + to_string(split);
	for (auto p : prediction)
		s += "\t" + to_string(p);
	return 0;
}

// printing all samples with prefix appearing in the begining of each line
//.......................................................................................
void MedSample::print(const string prefix) {
	MLOG("%s :: id %d time %d outcomeTime %d outcome %f split %d prediction(%d)", prefix.c_str(), id, time, outcomeTime, outcome, split, prediction.size());
	if (prediction.size() > 0)
		for (auto pred : prediction)
			MLOG(" %f", pred);
	MLOG("\n");
}


//=======================================================================================
// MedSamples
//=======================================================================================
// Extract predictions from MedFeatures and insert to corresponding samples
// Samples in MedFeatures are assumed to be of the same size and order as in MedSamples
//.......................................................................................
int MedSamples::insert_preds(MedFeatures& features) {

	size_t size = 0;
	for (MedIdSamples& idSampels : idSamples)
		size += idSampels.samples.size();

	if (features.samples.size() != size) {
		MERR("Size mismatch between features and samples (%d vs %d)\n",features.samples.size(),size);
		return -1;
	}

	int idx = 0;
	for (MedIdSamples& idSample : idSamples) {
		for (unsigned int i = 0; i < idSample.samples.size(); i++) 
			idSample.samples[i].prediction = features.samples[idx++].prediction;
	}

	return 0;
}

// Get all patient ids
//.......................................................................................
void MedSamples::get_ids(vector<int>& ids) {

	ids.resize(idSamples.size());
	for (unsigned int i = 0; i < idSamples.size(); i++)
		ids[i] = idSamples[i].id;

}

// Extract a single vector of concatanated (vectors of) predictions
//.......................................................................................
void MedSamples::get_preds(vector<float>& preds) {
	for (auto& idSample : idSamples) 
		for (auto& sample : idSample.samples) 
			for (int i = 0; i < sample.prediction.size(); i++)
				preds.push_back(sample.prediction[i]);
}

// Extract a vector of all outcomes
//.......................................................................................
void MedSamples::get_y(vector<float>& y) {
	for (auto& idSample : idSamples)
		for (auto& sample : idSample.samples)
			y.push_back(sample.outcome);
}

// gets a list of all categories (different values) appearing in the outcome
//.......................................................................................
void MedSamples::get_categs(vector<float>& categs) 
{
	map<float, int> categ_inside;
	categs.clear();

	// Collect categories
	for (auto &id : idSamples)
		for (auto &rec : id.samples)
			categ_inside[rec.outcome] = 1;
	
	// Create a vector
	for (auto &it : categ_inside)
		categs.push_back(it.first);

}

// Helper function : get a vector of fields and generate the fields' positions map 
//.......................................................................................
int extract_field_pos_from_header(vector<string> field_names, map <string, int> & pos) {
	pos["id"] = -1;
	pos["date"] = -1;
	pos["outcome"] = -1;
	pos["outcome_date"] = -1;
	pos["pred"] = -1;
	pos["split"] = -1;

	for (int i = 0; i < field_names.size(); i++) {
		if (field_names[i] == "id" || field_names[i] == "pid")
			pos["id"] = i;
		else if (field_names[i] == "date" || field_names[i] == "time")
			pos["date"] = i;
		else if (field_names[i] == "outcome")
			pos["outcome"] = i;
		else if (field_names[i] == "outcomeTime" || field_names[i] == "outcome_date")
			pos["outcome_date"] = i;
		else if (field_names[i] == "prediction" || field_names[i] == "pred")
			pos["pred"] = i;
		else if (field_names[i] == "split")
			pos["split"] = i;
		else MWARN("WARNING: header line contains [%s] which is not part of MedSample\n",field_names[i].c_str());
	}
	for (auto& e : pos)
		if (e.second == -1)
			MWARN("WARNING: header line does not contain [%s]\n", e.first.c_str());
		else MLOG("header line contains [%s] at column [%d]\n", e.first.c_str(), e.second);

	return 0;
}

// read from text file.
// If the line starting with EVENT_FIELDS (followed by tabe-delimeted field names : id,date,outcome,outcome_date,split,pred) appears before the data lines, it is used to determine
// fields positions, otherwise - old or new formats are used.
//-------------------------------------------------------------------------------------------
int MedSamples::read_from_file(const string &fname)
{
	ifstream inf(fname);

	MLOG("MedSamples: reading %s\n", fname.c_str());
	if (!inf) {
		MERR("MedSamples: can't open file %s for read\n", fname.c_str());
		return -1;
	}

	string curr_line;

	int samples = 0, read_records = 0, skipped_records = 0;
	idSamples.clear();
	int curr_id = -1;
	unordered_set<int> seen_ids;
	map<string, int> pos;
	while (getline(inf, curr_line)) {
		//MLOG("--> %s\n",curr_line.c_str());
		if ((curr_line.size() > 1) && (curr_line[0] != '#')) {
			if (curr_line[curr_line.size() - 1] == '\r')
				curr_line.erase(curr_line.size() - 1);
			read_records++;
			vector<string> fields;
			split(fields, curr_line, boost::is_any_of("\t"));
			if (fields.size() >= 2) {
				if (fields[0] == "NAME") MLOG("reading NAME = %s\n", fields[1].c_str());
				else if (fields[0] == "DESC") MLOG("reading DESC = %s\n", fields[1].c_str());
				else if (fields[0] == "TYPE") MLOG("reading TYPE = %s\n", fields[1].c_str());
				else if (fields[0] == "NCATEG")  MLOG("reading NCATEG = %s\n", fields[1].c_str());
				else if ((fields[0] == "EVENT_FIELDS" || fields[0] == "pid" || fields[0] == "id") && read_records == 1) {
					extract_field_pos_from_header(fields, pos);
					continue;
				}
				MedSample sample;
				 
				if (sample.parse_from_string(curr_line, pos) < 0) {
					MWARN("skipping [%s]\n", curr_line.c_str());
					skipped_records++;
					if (read_records > 30 && skipped_records > read_records / 2)
						MTHROW_AND_ERR("skipped %d/%d first records, exiting\n", skipped_records, read_records);
					continue;
				}
				if (sample.id != curr_id) {
					if (seen_ids.find(sample.id) != seen_ids.end()) {
						MERR("ERROR: Wrong MedSample format in line \"%s\"", curr_line.c_str());
						MTHROW_AND_ERR("Sample id [%d] records are not consecutive\n", sample.id);
					}
					seen_ids.insert(sample.id);
					// new idSample
					MedIdSamples mis;
					mis.id = sample.id;
					mis.split = sample.split;
					mis.samples.push_back(sample);
					curr_id = sample.id;
					idSamples.push_back(mis);
				}
				else if (sample.id == curr_id) {
					// another sample for the current MedIdSamples
					if (idSamples.back().id != sample.id || idSamples.back().split != sample.split) {
						MERR("ERROR: Wrong MedSample format in line \"%s\"", curr_line.c_str());
						MERR("Got conflicting split : %d,%d vs. %d,%d\n", idSamples.back().id, idSamples.back().split, sample.id, sample.split);
						return -1;
					}
					idSamples.back().samples.push_back(sample);
				}
				samples++;
			}
		}
	}
	MLOG("read [%d] samples for [%d] patient IDs. Skipped [%d] records\n", samples, idSamples.size(), skipped_records);
	sort_by_id_date();
	inf.close();
	return 0;
}

// Sort by id and then date
//.......................................................................................
void MedSamples::sort_by_id_date() {
	MLOG("sorting samples by id, date\n");
	sort(idSamples.begin(), idSamples.end(), comp_patient_id_time);
	for (auto& pat : idSamples)
		sort(pat.samples.begin(), pat.samples.end(), comp_sample_id_time);
}

// Make sure that : (1) every pid has one idSample at most and (2) everything is sorted
//.......................................................................................
void MedSamples::normalize() {
	
	// since order may be random, we need a map to collect by pid
	map<int, vector<MedSample>> pid_to_samples;
	map<int, int> pid_to_split;
	for (auto &ids : idSamples) {
		pid_to_split[ids.id] = ids.split;
		for (auto &s : ids.samples)
			pid_to_samples[s.id].push_back(s);
	}

	// copy back to idSamples and sort
	idSamples.clear();
	for (auto &vs : pid_to_samples) {
		MedIdSamples ids;
		ids.id = vs.first;
		ids.split = pid_to_split[ids.id];
		ids.samples = vs.second;
		sort(ids.samples.begin(), ids.samples.end(), comp_sample_id_time);
		idSamples.push_back(ids);
	}
}

// write to text file in new format
//.......................................................................................
int MedSamples::write_to_file(const string &fname)
{
	ofstream of(fname);

	MLOG("MedSamples: writing to %s\n", fname.c_str());
	if (!of) {
		MERR("MedSamples: can't open file %s for writing\n", fname.c_str());
		return -1;
	}
	int samples = 0;
	int buffer_write = 100000;
		
	//of << "EVENT_FIELDS" << '\t' << "id" << '\t' << "time" << '\t' << "outcome" << '\t' << "outcomeLength" <<
	//	'\t' << "outcomeTime" << '\t' << "split" << '\t' << "prediction" << endl;
	of << "EVENT_FIELDS" << '\t' << "id" << '\t' << "time" << '\t' << "outcome" << '\t' << "outcomeTime" << '\t' << "split" << '\t' << "prediction" << endl;

	int line = 0;
	for (auto &s: idSamples) {
		for (auto ss : s.samples) {
			samples++;
			string sout;
			ss.write_to_string(sout);
			//of << "EVENT" << '\t' << ss.id << '\t' << ss.time << '\t' << ss.outcome << '\t' << 100000 << '\t' <<
			//	ss.outcomeTime << '\t' << s.split << '\t' << ss.prediction.front() << endl;
			if (line >= buffer_write) {
				of << sout << endl;
				line = 0;
			}
			else {
				of << sout << "\n"; //no flush of buffer - much faster when writing large files
				++line;
			}
		}
	}

	MLOG("wrote [%d] samples for [%d] patient IDs\n", samples, idSamples.size());
	of.close(); //will flush buffer if needed
	return 0;
}

// Count samples
//.......................................................................................
int MedSamples::nSamples()
{
	int n = 0;
	for (auto& idSample : idSamples)
		n += (int)idSample.samples.size();

	return n;
}

// API's for online insertions : main use case is a single time point for prediction per pid
//.......................................................................................
int MedSamples::insertRec(int pid, int time, float outcome, int outcomeTime) 
{
	MedIdSamples sample;

	sample.id = pid;
	sample.split = -1;
	MedSample s;
	s.id = pid;
	s.time = time;
	s.outcome = outcome;
	s.outcomeTime = outcomeTime;
	sample.samples.push_back(s);
	idSamples.push_back(sample);
	return 0;
}

//.......................................................................................
int MedSamples::insertRec(int pid, int time, float outcome, int outcomeTime, float pred)
{
	MedIdSamples sample;

	sample.id = pid;
	sample.split = -1;
	MedSample s;
	s.id = pid;
	s.time = time;
	s.outcome = outcome;
	s.outcomeTime = outcomeTime;
	s.prediction.push_back(pred);
	sample.samples.push_back(s);
	idSamples.push_back(sample);
	return 0;
}

// Get all MedSamples as a single vector
//.......................................................................................
void MedSamples::export_to_sample_vec(vector<MedSample> &vec_samples)
{
	vec_samples.clear();
	for (auto &s: idSamples) {
		for (auto &samp : s.samples) {
			vec_samples.push_back(samp);
		}
	}
}

//.......................................................................................
void MedSamples::dilute(float prob)
{
	if (prob >= 1)
		return;

	vector<MedIdSamples> NewidSamples;

	for (auto &id : idSamples) {
		MedIdSamples mid;
		mid.id = id.id;
		mid.split = id.split;
		for (auto &s : id.samples)
			if (rand_1() < prob)
				mid.samples.push_back(s);
		if (mid.samples.size() > 0)
			NewidSamples.push_back(mid);
	}

	idSamples = NewidSamples;
}
