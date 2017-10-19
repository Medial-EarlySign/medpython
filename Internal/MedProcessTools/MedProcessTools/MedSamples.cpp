#include "MedSamples.h"
#include "MedProcessTools/MedProcessTools/MedFeatures.h"
#include "Logger/Logger/Logger.h"

#define LOCAL_SECTION MED_SAMPLES_CV
#define LOCAL_LEVEL	LOG_DEF_LEVEL

//=======================================================================================
// MedSample
//=======================================================================================

int MedSample::parse_from_string(string &s, map <string, int> & pos) {
	if (pos.size() == 0)
		return parse_from_string(s);
	vector<string> fields; 
	if (fields.size() == 0)
		return -1;
	boost::split(fields, s, boost::is_any_of("\t\n\r"));
	if (pos["id"] != -1)
		id = stoi(fields[pos["id"]]);
	if (pos["date"] != -1)
		time = stoi(fields[pos["time"]]);
	if (pos["outcome"] != -1)
		outcome = stof(fields[pos["outcome"]]);
	if (pos["outcome_date"] != -1)
		outcomeTime = stoi(fields[pos["outcome_date"]]);
	if (pos["split"] != -1)
		split = stoi(fields[pos["split"]]);
	if (pos["pred"] != -1)
		prediction.push_back(stof(fields[pos["pred"]]));
	return 0;
}

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

//.......................................................................................
void MedSample::print(const string prefix) {
	MLOG("%s :: id %d time %d outcomeTime %d outcome %f split %d prediction(%d)", prefix.c_str(), id, time, outcomeTime, outcome, split, prediction.size());
	if (prediction.size() > 0)
		for (auto pred : prediction)
			MLOG(" %f", pred);
	MLOG("\n");
}

//=======================================================================================
// MedIdSample
//=======================================================================================


//=======================================================================================
// MedSamples
//=======================================================================================
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

//.......................................................................................
void MedSamples::get_ids(vector<int>& ids) {

	ids.resize(idSamples.size());
	for (unsigned int i = 0; i < idSamples.size(); i++)
		ids[i] = idSamples[i].id;

}

//.......................................................................................
void MedSamples::get_preds(vector<float>& preds) {
	for (auto& idSample : idSamples) 
		for (auto& sample : idSample.samples) 
			for (int i = 0; i < sample.prediction.size(); i++)
				preds.push_back(sample.prediction[i]);
}

//.......................................................................................
void MedSamples::get_y(vector<float>& y) {
	for (auto& idSample : idSamples)
		for (auto& sample : idSample.samples)
			y.push_back(sample.outcome);
}

//.......................................................................................
void MedSamples::get_categs(vector<float>& categs) 
{
	map<float, int> categ_inside;
	categs.clear();

	for (auto &id : idSamples)
		for (auto &rec : id.samples)
			categ_inside[rec.outcome] = 1;
	
	for (auto &it : categ_inside)
		categs.push_back(it.first);

}

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

	int samples = 0, skipped_records = 0;
	idSamples.clear();
	int curr_id = -1;
	unordered_set<int> seen_ids;
	while (getline(inf, curr_line)) {
		//MLOG("--> %s\n",curr_line.c_str());
		if ((curr_line.size() > 1) && (curr_line[0] != '#')) {
			if (curr_line[curr_line.size() - 1] == '\r')
				curr_line.erase(curr_line.size() - 1);

			vector<string> fields;
			split(fields, curr_line, boost::is_any_of("\t"));
			map<string, int> pos;
			if (fields.size() >= 2) {

				if (fields[0] == "NAME") MLOG("reading NAME = %s\n", fields[1].c_str());
				if (fields[0] == "DESC") MLOG("reading DESC = %s\n", fields[1].c_str());
				if (fields[0] == "TYPE") MLOG("reading TYPE = %s\n", fields[1].c_str());
				if (fields[0] == "NCATEG")  MLOG("reading NCATEG = %s\n", fields[1].c_str());
				if (fields[0] == "EVENT_FIELDS") {
					extract_field_pos_from_header(fields, pos);
					continue;
				}
				MedSample sample;
				 
				if (sample.parse_from_string(curr_line, pos) < 0) {
					MWARN("skipping [%s]\n", curr_line.c_str());
					skipped_records++;
					continue;
				}
				if (sample.id != curr_id) {
					if (seen_ids.find(sample.id) != seen_ids.end())
						MTHROW_AND_ERR("Sample id [%d] records are not consecutive\n", sample.id);
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

void MedSamples::sort_by_id_date() {
	MLOG("sorting samples by id, date\n");
	sort(idSamples.begin(), idSamples.end(), comp_patient_id_time);
	for (auto& pat : idSamples)
		sort(pat.samples.begin(), pat.samples.end(), comp_sample_id_time);
}

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

	//of << "EVENT_FIELDS" << '\t' << "id" << '\t' << "time" << '\t' << "outcome" << '\t' << "outcomeLength" <<
	//	'\t' << "outcomeTime" << '\t' << "split" << '\t' << "prediction" << endl;
	of << "EVENT_FIELDS" << '\t' << "id" << '\t' << "time" << '\t' << "outcome" << '\t' << "outcomeTime" << '\t' << "split" << '\t' << "prediction" << endl;

	for (auto &s: idSamples) {
		for (auto ss : s.samples) {
			samples++;
			string sout;
			ss.write_to_string(sout);
			//of << "EVENT" << '\t' << ss.id << '\t' << ss.time << '\t' << ss.outcome << '\t' << 100000 << '\t' <<
			//	ss.outcomeTime << '\t' << s.split << '\t' << ss.prediction.front() << endl;
			of << sout << endl;
		}
	}

	MLOG("wrote [%d] samples for [%d] patient IDs\n", samples, idSamples.size());
	of.close();
	return 0;
}

// Count
//.......................................................................................
int MedSamples::nSamples()
{
	int n = 0;
	for (auto& idSample : idSamples)
		n += (int)idSample.samples.size();

	return n;
}

//.......................................................................................
int MedSamples::insertRec(int pid, int time, float outcome, int outcomeTime) 
{
	MedIdSamples sample;

	sample.id = pid;
	MedSample s;
	s.id = pid;
	s.time = time;
	s.outcome = outcome;
	s.outcomeTime = outcomeTime;
	sample.samples.push_back(s);
	idSamples.push_back(sample);
	return 0;
}