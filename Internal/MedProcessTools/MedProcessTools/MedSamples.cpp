#include "MedSamples.h"
#include "MedProcessTools/MedProcessTools/MedFeatures.h"
#include "Logger/Logger/Logger.h"
#include <MedUtils/MedUtils/MedGenUtils.h>
#include <MedUtils/MedUtils/MedUtils.h>
#include <boost/crc.hpp>
#include <random>
#include <algorithm>

#define LOCAL_SECTION MED_SAMPLES_CV
#define LOCAL_LEVEL	LOG_DEF_LEVEL

//=======================================================================================
// MedSample
//=======================================================================================
// Get sample from tab-delimited string, where pos indicate the position of each field (fields are id,date,outcome,outcome_date,split) in addition to pred_pos vector and attr_pos map
//.......................................................................................
int MedSample::parse_from_string(string &s, map <string, int> & pos, vector<int>& pred_pos, map<string, int>& attr_pos, map<string, int>& str_attr_pos, int time_unit, int raw_format) {
	if (pos.size() == 0)
		return parse_from_string(s, time_unit);
	vector<string> fields;
	boost::split(fields, s, boost::is_any_of("\t\n\r"));
	if (fields.size() == 0)
		return -1;
	try {
		if (pos["id"] != -1)
			id = (int)stod(fields[pos["id"]]);
		if (pos["date"] != -1) {
			if (raw_format)
				time = stoi(fields[pos["date"]]);
			else
				time = med_time_converter.convert_datetime_safe(time_unit, fields[pos["date"]], 2);
		}
		if (pos["outcome"] != -1)
			outcome = stof(fields[pos["outcome"]]);
		if (pos["outcome_date"] != -1) {
			if (raw_format)
				outcomeTime = stoi(fields[pos["outcome_date"]]);
			else
				outcomeTime = med_time_converter.convert_datetime_safe(time_unit, fields[pos["outcome_date"]], 1);
		}
		if (pos["split"] != -1 && fields.size() > pos["split"])
			split = stoi(fields[pos["split"]]);

		for (int pos : pred_pos) {
			if (pos != -1 && fields.size() > pos)
				prediction.push_back(stof(fields[pos]));
		}

		for (auto& attr : attr_pos) {
			if (attr.second != -1 && fields.size() > attr.second)
				attributes[attr.first] = stof(fields[attr.second]);
		}

		for (auto& attr : str_attr_pos) {
			if (attr.second != -1 && fields.size() > attr.second)
				str_attributes[attr.first] = fields[attr.second];
		}

		return 0;
	}
	catch (std::invalid_argument e) {
		MLOG("could not parse [%s]\n", s.c_str());
		return -1;
	}
}

// Get sample from tab-delimited string, in old or new format (<split> and <prediction> optional, <predictions> can be several numbers (tab delimited))
//.......................................................................................
int MedSample::parse_from_string(string &s, int time_unit)
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
		time = med_time_converter.convert_datetime_safe(time_unit, fields[2], 2);
		outcome = stof(fields[3]);
		//outcomeTime = stoi(fields[5]);
		outcomeTime = med_time_converter.convert_datetime_safe(time_unit, fields[5], 2);
		if (fields.size() >= 7)
			split = stoi(fields[6]);
		if (fields.size() >= 8) {
			for (int i = 7; i < fields.size(); i++)
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
		time = med_time_converter.convert_datetime_safe(time_unit, fields[2], 2);
		outcome = stof(fields[3]);
		//outcomeTime = stoi(fields[4]);
		outcomeTime = med_time_converter.convert_datetime_safe(time_unit, fields[4], 2);
		if (fields.size() >= 6)
			split = stoi(fields[5]);
		if (fields.size() >= 7) {
			for (int i = 6; i < fields.size(); i++)
				prediction.push_back(stof(fields[i]));
		}
		return 0;
	}

	return -1;

}


void MedSample::write_to_string(string &s, int time_unit)
{
	vector<string> my_attributes, my_str_attributes;
	for (auto& attr : attributes)
		my_attributes.push_back(attr.first);
	for (auto& attr : str_attributes)
		my_str_attributes.push_back(attr.first);
	write_to_string(s, my_attributes, my_str_attributes, time_unit);
}

// Write to string in new format
//.......................................................................................
void MedSample::write_to_string(string &s, const vector<string>& attr, const vector<string>& str_attr, int time_unit)
{
	s = "";
	s += "SAMPLE\t" + to_string(id) + "\t" + med_time_converter.convert_times_S(time_unit, MedTime::DateTimeString, time) + "\t" + to_string(outcome) + "\t"
		+ med_time_converter.convert_times_S(time_unit, MedTime::DateTimeString, outcomeTime);
	s += "\t" + to_string(split);
	for (auto p : prediction)
		s += "\t" + to_string(p);
	for (const string& a : attr)
		s += "\t" + to_string(attributes[a]);
	for (const string& a : str_attr)
		s += "\t" + str_attributes[a];
	return;
}

// printing all samples with prefix appearing in the begining of each line
//.......................................................................................
void MedSample::print(const string prefix) {
	MLOG("%s :: id %d time %d outcomeTime %d outcome %f split %d prediction(%d)", prefix.c_str(), id, time, outcomeTime, outcome, split, prediction.size());
	if (prediction.size() > 0)
		for (auto pred : prediction)
			MLOG(" %f", pred);
	for (auto& attr : attributes)
		MLOG("%s=%f", attr.first.c_str(), attr.second);
	for (auto& attr : str_attributes)
		MLOG("%s=%s", attr.first.c_str(), attr.second.c_str());
	MLOG("\n");
}

//=======================================================================================
// MedIdSamples
//=======================================================================================
// Comparison function : mode 0 requires equal id/time, mode 1 requires equal outcome info, mode 2 also compares split, attributes and prediction
//.......................................................................................
bool MedIdSamples::same_as(MedIdSamples &other, int mode) {
	if (other.samples.size() != samples.size())
		return false;

	for (unsigned int i = 0; i < samples.size(); i++) {

		if (samples[i].id != other.samples[i].id || samples[i].time != other.samples[i].time)
			return false;

		if (mode > 0 && (samples[i].outcome != other.samples[i].outcome && samples[i].outcomeTime != other.samples[i].outcomeTime))
			return false;

		if (mode > 1) {
			if (samples[i].split != other.samples[i].split || samples[i].prediction.size() != other.samples[i].prediction.size() || samples[i].attributes.size() != other.samples[i].attributes.size())
				return false;
			for (unsigned int j = 0; j < samples[i].prediction.size(); j++) {
				if (samples[i].prediction[j] != other.samples[i].prediction[j])
					return false;
			}
			for (auto& attr : samples[i].attributes) {
				if (other.samples[i].attributes.find(attr.first) == other.samples[i].attributes.end() || other.samples[i].attributes[attr.first] != attr.second)
					return false;
			}

			for (auto& attr : samples[i].str_attributes) {
				if (other.samples[i].str_attributes.find(attr.first) == other.samples[i].str_attributes.end() || other.samples[i].str_attributes[attr.first] != attr.second)
					return false;
			}
		}
	}

	return true;
}

//=======================================================================================
// MedSamples
//=======================================================================================
// Extract predictions from MedFeatures and insert to corresponding samples
// Samples in MedFeatures are assumed to be of the same size and order as in MedSamples
// Return -1 if samples and features do not match in length, 0 upon success
//.......................................................................................
int MedSamples::insert_preds(MedFeatures& features) {

	size_t size = (size_t)nSamples();
	if (features.samples.size() != size) {
		MERR("Size mismatch between features and samples (%d vs %d)\n", features.samples.size(), size);
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
int extract_field_pos_from_header(vector<string> field_names, map <string, int> & pos, vector<int>& pred_pos, map<string, int>& attr_pos, map<string, int>& str_attr_pos) {
	pos["id"] = -1;
	pos["date"] = -1;
	pos["outcome"] = -1;
	pos["outcome_date"] = -1;
	pos["split"] = -1;

	vector<string> unknown_fields;
	for (int i = 0; i < field_names.size(); i++) {
		if (field_names[i] == "id" || field_names[i] == "pid")
			pos["id"] = i;
		else if (field_names[i] == "date" || field_names[i] == "time")
			pos["date"] = i;
		else if (field_names[i] == "outcome")
			pos["outcome"] = i;
		else if (field_names[i] == "outcomeTime" || field_names[i] == "outcome_date")
			pos["outcome_date"] = i;
		else if (field_names[i] == "split")
			pos["split"] = i;
		else if (field_names[i] == "prediction" || field_names[i] == "pred" || field_names[i].substr(0, 5) == "pred_") // Note that we don't check that pred_# are actually ordered
			pred_pos.push_back(i);
		else if (field_names[i].substr(0, 5) == "attr_") {
			string attr_name = field_names[i].substr(5, field_names[i].length() - 5);
			attr_pos[attr_name] = i;
		}
		else if (field_names[i].substr(0, 5) == "str_attr_") {
			string attr_name = field_names[i].substr(9, field_names[i].length() - 9);
			str_attr_pos[attr_name] = i;
		}
		else unknown_fields.push_back(field_names[i]);
	}
	if (unknown_fields.size() > 0) {
		string warning = "WARNING: header line contains unused fields [";
		for (string u : unknown_fields)
			warning += u + ",";
		warning += "]";
		MWARN("%s\n", warning.c_str());
	}
	for (auto& e : pos)
		if (e.second == -1)
			MWARN("[%s]=unspecified, ", e.first.c_str());
		else MLOG("[%s]=%d, ", e.first.c_str(), e.second);
		MLOG("\n");
		return 0;
}

// read from text file.
// If the line starting with EVENT_FIELDS (followed by tabe-delimeted field names : id,date,outcome,outcome_date,split,preds,attr) appears before the data lines, it is used to determine
// fields positions, otherwise - old or new formats are used. Return -1 upon failure to open file
//-------------------------------------------------------------------------------------------
int MedSamples::read_from_file(const string &fname, bool sort_rows)
{

	ifstream inf(fname);

	MLOG("MedSamples: reading %s\n", fname.c_str());
	if (!inf) {
		MTHROW_AND_ERR("MedSamples: can't open file %s for read\n", fname.c_str());
	}

	string curr_line;

	int samples = 0, read_records = 0, skipped_records = 0;
	idSamples.clear();
	int curr_id = -1;
	unordered_set<int> seen_ids;
	map<string, int> pos;
	vector<int> pred_pos;
	map<string, int> attr_pos;
	map<string, int> str_attr_pos;
	time_unit = global_default_time_unit;

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
				else if (fields[0] == "RAW_FORMAT") {
					MLOG("reading RAW_FORMAT = %s\n", fields[1].c_str());
					raw_format = stoi(fields[1]);
					MLOG("raw_format is %d\n", raw_format);
				}
				else if (fields[0] == "TIME_UNIT") {
					MLOG("reading TIME_UNIT = %s\n", fields[1].c_str());
					time_unit = med_time_converter.string_to_type(fields[1]);
					MLOG("time unit is %d\n", time_unit);
				}
				else if ((fields[0] == "EVENT_FIELDS" || fields[0] == "pid" || fields[0] == "id") && read_records == 1) {
					extract_field_pos_from_header(fields, pos, pred_pos, attr_pos, str_attr_pos);
					continue;
				}

				else {
					MedSample sample;

					if (sample.parse_from_string(curr_line, pos, pred_pos, attr_pos, str_attr_pos, time_unit, raw_format) < 0) {
						MWARN("skipping [%s]\n", curr_line.c_str());
						skipped_records++;
						if (read_records > 30 && skipped_records > read_records / 2)
							MTHROW_AND_ERR("skipped %d/%d first records, exiting\n", skipped_records, read_records);
						continue;
					}

					if (0 && read_records < 10) { // use for debug
						MLOG("Read samples line %s ... parsed it as : pid %d time %d outcome %f outcomeTime %d\n", curr_line.c_str(), sample.id, sample.time, sample.outcome, sample.outcomeTime);
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
	}
	MLOG("read [%d] samples for [%d] patient IDs. Skipped [%d] records\n", samples, idSamples.size(), skipped_records);
	if (sort_rows)
		sort_by_id_date();
	inf.close();
	return 0;
}

// Get predictions vector size. Return -1 if not-consistent
//-------------------------------------------------------------------------------------------
int MedSamples::get_predictions_size(int& nPreds) {

	nPreds = -1;
	for (auto &s : idSamples) {
		for (auto& ss : s.samples) {
			if (nPreds == -1)
				nPreds = (int)ss.prediction.size();
			else if (ss.prediction.size() != nPreds)
				return -1;
		}
	}
	return 0;

}


// Get all attributes 
//-------------------------------------------------------------------------------------------
int MedSamples::get_all_attributes(vector<string>& attributes, vector<string>& str_attributes) {
	attributes.clear();
	str_attributes.clear();
	bool first = true;
	for (auto &s : idSamples) {
		for (auto& ss : s.samples) {
			if (first) {
				for (auto& attr : ss.attributes)
					attributes.push_back(attr.first);
				for (auto& attr : ss.str_attributes)
					str_attributes.push_back(attr.first);
				first = false;
			}
			else {
				vector<string> my_attributes;
				for (auto& attr : ss.attributes)
					my_attributes.push_back(attr.first);
				if (attributes != my_attributes)
					MTHROW_AND_ERR("attributes are not the same across all samples");

				my_attributes.clear();
				for (auto& attr : ss.str_attributes)
					my_attributes.push_back(attr.first);
				if (str_attributes != my_attributes)
					MTHROW_AND_ERR("str_attributes are not the same across all samples");
			}
		}
	}
	return 0;
}

// write to text file in new format. 
// Return -1 upon failure to open file.
// Return -2 upon prediction-length inconsistency
// Return -3 upon attributes inconsistency
//.......................................................................................
int MedSamples::write_to_file(const string &fname)
{
	ofstream of(fname);

	int nPreds;
	if (get_predictions_size(nPreds) < 0) {
		MERR("MedSampels: Predictions vectors sizes inconsistent\n");
		return -2;
	}
	vector<string> attributes, str_attributes;
	if (get_all_attributes(attributes, str_attributes) < 0) {
		MERR("MedSamples: Attributes sets inconsistency\n");
		return -3;
	}
	MLOG("attributes size %d\n", attributes.size());

	MLOG("MedSamples: writing to %s\n", fname.c_str());
	if (!of) {
		MERR("MedSamples: can't open file %s for writing\n", fname.c_str());
		return -1;
	}
	int samples = 0;
	int buffer_write = 0;

	of << "EVENT_FIELDS" << '\t' << "id" << '\t' << "time" << '\t' << "outcome" << '\t' << "outcomeTime" << '\t' << "split";
	for (int i = 0; i < nPreds; i++)
		of << "\tpred_" << i;
	for (string name : attributes)
		of << "\tattr_" << name;
	for (string name : str_attributes)
		of << "\tstr_attr_" << name;
	of << "\n";

	//removing this line since it creates something which isn't a csv, and I can't read it with python later...
	//of << "TIME_UNIT" << "\t" << med_time_converter.type_to_string(time_unit) << "\n";

	int line = 0;
	for (auto &s : idSamples) {
		for (auto ss : s.samples) {
			samples++;
			string sout;
			ss.write_to_string(sout, attributes, str_attributes, time_unit);
			//of << "EVENT" << '\t' << ss.id << '\t' << ss.time << '\t' << ss.outcome << '\t' << 100000 << '\t' <<
			//	ss.outcomeTime << '\t' << s.split << '\t' << ss.prediction.front() << endl;
			if (buffer_write > 0 && line >= buffer_write) {
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

// Count samples
//.......................................................................................
int MedSamples::nSamples() const
{
	int n = 0;
	for (auto& idSample : idSamples)
		n += (int)idSample.samples.size();

	return n;
}

// API's for online insertions : main use case is a single time point for prediction per pid
//.......................................................................................
void MedSamples::insertRec(int pid, int time, float outcome, int outcomeTime)
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
	return;
}

//.......................................................................................
void MedSamples::insertRec(int pid, int time, float outcome, int outcomeTime, float pred)
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
	return;
}

// Get all MedSamples as a single vector
//.......................................................................................
void MedSamples::export_to_sample_vec(vector<MedSample> &vec_samples)
{
	vec_samples.clear();
	for (auto &s : idSamples) {
		for (auto &samp : s.samples) {
			vec_samples.push_back(samp);
		}
	}
}

// Create a MedSamples object from a vector of MedSample
//.......................................................................................
void MedSamples::import_from_sample_vec(vector<MedSample> &vec_samples, bool allow_split_inconsistency) {

	idSamples.clear();
	map<int, int> id2idx;
	map<int, int> id2split;

	for (MedSample& sample : vec_samples) {
		if (id2idx.find(sample.id) == id2idx.end()) {
			id2idx[sample.id] = (int)idSamples.size();

			idSamples.resize(idSamples.size() + 1);
			idSamples.back().id = sample.id;
			idSamples.back().split = sample.split;
		}

		int idx = id2idx[sample.id];
		if (!allow_split_inconsistency && idSamples[idx].split != sample.split)
			MTHROW_AND_ERR("Split incosistency for pid=%d\n", sample.id);
		idSamples[idx].samples.push_back(sample);
	}

	// Sort
	sort_by_id_date();
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

// Comparison function : mode 0 requires equal id/time, mode 1 requires equal outcome info, mode 2 also compares split and prediction
//.......................................................................................
bool MedSamples::same_as(MedSamples &other, int mode) {
	if (other.idSamples.size() != idSamples.size())
		return false;

	for (unsigned int i = 0; i < idSamples.size(); i++) {
		if (!idSamples[i].same_as(other.idSamples[i], mode))
			return false;
	}

	return true;
}

void MedSamples::flatten(vector<MedSample> &flat) const {
	for (size_t i = 0; i < idSamples.size(); ++i)
		for (size_t j = 0; j < idSamples[i].samples.size(); ++j)
			flat.push_back(idSamples[i].samples[j]);
}

void medial::print::print_samples_stats(const vector<MedSample> &samples, const string &log_file) {
	ofstream fo;
	if (!log_file.empty()) {
		fo.open(log_file, ios::app);
		if (!fo.good())
			MWARN("Warning: can log into file %s\n", log_file.c_str());
	}

	map<float, int> histCounts, histCountAll;
	vector<unordered_set<int>> pid_index(2);
	for (size_t k = 0; k < samples.size(); ++k)
	{
		if (pid_index[samples[k].outcome > 0].find(samples[k].id) == pid_index[samples[k].outcome > 0].end()) {
			if (histCounts.find(samples[k].outcome) == histCounts.end()) {
				histCounts[samples[k].outcome] = 0;
			}
			++histCounts[samples[k].outcome];
		}
		if (histCountAll.find(samples[k].outcome) == histCountAll.end()) {
			histCountAll[samples[k].outcome] = 0;
		}
		++histCountAll[samples[k].outcome];
		pid_index[samples[k].outcome > 0].insert(samples[k].id);
	}
	int total = 0, total_all = 0;
	for (auto it = histCounts.begin(); it != histCounts.end(); ++it)
		total += it->second;
	for (auto it = histCountAll.begin(); it != histCountAll.end(); ++it)
		total_all += it->second;

	log_with_file(fo, "Samples has %d records. for uniq_pids = [", (int)samples.size());
	auto iter = histCounts.begin();
	if (!histCounts.empty()) {
		log_with_file(fo, "%d=%d(%2.2f%%)", (int)iter->first, iter->second,
			100.0 * iter->second / float(total));
		++iter;
	}
	for (; iter != histCounts.end(); ++iter)
		log_with_file(fo, ", %d=%d(%2.2f%%)", (int)iter->first, iter->second,
			100.0 * iter->second / float(total));

	log_with_file(fo, "] All = [");
	iter = histCountAll.begin();
	if (!histCountAll.empty()) {
		log_with_file(fo, "%d=%d(%2.2f%%)", (int)iter->first, iter->second,
			100.0 * iter->second / float(total_all));
		++iter;
	}

	for (; iter != histCountAll.end(); ++iter)
		log_with_file(fo, ", %d=%d(%2.2f%%)", (int)iter->first, iter->second, iter->second,
			100.0 * iter->second / float(total_all));
	log_with_file(fo, "]\n");
	if (fo.good())
		fo.close();
}

void medial::print::print_samples_stats(const MedSamples &samples, const string &log_file) {
	vector<MedSample> smps;
	for (size_t i = 0; i < samples.idSamples.size(); ++i)
		smps.insert(smps.end(), samples.idSamples[i].samples.begin(), samples.idSamples[i].samples.end());
	print_samples_stats(smps, log_file);
}

void medial::print::print_by(const vector<MedSample> &data_records, const vector<string> &groups,
	bool unique_ids, const string &log_file) {
	ofstream fo;
	if (data_records.size() != groups.size())
		MTHROW_AND_ERR("Error in medial::print::print_b - data_records and groups should be same size\n");
	if (!log_file.empty()) {
		fo.open(log_file, ios::app);
		if (!fo.good())
			MWARN("Warning: can log into file %s\n", log_file.c_str());
	}

	unordered_map<string, int> count_0, count_1;
	vector<string> all_groups;
	unordered_set<string> seen_group;
	vector<unordered_map<string, unordered_set<int>>> year_to_seen_pid(2); //of_predicition
	for (size_t i = 0; i < data_records.size(); ++i)
	{
		int label = int(data_records[i].outcome > 0);

		if ((label > 0) && seen_group.find(groups[i]) == seen_group.end()) {
			all_groups.push_back(groups[i]);
			seen_group.insert(groups[i]);
		}

		if (label > 0) {
			if (!unique_ids ||
				year_to_seen_pid[1][groups[i]].find(data_records[i].id) == year_to_seen_pid[1][groups[i]].end()) {
				++count_1[groups[i]];
				year_to_seen_pid[1][groups[i]].insert(data_records[i].id);
			}
		}
		else {
			if (!unique_ids ||
				year_to_seen_pid[0][groups[i]].find(data_records[i].id) == year_to_seen_pid[0][groups[i]].end()) {
				++count_0[groups[i]];
				year_to_seen_pid[0][groups[i]].insert(data_records[i].id);
			}
		}
	}

	unordered_map<string, int> group_total;
	unordered_map<string, float>  group_ratio;
	int i = 0;
	sort(all_groups.begin(), all_groups.end());
	log_with_file(fo, "Group"  "\t"  "Controls"  "\t"  "Cases"  "\t"  "outcome_percentage\n");
	for (const string &grp : all_groups)
	{
		group_total[grp] = count_0[grp] + count_1[grp];
		group_ratio[grp] = count_1[grp] / float(count_0[grp] + count_1[grp]);
		++i;
		log_with_file(fo, "%s\t%d\t%d\t%2.3f%%\n", grp.c_str(), count_0[grp], count_1[grp],
			100 * count_1[grp] / float(count_1[grp] + count_0[grp]));
	}
	//special case for binary: show work point and AUC:
	if (all_groups.size() == 2) {
		float fpr = count_0[all_groups[1]] / float(count_0[all_groups[0]] + count_0[all_groups[1]]);
		float tpr = count_1[all_groups[1]] / float(count_1[all_groups[0]] + count_1[all_groups[1]]);
		float prior = (count_1[all_groups[0]] + count_1[all_groups[1]]) / float(data_records.size());
		float auc = 0.5*(fpr*tpr + (1 - fpr)*(1 + tpr));
		float lift = 0;
		if (group_ratio[all_groups[0]] > 0)
			lift = group_ratio[all_groups[1]] / group_ratio[all_groups[0]];
		log_with_file(fo, "FPR=%2.4f%%, TPR=%2.4f%%, AUC=%2.4f, prior=%2.4f%%, lift_between=%2.4f\n",
			100.0*fpr, 100.0*tpr, auc, 100.0*prior, lift);
	}

	if (fo.good())
		fo.close();
}

void medial::print::print_by(const MedSamples &data_records, const vector<string> &groups,
	bool unique_ids, const string &log_file) {
	vector<MedSample> vec;
	for (size_t i = 0; i < data_records.idSamples.size(); ++i)
		vec.insert(vec.end(), data_records.idSamples[i].samples.begin(), data_records.idSamples[i].samples.end());
	print_by(vec, groups, unique_ids, log_file);
}

void medial::print::print_by_year(const vector<MedSample> &data_records, int year_bin_size, bool unique_ids, bool take_prediction_time, const string &log_file) {
	ofstream fo;
	if (!log_file.empty()) {
		fo.open(log_file, ios::app);
		if (!fo.good())
			MWARN("Warning: can log into file %s\n", log_file.c_str());
	}

	unordered_map<int, int> count_0, count_1;
	vector<int> all_years;
	unordered_set<int> seen_year;
	vector<unordered_map<int, unordered_set<int>>> year_to_seen_pid(2); //of_predicition
	for (size_t i = 0; i < data_records.size(); ++i)
	{
		//int year = int(year_bin_size * round((it->registry.date / 10000) / year_bin_size));
		int label = int(data_records[i].outcome > 0);
		int tm = data_records[i].time;
		if (!take_prediction_time)
			tm = data_records[i].outcomeTime;
		int prediction_year = int(year_bin_size*round(tm / 10000 / year_bin_size));

		if ((label > 0) && seen_year.find(prediction_year) == seen_year.end()) {
			all_years.push_back(prediction_year);
			seen_year.insert(prediction_year);
		}

		if (label > 0) {
			if (!unique_ids ||
				year_to_seen_pid[1][prediction_year].find(data_records[i].id) == year_to_seen_pid[1][prediction_year].end()) {
				++count_1[prediction_year];
				year_to_seen_pid[1][prediction_year].insert(data_records[i].id);
			}
		}
		else {
			if (!unique_ids ||
				year_to_seen_pid[0][prediction_year].find(data_records[i].id) == year_to_seen_pid[0][prediction_year].end()) {
				++count_0[prediction_year];
				year_to_seen_pid[0][prediction_year].insert(data_records[i].id);
			}
		}
	}

	unordered_map<int, int> year_total;
	unordered_map<int, float> year_ratio;
	int i = 0;
	sort(all_years.begin(), all_years.end());
	if (take_prediction_time)
		log_with_file(fo, "Printing by prediction time...\n");
	else
		log_with_file(fo, "Printing by outcome time...\n");
	log_with_file(fo, "Year"  "\t"  "Count_0"  "\t"  "Count_1"  "\t"  "percentage\n");
	for (int year : all_years)
	{
		year_total[year] = count_0[year] + count_1[year];
		year_ratio[year] = count_1[year] / float(count_0[year] + count_1[year]);
		++i;
		log_with_file(fo, "%d\t%d\t%d\t%2.3f%%\n", year, count_0[year], count_1[year],
			100 * count_1[year] / float(count_1[year] + count_0[year]));
	}

	if (fo.good())
		fo.close();
}

void medial::print::print_by_year(const MedSamples &data_records, int year_bin_size, bool unique_ids, bool take_prediction_time, const string &log_file) {
	vector<MedSample> vec;
	for (size_t i = 0; i < data_records.idSamples.size(); ++i)
		vec.insert(vec.end(), data_records.idSamples[i].samples.begin(), data_records.idSamples[i].samples.end());
	print_by_year(vec, year_bin_size, unique_ids, take_prediction_time, log_file);
}

void medial::process::down_sample(MedSamples &samples, double take_ratio, bool with_repeats) {
	int tot_samples = samples.nSamples();
	//int tot_samples = (int)samples.idSamples.size();
	vector<int> pids_index;
	pids_index.reserve(tot_samples);
	for (size_t i = 0; i < samples.idSamples.size(); ++i)
		for (size_t j = 0; j < samples.idSamples[i].samples.size(); ++j)
			pids_index.push_back((int)i);

	int final_cnt = int(take_ratio * tot_samples);
	if (take_ratio >= 1 || final_cnt == 0) {
		return;
	}

	vector<int> all_selected_indexes(final_cnt);
	vector<bool> seen_index(tot_samples);
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<> dist_gen(0, tot_samples - 1);
	MedSamples filterd;
	filterd.time_unit = samples.time_unit;
	vector<vector<MedSample>> new_samples((int)samples.idSamples.size());
	for (size_t k = 0; k < final_cnt; ++k) //for 0 and 1:
	{
		int num_ind = dist_gen(gen);
		if (!with_repeats) {
			while (seen_index[num_ind])
				num_ind = dist_gen(gen);
			seen_index[num_ind] = true;
		}
		int index_i = pids_index[num_ind];
		int index_j = 0;
		while (index_j < samples.idSamples[index_i].samples.size() &&
			num_ind - index_j - 1 >= 0 &&
			pids_index[num_ind - index_j - 1] == index_i)
			++index_j;

		new_samples[index_i].push_back(samples.idSamples[index_i].samples[index_j]);
	}
	for (size_t i = 0; i < new_samples.size(); ++i)
		if (!new_samples[i].empty()) {
			MedIdSamples smp(new_samples[i].front().id);
			smp.samples.swap(new_samples[i]);
			filterd.idSamples.push_back(smp);
		}
	samples.idSamples.swap(filterd.idSamples);
	samples.sort_by_id_date();
	medial::print::print_samples_stats(samples);
}

void medial::process::down_sample_by_pid(MedSamples &samples, double take_ratio, bool with_repeats) {
	if (take_ratio >= 1)
		return;
	unordered_map<int, vector<int>> pid_to_inds;
	vector<int> id_to_pid;
	for (size_t i = 0; i < samples.idSamples.size(); ++i) {
		if (pid_to_inds.find(samples.idSamples[i].id) == pid_to_inds.end())
			id_to_pid.push_back(samples.idSamples[i].id);
		pid_to_inds[samples.idSamples[i].id].push_back((int)i);
	}
	int final_cnt = int(take_ratio * id_to_pid.size());

	vector<int> all_selected_indexes(final_cnt);
	vector<bool> seen_index((int)id_to_pid.size());
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<> dist_gen(0, (int)id_to_pid.size() - 1);
	MedSamples filterd;
	filterd.time_unit = samples.time_unit;
	for (size_t k = 0; k < final_cnt; ++k) //for 0 and 1:
	{
		int num_ind = dist_gen(gen);
		if (!with_repeats) {
			while (seen_index[num_ind])
				num_ind = dist_gen(gen);
			seen_index[num_ind] = true;
		}
		int pid = id_to_pid[num_ind];
		vector<int> take_inds = pid_to_inds.at(pid);
		for (int ind : take_inds)
			//#pragma omp critical
			filterd.idSamples.push_back(samples.idSamples[ind]);
	}
	samples.idSamples.swap(filterd.idSamples);
	samples.sort_by_id_date();
}

double medial::stats::kaplan_meir_on_samples(const vector<MedSample> &incidence_samples, int time_unit, int time_period, const vector<int> *filtered_idx) {
	vector<int> sorted_times;
	vector<bool> all_times(time_period + 1);
	vector<vector<int>> times_indexes;
	sorted_times.reserve(time_period + 1);
	vector<int> final_filter;
	const vector<int> *p_filter = filtered_idx;
	if (filtered_idx == NULL) {
		for (int i = 0; i < incidence_samples.size(); ++i)
			final_filter.push_back(i);
		p_filter = &final_filter;
	}

	double controls = 0, cases = 0, prob = 1;
	double curr_total_ctrls = (double)p_filter->size();

	for (int idx : *p_filter) {
		int time_diff =
			med_time_converter.convert_times(time_unit, global_default_windows_time_unit, incidence_samples[idx].outcomeTime) -
			med_time_converter.convert_times(time_unit, global_default_windows_time_unit, incidence_samples[idx].time);
		if (time_diff > time_period)
			time_diff = time_period;
		if (time_diff < 0)
			continue;
		if (!all_times[time_diff]) {
			sorted_times.push_back(time_diff);
			all_times[time_diff] = true;
		}
	}


	sort(sorted_times.begin(), sorted_times.end());
	times_indexes.resize(sorted_times.size());
	bool warn_show_neg = false, warn_case = false;
	for (int idx : *p_filter) {
		int time_diff =
			med_time_converter.convert_times(time_unit, global_default_windows_time_unit, incidence_samples[idx].outcomeTime) -
			med_time_converter.convert_times(time_unit, global_default_windows_time_unit, incidence_samples[idx].time);
		int original_time = time_diff;
		if (time_diff > time_period)
			time_diff = time_period;
		if (time_diff < 0) {
			if (!warn_show_neg)
				MWARN("Warning - kaplan_meir_on_samples: got negative time. time=%d, outcomeTime=%d\n",
					incidence_samples[idx].time, incidence_samples[idx].outcomeTime);
			warn_show_neg = true;
			continue;
		}
		int ind = medial::process::binary_search_index(sorted_times.data(),
			sorted_times.data() + sorted_times.size() - 1, time_diff);

		if (incidence_samples[idx].outcome <= 0 || original_time <= time_period)
			times_indexes[ind].push_back(idx);
		else {
			if (!warn_case) {
				MWARN("Warning - kaplan_meir_on_samples: got case beyond period time: %d on time %d, outcomeTime %d\n",
					incidence_samples[idx].id, incidence_samples[idx].time, incidence_samples[idx].outcomeTime);
				warn_case = true;
			}
		}
	}

	for (size_t sort_ind = 0; sort_ind < sorted_times.size(); ++sort_ind) {
		const vector<int> &index_order = times_indexes[sort_ind];
		for (int p_i_j : index_order) {
			//keep update kaplan meir in time point
			if (incidence_samples[p_i_j].outcome > 0)
				++cases;
			else
				++controls;
		}
		//reset kaplan meir - flash last time prob

		if (curr_total_ctrls == cases) {
			MWARN("Warning medial::stats::kaplan_meir_on_samples for %d period, in time %d left with %d controls."
				" prob till here %2.4f%% => all controls changed to cases in next period, so stopped here\n",
				time_period, sorted_times[sort_ind], (int)curr_total_ctrls, 100 * (1 - prob));
			break;
		}
		if (curr_total_ctrls < 50)
			MWARN("Warning medial::stats::kaplan_meir_on_samples for %d period, in time %d left with %d controls. prob till here %2.4f%%\n",
				time_period, sorted_times[sort_ind], (int)curr_total_ctrls, 100 * (1 - prob));
		if (curr_total_ctrls > 0 && cases > 0)
			prob *= (curr_total_ctrls - cases) / curr_total_ctrls;
		//MLOG_D("Current Time= %d, total_controls=%d [controls=%d, cases=%d], curr_prob=%2.3f%%\n",
		//	sorted_times[sort_ind], (int)curr_total_ctrls, (int)controls, (int)cases, 100 * (1 - prob));
		curr_total_ctrls -= (controls + cases); //remove controls from current time-window - they are now censored, cases are no longer controls
		controls = 0; cases = 0;
	}
	prob = 1 - prob;

	return prob;
}

double medial::stats::kaplan_meir_on_samples(const MedSamples &incidence_samples, int time_period, const vector<pair<int, int>> *filtered_idx) {
	vector<MedSample> final_samples;

	vector<pair<int, int>> final_filter;
	const vector<pair<int, int>> *p_filter = filtered_idx;
	if (filtered_idx == NULL) {
		for (int i = 0; i < incidence_samples.idSamples.size(); ++i)
			for (int j = 0; j < incidence_samples.idSamples[i].samples.size(); ++j)
				final_filter.push_back(pair<int, int>(i, j));
		p_filter = &final_filter;
	}
	for (const pair<int, int> &idx : *p_filter)
		final_samples.push_back(incidence_samples.idSamples[idx.first].samples[idx.second]);

	return kaplan_meir_on_samples(final_samples, incidence_samples.time_unit, time_period);
}