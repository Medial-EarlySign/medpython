#include "MedSamples.h"
#include "MedProcessTools/MedProcessTools/MedFeatures.h"
#include "Logger/Logger/Logger.h"

#define LOCAL_SECTION MED_SAMPLES_CV
#define LOCAL_LEVEL	LOG_DEF_LEVEL

//=======================================================================================
// MedSample
//=======================================================================================
// (De)Serialization

//.......................................................................................
size_t MedSample::get_size() {
	return MedSerialize::get_size(id, time, outcomeTime, outcome, prediction);
}

//.......................................................................................
size_t MedSample::serialize(unsigned char *blob) {
	return MedSerialize::serialize(blob, id, time, outcomeTime, outcome, prediction);
}

//.......................................................................................
size_t MedSample::deserialize(unsigned char *blob) {
	return MedSerialize::deserialize(blob, id, time, outcomeTime, outcome, prediction);
}

//.......................................................................................
void MedSample::print(const string prefix) {
	MLOG("%s :: id %d time %d outcomeTime %d outcome %f prediction(%d)", prefix.c_str(), id, time, outcomeTime, outcome, prediction.size());
	if (prediction.size() > 0)
		for (auto pred : prediction)
			MLOG(" %f", pred);
	MLOG("\n");
}

//=======================================================================================
// MedIdSample
//=======================================================================================
// De(Serialize)
//.......................................................................................
size_t MedIdSamples::get_size() {
	return MedSerialize::get_size(id, split, samples);
}

//.......................................................................................
size_t MedIdSamples::serialize(unsigned char *blob) {
	return MedSerialize::serialize(blob, id, split, samples);
}

//.......................................................................................
size_t MedIdSamples::deserialize(unsigned char *blob) {
	return MedSerialize::deserialize(blob, id, split, samples);
}

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

void MedSamples::get_preds(vector<float>& preds) {
	for (auto& idSample : idSamples) 
		for (auto& sample : idSample.samples) 
			for (int i = 0; i < sample.prediction.size(); i++)
				preds.push_back(sample.prediction[i]);
}

void MedSamples::get_y(vector<float>& y) {
	for (auto& idSample : idSamples)
		for (auto& sample : idSample.samples)
			for (int i = 0; i < sample.prediction.size(); i++)
				y.push_back(sample.outcome);
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
	unordered_set<int> allIds;

	int samples = 0;
	while (getline(inf, curr_line)) {
		//MLOG("--> %s\n",curr_line.c_str());
		if ((curr_line.size() > 1) && (curr_line[0] != '#')) {

			if (curr_line[curr_line.size() - 1] == '\r')
				curr_line.erase(curr_line.size() - 1);

			vector<string> fields;
			split(fields, curr_line, boost::is_any_of("\t"));

			if (fields.size() >= 2) {

				if (fields[0] == "NAME") MLOG("reading NAME = %s\n", fields[1].c_str());
				if (fields[0] == "DESC") MLOG("reading DESC = %s\n", fields[1].c_str());
				if (fields[0] == "TYPE") MLOG("reading TYPE = %s\n", fields[1].c_str());
				if (fields[0] == "NCATEG")  MLOG("reading NCATEG = %s\n", fields[1].c_str());
				if (fields[0] == "EVENT_FIELDS") {
					assert(fields[1] == "id");
					assert(fields[2] == "time");
					assert(fields[3] == "outcome");
					assert(fields[4] == "outcomeLength");
					assert(fields[5] == "outcomeTime");
					assert(fields[6] == "split");
					assert(fields[7] == "prediction");
				}
				if (fields[0] == "EVENT" && fields.size() >= 4) {
					//MLOG("-->### %s %s (%d)\n",fields[0].c_str(),fields[1].c_str(),out.size());
					MedSample sample;
					int id = sample.id = stoi(fields[1]);
					sample.time = stoi(fields[2]);
					sample.outcome = stof(fields[3]);
					int split = -1;
					if (fields.size() > 5) 
						sample.outcomeTime = stoi(fields[5]);
					if (fields.size() > 6)
						split = stoi(fields[6]);
					if (fields.size() > 7)
						sample.prediction.push_back(stof(fields[7]));

					if (idSamples.empty() || (id != idSamples.back().id && allIds.find(id) == allIds.end())) {
						MedIdSamples newIdSamples(id);
						newIdSamples.samples.push_back(sample);
						idSamples.push_back(newIdSamples);
						allIds.insert(id);
					}
					else if (id == idSamples.back().id) {
						idSamples.back().samples.push_back(sample);
					}
					else {
						MERR("Id %d appears not consecutively\n", id);
						return -1;
					}
					samples++;
				}
			}
		}
	}
	MLOG("read [%d] samples for [%d] patient IDs\n", samples, idSamples.size());
	inf.close();
	return 0;
}

int MedSamples::write_to_file(const string &fname)
{
	ofstream of(fname);

	MLOG("MedSamples: writing to %s\n", fname.c_str());
	if (!of) {
		MERR("MedSamples: can't open file %s for writing\n", fname.c_str());
		return -1;
	}
	int samples = 0;

	of << "EVENT_FIELDS" << '\t' << "id" << '\t' << "time" << '\t' << "outcome" << '\t' << "outcomeLength" <<
		'\t' << "outcomeTime" << '\t' << "split" << '\t' << "prediction" << endl;

	for (auto &s: idSamples) {
		for (auto ss : s.samples) {
			samples++;
			of << "EVENT" << '\t' << ss.id << '\t' << ss.time << '\t' << ss.outcome << '\t' << 100000 << '\t' <<
				ss.outcomeTime << '\t' << s.split << '\t' << ss.prediction.front() << endl;
		}
	}

	MLOG("wrote [%d] samples for [%d] patient IDs\n", samples, idSamples.size());
	of.close();
	return 0;
}

// De(Serialize)
//.......................................................................................
size_t MedSamples::get_size() {
	return MedSerialize::get_size(time_unit, idSamples);
}

//.......................................................................................
size_t MedSamples::serialize(unsigned char *blob) {
	return MedSerialize::serialize(blob, time_unit, idSamples);
}

//.......................................................................................
size_t MedSamples::deserialize(unsigned char *blob) {
	return MedSerialize::deserialize(blob, time_unit, idSamples);
}


