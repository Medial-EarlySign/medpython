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

	size_t size = 0;
	// id, date, outcome-date
	size += 3 * sizeof(int);

	// Outcome 
	size += sizeof(float);

	// Prediction
	size += sizeof(size_t);
	size += sizeof(float) * prediction.size();

	return size;
}

//.......................................................................................
size_t MedSample::serialize(unsigned char *blob) {

	size_t ptr = 0;
	memcpy(blob + ptr, &id, sizeof(int)); ptr += sizeof(int);
	memcpy(blob + ptr, &date, sizeof(int)); ptr += sizeof(int);
	memcpy(blob + ptr, &outcomeDate, sizeof(int)); ptr += sizeof(int);
	memcpy(blob + ptr, &outcome, sizeof(float)); ptr += sizeof(float);

	size_t nPreds = prediction.size();
	memcpy(blob + ptr, &nPreds, sizeof(size_t)); ptr += sizeof(size_t);

	memcpy(blob + ptr, &(prediction[0]), nPreds*sizeof(float)); ptr += nPreds*sizeof(float);

	return ptr;
}

//.......................................................................................
size_t MedSample::deserialize(unsigned char *blob) {

	size_t ptr = 0;
	memcpy(&id, blob + ptr, sizeof(int)); ptr += sizeof(int);
	memcpy(&date, blob + ptr, sizeof(int)); ptr += sizeof(int);
	memcpy(&outcomeDate, blob + ptr, sizeof(int)); ptr += sizeof(int);
	memcpy(&outcome, blob + ptr, sizeof(float)); ptr += sizeof(float);

	int nPreds; 
	memcpy(&nPreds, blob + ptr, sizeof(size_t)); ptr += sizeof(size_t);

	prediction.resize(nPreds);
	memcpy(&(prediction[0]), blob + ptr, nPreds*sizeof(float)); ptr += nPreds*sizeof(float);

	return ptr;

}


//=======================================================================================
// MedIdSample
//=======================================================================================
// De(Serialize)
//.......................................................................................
size_t MedIdSamples::get_size() {

	size_t size = 0;

	// Id
	size += sizeof(int);

	// Samples
	size += sizeof(size_t);
	for (auto& sample : samples)
		size += sample.get_size();

	return size;
}

//.......................................................................................
size_t MedIdSamples::serialize(unsigned char *blob) {
	
	size_t ptr = 0;
	memcpy(blob + ptr, &id, sizeof(int)); ptr += sizeof(int);

	size_t nSamples = samples.size();
	memcpy(blob + ptr, &nSamples, sizeof(size_t)); ptr += sizeof(size_t);

	for (auto& sample : samples)
		ptr += sample.serialize(blob + ptr);

	return ptr;
}

//.......................................................................................
size_t MedIdSamples::deserialize(unsigned char *blob) {

	size_t ptr = 0;
	memcpy(&id, blob + ptr, sizeof(int)); ptr += sizeof(int);

	size_t nSamples;
	memcpy(&nSamples, blob + ptr, sizeof(size_t)); ptr += sizeof(size_t);

	samples.resize(nSamples);

	for (auto& sample : samples)
		ptr += sample.serialize(blob + ptr);

	return ptr;
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

				if (fields[0] == "EVENT" && fields.size() >= 4) {
					//MLOG("-->### %s %s (%d)\n",fields[0].c_str(),fields[1].c_str(),out.size());
					MedSample sample;
					int id = sample.id = stoi(fields[1]);
					sample.date = stoi(fields[2]);
					sample.outcome = stof(fields[3]);

					if (fields.size() > 5) 
						sample.outcomeDate = stoi(fields[5]);

					if (idSamples.empty() || (id != idSamples.back().id && allIds.find(id) == allIds.end())) {
						MedIdSamples newIdSamples;
						newIdSamples.id = id;
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
				}
			}
		}
	}

	inf.close();
	return 0;
}

// De(Serialize)
//.......................................................................................
size_t MedSamples::get_size() {

	size_t size = 0;
	// Samples
	size += sizeof(size_t);
	for (auto& sample : idSamples)
		size += sample.get_size();

	return size;
}

//.......................................................................................
size_t MedSamples::serialize(unsigned char *blob) {

	size_t ptr = 0;
	size_t nSamples = idSamples.size();
	memcpy(blob + ptr, &nSamples, sizeof(size_t)); ptr += sizeof(size_t);

	for (auto& sample : idSamples) 
		ptr += sample.serialize(blob + ptr);

	return ptr;
}

//.......................................................................................
size_t MedSamples::deserialize(unsigned char *blob) {

	size_t ptr = 0;

	size_t nSamples;
	memcpy(&nSamples, blob + ptr, sizeof(size_t)); ptr += sizeof(size_t);

	idSamples.resize(nSamples);
	for (auto& sample : idSamples)
		ptr += sample.deserialize(blob + ptr);

	return ptr;
}


