#define _CRT_SECURE_NO_WARNINGS

#include "MedFeatures.h"

#include "Logger/Logger/Logger.h"
#define LOCAL_SECTION LOG_MEDFEAT
#define LOCAL_LEVEL	LOG_DEF_LEVEL

//=======================================================================================
// MedFeatures
//=======================================================================================
//.......................................................................................
void MedFeatures::get_feature_names(vector<string>& names) {

	names.resize(data.size());

	int i = 0;
	for (auto& rec : data)
		names[i++] = rec.first;
}

//.......................................................................................
void MedFeatures::get_as_matrix(MedMat<float>& mat) {

	int ncols = (int) data.size();
	int nrows = (int)samples.size();

	mat.resize(nrows, ncols);
	
	int i = 0;
	for (auto& rec : data) {
		for (int j = 0; j < nrows; j++)
			mat(j, i) = rec.second[j];
		i++;
	}

	// Normalization flag
	mat.normalized_flag = true;
	for (auto& attr : attributes)
		mat.normalized_flag &= (int)attr.second.normalized;
}

//.......................................................................................
void MedFeatures::append_samples(MedIdSamples& in_samples) {
	samples.insert(samples.end(), in_samples.samples.begin(), in_samples.samples.end());
}

//.......................................................................................
void MedFeatures::insert_samples(MedIdSamples& in_samples, int index) {

	for (unsigned int i = 0; i < in_samples.samples.size(); i++)
		samples[index + i] = in_samples.samples[i];
}

// (De)Serialization
//.......................................................................................
size_t MedFeatures::get_size() {

	size_t size = 0;

	// Number of Features and samples
	size += 2*sizeof(int);

	// Data
	for (auto& rec : data) {
		// name
		size += sizeof(size_t);
		size += rec.first.length() + 1;
		// data
		size += rec.second.size() * sizeof(float);
	}

	// Samples
	size += samples.size() * sizeof(MedSample);

	return size;
}

extern char signalName_c[MAX_NAME_LEN + 1];

//.......................................................................................
size_t  MedFeatures::serialize(unsigned char *blob) {

	size_t ptr = 0;
	
	// Number of Features
	int nFeatures = (int) data.size();
	memcpy(blob + ptr, &nFeatures, sizeof(int)); ptr += sizeof(int);

	// Features size
	int nSamples = (int) samples.size();
	memcpy(blob + ptr, &nSamples, sizeof(int)); ptr += sizeof(int);


	// Data
	for (auto& rec : data) {
		// name
		size_t nameLen = rec.first.length();
		assert(nameLen < MAX_NAME_LEN);

		strcpy(signalName_c, rec.first.c_str());

		memcpy(blob + ptr, &nameLen, sizeof(size_t)); ptr += sizeof(size_t);
		memcpy(blob + ptr, signalName_c, nameLen + 1); ptr += nameLen + 1;

		// data
		memcpy(blob + ptr, &(rec.second[0]), nSamples * sizeof(float)); ptr += nSamples * sizeof(float);
	}

	// Samples
	for (unsigned int i = 0; i < samples.size(); i++)
		ptr += samples[i].serialize(blob+ptr);

	return ptr;
}

//.......................................................................................
size_t MedFeatures::deserialize(unsigned char *blob) {

	size_t ptr = 0;

	// Number of Features
	int nFeatures;
	memcpy(&nFeatures, blob + ptr, sizeof(int)); ptr += sizeof(int);

	// Number of Samples
	int nSamples;
	memcpy(&nSamples, blob + ptr, sizeof(int)); ptr += sizeof(int);

	// Data
	for (int i = 0; i < nFeatures; i++) {
		// name
		size_t nameLen;
		memcpy(&nameLen, blob + ptr, sizeof(size_t)); ptr += sizeof(size_t);
		assert(nameLen < MAX_NAME_LEN);

		memcpy(signalName_c, blob + ptr, nameLen + 1); ptr += nameLen + 1;

		// data
		data[signalName_c] = vector<float>(nSamples);
		memcpy(blob + ptr, &(data[signalName_c][0]), nSamples * sizeof(float)); ptr += nSamples * sizeof(float);
	}

	// Samples
	samples.resize(nSamples);
	for (unsigned int i = 0; i < samples.size(); i++)
		ptr += samples[i].deserialize(blob + ptr);

	return ptr;
}

//.......................................................................................
void MedFeatures::init_pid_pos_len()
{
	int curr_pid = -1, curr_len = 0 , curr_pos = -1;

	int i = -1;
	for (auto& sample : samples) {
		i++;
		if (curr_pid != sample.id) {
			if (curr_len > 0)
				pid_pos_len[curr_pid] = make_pair(curr_pos, curr_len);
			curr_pid = sample.id;
			curr_pos = i;
			curr_len = 0;
		}
		curr_len++;
	}
	// last one
	pid_pos_len[curr_pid] = make_pair(curr_pos, curr_len);
}

//.......................................................................................
unsigned int MedFeatures::get_crc()
{
	int i, j;
	int byte, crc;
	int mask;

	i = 0;
	crc = 0xFFFFFFFF;
	for (auto& sig_data : data) {
		unsigned char *msg = (unsigned char *)&sig_data.second[0];
		int len = (int)(sig_data.second.size() * sizeof(float));
		for (i=0; i<len; i++) {
			byte = msg[i];            // Get next byte.
			crc = crc ^ byte;
			for (j = 7; j >= 0; j--) {    // Do eight times.
				mask = -(crc & 1);
				crc = (crc >> 1) ^ (0xEDB88320 & mask);
			}
		}
	}

	return (unsigned int)(~crc);
}

//.......................................................................................
void MedFeatures::print_csv()
{
	for (auto &vec : data) {
		MLOG("%s :: ", vec.first.c_str());
		for (auto v : vec.second)
			MLOG("%f,", v);
		MLOG("\n");
	}
}