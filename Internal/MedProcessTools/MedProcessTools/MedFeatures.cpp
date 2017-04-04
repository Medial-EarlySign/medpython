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
	
	vector<float *> datap;
	for (auto& rec : data)
		datap.push_back((float *)(&rec.second[0]));

#pragma omp parallel for schedule(dynamic)
	for (int i=0; i<(int)datap.size(); i++) {
		for (int j = 0; j < nrows; j++)
			mat(j, i) = datap[i][j];
	}

	// Normalization flag
	mat.normalized_flag = true;
	for (auto& attr : attributes)
		mat.normalized_flag &= (int)attr.second.normalized;
	for (auto& attr : attributes)
		mat.normalized_flag &= (int)attr.second.imputed;

	vector<string> feat_names;
	get_feature_names(feat_names);
	mat.signals.insert(mat.signals.end(), feat_names.begin(), feat_names.end());
	for (auto& ss : samples) {
		RecordData rd;
		rd.time = 0L;  rd.weight = rd.pred = rd.label = 0.0; rd.split = 0;
		rd.id = ss.id;
		rd.date = ss.time;
		mat.recordsMetadata.push_back(rd);
	}
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
	return MedSerialize::get_size(data, samples);
}

//.......................................................................................
size_t  MedFeatures::serialize(unsigned char *blob) {
	return MedSerialize::serialize(blob, data, samples);
}

//.......................................................................................
size_t MedFeatures::deserialize(unsigned char *blob) {
	return MedSerialize::deserialize(blob, data, samples);
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

//.......................................................................................
int MedFeatures::write_as_csv_mat(const string &csv_fname)
{
	ofstream out_f;

	out_f.open(csv_fname);

	if (!out_f.is_open()) {
		MERR("ERROR: MedFeatures::write_as_csv_mat() :: Can't open file %s for writing\n", csv_fname.c_str());
		return -1;
	}

	vector<string> col_names;
	get_feature_names(col_names);
	int n_preds = 0;

	// names line
	// serial
	out_f << "serial";
	// samples
	out_f << ",id,time,outcome,outcome_time";
//	MLOG("samples.size %d , preds %d\n", samples.size(), samples[0].prediction.size());
	if (samples.size() > 0 && samples[0].prediction.size() > 0) {
		n_preds = (int)samples[0].prediction.size();
		//MLOG("n_preds = %d\n", n_preds);
		for (int j=0; j<n_preds; j++)
			out_f << ",pred_"+to_string(j);
	}
	// names of features
	for (int j=0; j<col_names.size(); j++)
		out_f << "," + col_names[j];
	out_f << "\n";

	for (int i=0; i<samples.size(); i++) {

		// serial
		out_f << to_string(i);

		// sample
		out_f << "," + to_string(samples[i].id);
		out_f << "," + to_string(samples[i].time);
		out_f << "," + to_string(samples[i].outcome);
		out_f << "," + to_string(samples[i].outcomeTime);
		for (int j=0; j<n_preds; j++)
			out_f << "," + to_string(samples[i].prediction[j]);

		// features
		for (int j=0; j<col_names.size(); j++)
			out_f << "," + to_string(data[col_names[j]][i]);
		out_f << "\n";
	}

	out_f.close();
	return 0;
}