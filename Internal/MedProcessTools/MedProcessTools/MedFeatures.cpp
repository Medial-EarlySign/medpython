#define _CRT_SECURE_NO_WARNINGS

#include "MedFeatures.h"

#include "Logger/Logger/Logger.h"
#define LOCAL_SECTION LOG_MEDFEAT
#define LOCAL_LEVEL	LOG_DEF_LEVEL

int MedFeatures::global_serial_id_cnt = 0;

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
		for (int j = 0; j < nrows; j++) {
			if (!isfinite(datap[i][j])) {
				vector<string> fnames;
				get_feature_names(fnames);
				MTHROW_AND_ERR("nan in col [%s] in record [%d]", fnames[i].c_str(), j);
			}
			mat(j, i) = datap[i][j];
		}
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
		rd.time = 0L;  rd.weight = 0.0;  rd.label = ss.outcome; ; rd.split = ss.split;
		if (ss.prediction.size() == 1)
			rd.pred = ss.prediction[0];
		else rd.pred = 0.0;
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

//.......................................................................................
int MedFeatures::read_from_csv_mat(const string &csv_fname)
{
	if (!file_exists(csv_fname)) {
		fprintf(stderr, "File %s doesn't exist\n", csv_fname.c_str());
		throw exception();
	}
	fprintf(stderr, "reading data from %s\n", csv_fname.c_str());
	ifstream inf;
	inf.open(csv_fname, ios::in);
	if (!inf) {
		cerr << "can not open file\n";
		throw exception();
	}
	
	int ncols = -1;
	string curr_line;
	vector<string> names; 
	while (getline(inf, curr_line)) {
		boost::trim(curr_line);
		vector<string> fields;
		boost::split(fields, curr_line, boost::is_any_of(","));
		if (ncols == -1) {
			assert(fields[0].compare("serial") == 0);
			assert(fields[1].compare("id") == 0);
			assert(fields[2].compare("time") == 0);
			assert(fields[3].compare("outcome") == 0);
			assert(fields[4].compare("outcome_time") == 0);

			for (int i = 5; i < fields.size(); i++) {
				data[fields[i]] = vector<float>();
				attributes[fields[i]].normalized = attributes[fields[i]].imputed = false;
				names.push_back(fields[i]);
			}
			ncols = fields.size();
		}
		else {
			if (fields.size() != ncols) {
				string msg = "expected " + to_string(ncols) + " fields, got " + to_string((int)fields.size()) + "fields in line: " + curr_line.c_str() + "\n";
				throw runtime_error(msg.c_str());
			}

			MedSample newSample;
			newSample.id = stoi(fields[1]);
			newSample.time = stoi(fields[2]);
			newSample.outcome = stof(fields[3]);
			newSample.outcomeTime = stoi(fields[4]);
			samples.push_back(newSample);

			for (int i = 0; i < names.size(); i++)
				data[names[i]].push_back(stof(fields[5 + i]));
		}
	}

	inf.close();
	return 0;
}


// Filter set of features
//.......................................................................................
int MedFeatures::filter(unordered_set<string>& selectedFeatures) {

	// Sanity
	for (string feature : selectedFeatures) {
		if (data.find(feature) == data.end()) {
			MERR("Cannot find feature %s in Matrix\n", feature.c_str());
			return -1;
		}
	}

	// Cleaning
	vector<string> removedFeatures;
	for (auto& rec : data) {
		string feature = rec.first;
		if (selectedFeatures.find(feature) == selectedFeatures.end())
			removedFeatures.push_back(feature);
	}

	for (string& feature : removedFeatures) {
		data.erase(feature);
		attributes.erase(feature);
	}

	return 0;
}