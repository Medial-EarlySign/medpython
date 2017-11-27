#define _CRT_SECURE_NO_WARNINGS

#include "MedFeatures.h"

#include "Logger/Logger/Logger.h"
#define LOCAL_SECTION LOG_MEDFEAT
#define LOCAL_LEVEL	LOG_DEF_LEVEL

int MedFeatures::global_serial_id_cnt = 0;

//=======================================================================================
// MedFeatures
//=======================================================================================
// Get a vector of feature names
//.......................................................................................
void MedFeatures::get_feature_names(vector<string>& names) {

	names.resize(data.size());

	int i = 0;
	for (auto& rec : data)
		names[i++] = rec.first;
}

// Get data (+attributes) as matrix
//.......................................................................................
void MedFeatures::get_as_matrix(MedMat<float>& mat) {
	vector<string> dummy_names;
	get_as_matrix(mat, dummy_names);
}

// Get subset of data (+attributes) as matrix : Only features in 'names'
//.......................................................................................
void MedFeatures::get_as_matrix(MedMat<float>& mat, vector<string>& names) {

	// Which Features to take ?
	vector<string> namesToTake;
	if (names.size())
		namesToTake = names;
	else
		get_feature_names(namesToTake);


	int ncols = (int)namesToTake.size();
	int nrows = (int)samples.size();

	mat.resize(nrows, ncols);

	vector<float *> datap;
	for (string& name : namesToTake)
		datap.push_back((float *)(&data[name][0]));

#pragma omp parallel for schedule(dynamic)
	for (int i=0; i<(int)datap.size(); i++) {
		for (int j = 0; j < nrows; j++) {
			if (!isfinite(datap[i][j])) {
				MTHROW_AND_ERR("nan in col [%s] in record [%d]", namesToTake[i].c_str(), j);
			}
			mat(j, i) = datap[i][j];
		}
	}

	// Normalization flag
	mat.normalized_flag = true;
	for (string& name : namesToTake)
		mat.normalized_flag &= (int)attributes[name].normalized;
	for (string& name : namesToTake)
		mat.normalized_flag &= (int)attributes[name].imputed;

	mat.signals.insert(mat.signals.end(), namesToTake.begin(), namesToTake.end());
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

// Get subset of data (+attributes) as matrix: Only features in 'names' and rows in 'idx'
//.......................................................................................
void MedFeatures::get_as_matrix(MedMat<float> &mat, const vector<string> &names, vector<int> &idx)
{
	// Which Features to take ?
	vector<string> namesToTake;
	if (names.size())
		namesToTake = names;
	else
		get_feature_names(namesToTake);


	int ncols = (int)namesToTake.size();
	int nrows = (int)idx.size();

	mat.resize(nrows, ncols);

	vector<float *> datap;
	for (string& name : namesToTake)
		datap.push_back((float *)(&data[name][0]));

#pragma omp parallel for schedule(dynamic)
	for (int i=0; i<(int)datap.size(); i++) {
		for (int j = 0; j < nrows; j++) {
			int jj = idx[j];
			if (!isfinite(datap[i][jj])) {
				MTHROW_AND_ERR("nan in col [%s] in record [%d]", namesToTake[i].c_str(), jj);
			}
			mat(j, i) = datap[i][jj];
		}
	}

	// Normalization flag
	mat.normalized_flag = true;
	for (string& name : namesToTake)
		mat.normalized_flag &= (int)attributes[name].normalized;
	for (string& name : namesToTake)
		mat.normalized_flag &= (int)attributes[name].imputed;

	mat.signals.insert(mat.signals.end(), namesToTake.begin(), namesToTake.end());
	for (int i=0; i<nrows; i++) {
	//for (auto& ss : samples) {
		auto &ss = samples[idx[i]];
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

// Append samples at end of samples vector
//.......................................................................................
void MedFeatures::append_samples(MedIdSamples& in_samples) {
	samples.insert(samples.end(), in_samples.samples.begin(), in_samples.samples.end());
}

// Insert samples at position idex, assuming samples vector is properly allocated
//.......................................................................................
void MedFeatures::insert_samples(MedIdSamples& in_samples, int index) {

	for (unsigned int i = 0; i < in_samples.samples.size(); i++)
		samples[index + i] = in_samples.samples[i];
}

// initialize pid_pos_len vector according to samples
//.......................................................................................
void MedFeatures::init_pid_pos_len()
{
	int curr_pid = -1, curr_len = 0 , curr_pos = -1;

	int i = -1;
	for (auto& sample : samples) {
		i++;

		// New pid
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

// Calculate a crc for the data (used for debugging mainly)
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

// MLOG data in csv format 
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

// Write features (samples + weights + data) as csv with a header line
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

	// header line
	out_f << "serial"; // serial
	if (weights.size()) out_f << ",weight"; // Weight (if given)
	out_f << ",id,time,outcome,outcome_time,split"; // samples

	// Predictions
	if (samples.size() > 0 && samples[0].prediction.size() > 0) {
		n_preds = (int)samples[0].prediction.size();
		for (int j=0; j<n_preds; j++)
			out_f << ",pred_"+to_string(j);
	}
	
	// names of features
	for (int j=0; j<col_names.size(); j++)
		out_f << "," + col_names[j];
	out_f << "\n";

	// data
	for (int i=0; i<samples.size(); i++) {

		out_f << to_string(i); // serial
		if (weights.size()) out_f << "," + to_string(weights[i]); // Weights

		// sample
		out_f << "," + to_string(samples[i].id);
		out_f << "," + to_string(samples[i].time);
		out_f << "," + to_string(samples[i].outcome);
		out_f << "," + to_string(samples[i].outcomeTime);
		out_f << "," + to_string(samples[i].split);

		// predictions
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

// Read features (samples + weights + data) from a csv file with a header line
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
	int weighted = 0;
	while (getline(inf, curr_line)) {
		boost::trim(curr_line);
		vector<string> fields;
		boost::split(fields, curr_line, boost::is_any_of(","));
		int idx = 0;
		if (ncols == -1) { // Header line	
			assert(fields[idx++].compare("serial") == 0); 
			if (fields[idx].compare("weight") == 0) {
				weighted = 1; idx++;
			}
			assert(fields[idx++].compare("id") == 0);
			assert(fields[idx++].compare("time") == 0);
			assert(fields[idx++].compare("outcome") == 0);
			assert(fields[idx++].compare("outcome_time") == 0);
			assert(fields[idx++].compare("split") == 0);

			for (int i = idx; i < fields.size(); i++) {
				data[fields[i]] = vector<float>();
				attributes[fields[i]].normalized = attributes[fields[i]].imputed = false;
				names.push_back(fields[i]);
			}
			ncols = (int)fields.size();
		}
		else { // Data lines
			if (fields.size() != ncols) {
				string msg = "expected " + to_string(ncols) + " fields, got " + to_string((int)fields.size()) + "fields in line: " + curr_line.c_str() + "\n";
				throw runtime_error(msg.c_str());
			}

			if (weighted)
				weights.push_back(stof(fields[idx++]));

			MedSample newSample;
			newSample.id = stoi(fields[idx++]);
			newSample.time = stoi(fields[idx++]);
			newSample.outcome = stof(fields[idx++]);
			newSample.outcomeTime = stoi(fields[idx++]);
			newSample.split = stoi(fields[idx++]);
			samples.push_back(newSample);

			for (int i = 0; i < names.size(); i++)
				data[names[i]].push_back(stof(fields[idx + i]));
		}
	}

	inf.close();
	return 0;
}

// Filter data (and attributes) to include only selected features
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

// (De)Serialization
//.......................................................................................
size_t MedFeatures::get_size() {
	return MedSerialize::get_size(data, weights, samples, attributes);
}

//.......................................................................................
size_t  MedFeatures::serialize(unsigned char *blob) {
	return MedSerialize::serialize(blob, data, weights, samples, attributes);
}

//.......................................................................................
size_t MedFeatures::deserialize(unsigned char *blob) {
	return MedSerialize::deserialize(blob, data, weights, samples, attributes);
}