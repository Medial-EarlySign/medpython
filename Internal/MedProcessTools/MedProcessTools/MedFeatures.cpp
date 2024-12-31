#define _CRT_SECURE_NO_WARNINGS

#include <MedProcessTools/MedProcessTools/MedFeatures.h>
#include <MedUtils/MedUtils/MedUtils.h>
#include <MedUtils/MedUtils/MedGenUtils.h>
#include <random>
#include <omp.h>
#include <MedIO/MedIO/MedIO.h>
#include <MedUtils/MedUtils/MedGlobalRNG.h>

#include <Logger/Logger/Logger.h>
#define LOCAL_SECTION LOG_MEDFEAT
#define LOCAL_LEVEL	LOG_DEF_LEVEL

int MedFeatures::global_serial_id_cnt = 0;

void print_stats(const MedFeatures &features, const vector<int> &sel_idx, const vector<string> &group_values) {
	unordered_map<string, vector<int>> counts_stat;
	if (!sel_idx.empty()) {
		for (size_t i = 0; i < sel_idx.size(); ++i)
		{
			if (counts_stat[group_values[sel_idx[i]]].empty())
				counts_stat[group_values[sel_idx[i]]].resize(2);
			++counts_stat[group_values[sel_idx[i]]][features.samples[i].outcome > 0];
		}
	}
	else {
		for (size_t i = 0; i < group_values.size(); ++i)
		{
			if (counts_stat[group_values[i]].empty())
				counts_stat[group_values[i]].resize(2);
			++counts_stat[group_values[i]][features.samples[i].outcome > 0];
		}
	}
	MLOG("Group\tCount_0\tCount_1\tratio\n");
	vector<string> all_groups;
	all_groups.reserve(counts_stat.size());
	for (const auto &it : counts_stat)
		all_groups.push_back(it.first);
	sort(all_groups.begin(), all_groups.end());

	for (const string &grp : all_groups)
		MLOG("%s\t%d\t%d\t%2.5f\n", grp.c_str(),
			counts_stat[grp][0], counts_stat[grp][1], counts_stat[grp][1] / double(counts_stat[grp][1] + counts_stat[grp][0]));
}

//=======================================================================================
// MedFeatures
//=======================================================================================
// Get a vector of feature names
//.......................................................................................
void MedFeatures::get_feature_names(vector<string>& names) const {

	names.resize(data.size());

	int i = 0;
	for (auto& rec : data)
		names[i++] = rec.first;
}

// Get data (+attributes) as matrix
//.......................................................................................
void MedFeatures::get_as_matrix(MedMat<float>& mat) const {
	vector<string> dummy_names;
	get_as_matrix(mat, dummy_names);
}

// Get subset of data (+attributes) as matrix : Only features in 'names'
//.......................................................................................
void MedFeatures::get_as_matrix(MedMat<float>& mat, vector<string>& names) const {

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
		datap.push_back((float *)(&data.at(name)[0]));


	MedTimer time_me;
	time_me.start();


	// next algorithm is built to reduce cache misses on transpose
	// the idea is to split the work into batches of more or less rectangular matrices
	// instead of doing whole rows each time
	// we assume at the moment that the number of features is a good number to work with...
	// we also want each batch to handle ~1MB which is ~250k matrix size
	vector<int> batches;
	if (datap.size() > 0) {
		int nelements_in_batch = 250000;
		int n_batch = 1 + nelements_in_batch / (int)datap.size();
		if (n_batch > nrows) n_batch = nrows;
		int curr = 0;
		while (curr < nrows) {
			batches.push_back(curr);
			curr += n_batch;
		}
		if ((batches.size() == 0) || (batches.back() < nrows))
			batches.push_back(nrows);

#pragma omp parallel for schedule(dynamic)
		for (int b = 0; b < batches.size() - 1; b++) {

			for (int i = 0; i < (int)datap.size(); i++) {
				for (int j = batches[b]; j < batches[b + 1]; j++) {
					if (!isfinite(datap[i][j])) {
						MTHROW_AND_ERR("nan in col [%s] in record [%d]", namesToTake[i].c_str(), j);
					}
					mat(j, i) = datap[i][j];
				}
			}

		}
	}

	/*
	#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < (int)datap.size(); i++) {
			for (int j = 0; j < nrows; j++) {
				if (!isfinite(datap[i][j])) {
					MTHROW_AND_ERR("nan in col [%s] in record [%d]", namesToTake[i].c_str(), j);
				}
				mat(j, i) = datap[i][j];
			}
		}
	*/
	time_me.take_curr_time();
	MLOG_D("Matrix transpose time is %f sec\n", time_me.diff_sec());
	//Test:
	for (const string& name : namesToTake)
		if (attributes.find(name) == attributes.end())
			MTHROW_AND_ERR("Error feature \"%s\" is missing attribute information\n", name.c_str());
	// Normalization flag
	mat.normalized_flag = true;
	for (const string& name : namesToTake)
		mat.normalized_flag &= (int)attributes.at(name).normalized;
	for (const string& name : namesToTake)
		mat.normalized_flag &= (int)attributes.at(name).imputed;

	mat.signals.insert(mat.signals.end(), namesToTake.begin(), namesToTake.end());
	int index = 0;
	//mat.time_unit = time_unit;
	for (auto& ss : samples) {
		RecordData rd;
		rd.outcomeTime = (long)ss.outcomeTime; rd.label = ss.outcome; ; rd.split = ss.split; rd.weight = 0.0;
		if (index < weights.size())
			rd.weight = weights[index];
		if (ss.prediction.size() == 1)
			rd.pred = ss.prediction[0];
		else rd.pred = 0.0;
		rd.id = ss.id;
		rd.date = ss.time;
		mat.recordsMetadata.push_back(rd);
		++index;
	}
}

// Get subset of data (+attributes) as matrix: Only features in 'names' and rows in 'idx'
//.......................................................................................
void MedFeatures::get_as_matrix(MedMat<float> &mat, const vector<string> &names, vector<int> &idx) const
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
		datap.push_back((float *)(&data.at(name)[0]));

#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < (int)datap.size(); i++) {
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
		mat.normalized_flag &= (int)attributes.at(name).normalized;
	for (string& name : namesToTake)
		mat.normalized_flag &= (int)attributes.at(name).imputed;

	int index = 0;
	mat.signals.insert(mat.signals.end(), namesToTake.begin(), namesToTake.end());
	//mat.time_unit = time_unit;
	for (int i = 0; i < nrows; i++) {
		//for (auto& ss : samples) {
		auto &ss = samples[idx[i]];
		RecordData rd;
		rd.outcomeTime = (long)ss.outcomeTime;  rd.weight = 0.0;  rd.label = ss.outcome; ; rd.split = ss.split;
		if (index < weights.size())
			rd.weight = weights[index];
		if (ss.prediction.size() == 1)
			rd.pred = ss.prediction[0];
		else rd.pred = 0.0;
		rd.id = ss.id;
		rd.date = ss.time;
		mat.recordsMetadata.push_back(rd);
	}
}

void MedFeatures::set_as_matrix(const MedMat<float>& mat) {
	const vector<string> &namesToTake = mat.signals;

	for (int i = 0; i < (int)mat.ncols; ++i) {
		data[namesToTake[i]].resize(mat.nrows);
		mat.get_col(i, data[namesToTake[i]]);
	}

	// Normalization flag
	for (const string& name : namesToTake)
		attributes[name].normalized = mat.normalized_flag > 0;
	for (const string& name : namesToTake)
		attributes[name].imputed = mat.normalized_flag > 0;

	weights.reserve((int)mat.recordsMetadata.size());
	bool no_zero_weight = false;
	for (auto& rd : mat.recordsMetadata) {
		MedSample smp;
		smp.id = rd.id; smp.outcome = rd.label; smp.time = rd.date;
		smp.split = rd.split;
		smp.prediction.push_back(rd.pred);
		smp.outcomeTime = rd.outcomeTime;
		weights.push_back(rd.weight);
		if (!no_zero_weight)
			no_zero_weight = rd.weight > 0;

		samples.push_back(smp);
	}
	if (!no_zero_weight)
		weights.clear();
	init_pid_pos_len();
}

// Append samples at end of samples vector (used for generating samples set before generating features)
//.......................................................................................
void MedFeatures::append_samples(MedIdSamples& in_samples) {
	samples.insert(samples.end(), in_samples.samples.begin(), in_samples.samples.end());
}

// Insert samples at position idex, assuming samples vector is properly allocated  (used for generating samples set before generating features)
//.......................................................................................
void MedFeatures::insert_samples(MedIdSamples& in_samples, int index) {

	for (unsigned int i = 0; i < in_samples.samples.size(); i++)
		samples[index + i] = in_samples.samples[i];
}

// initialize pid_pos_len vector according to samples
//.......................................................................................
void MedFeatures::init_pid_pos_len()
{
	int curr_pid = -1, curr_len = 0, curr_pos = -1;

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
		for (i = 0; i < len; i++) {
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
void MedFeatures::print_csv() const
{
	for (auto &vec : data) {
		MLOG("%s :: ", vec.first.c_str());
		for (auto v : vec.second)
			MLOG("%f,", v);
		MLOG("\n");
	}
}

// Write features (samples + weights + data) as csv with a header line
// Return -1 upon failure to open file, 0 upon success
//.......................................................................................
void MedFeatures::write_csv_data(ofstream& out_f, bool write_attributes, vector<string>& col_names, int start_idx) const {

	for (int i = 0; i < samples.size(); i++) {

		out_f << to_string(i + start_idx); // serial
		if (weights.size()) out_f << "," << weights[i]; // Weights

														// sample
		string sample_str;

		samples[i].write_to_string(sample_str, time_unit, write_attributes, ",");
		boost::replace_all(sample_str, "SAMPLE,", "");

		out_f << "," << sample_str;

		// features
		for (int j = 0; j < col_names.size(); j++)
			out_f << "," << data.at(col_names[j])[i];
		out_f << "\n";

	}
}

int MedFeatures::write_as_csv_mat(const string &csv_fname, bool write_attributes) const
{
	ofstream out_f;

	// Sanity - if write_attributes is true, all samples must have the same attributes
	set<string> attr_names, str_attr_names;
	if (write_attributes && !samples.empty()) {
		vector<string> attributes_nm, attributes_str_names;
		samples[0].get_all_attributes(attributes_nm, attributes_str_names);

		int nAttr = (int)attributes_nm.size();
		int nStrAttr = (int)attributes_str_names.size();

		for (const string& attr : attributes_nm)
			attr_names.insert(attr);

		for (const string& attr : attributes_str_names)
			str_attr_names.insert(attr);

		for (unsigned int i = 1; i < samples.size(); i++) {
			vector<string> attributes_nm_2, attributes_str_names_2;
			samples[i].get_all_attributes(attributes_nm_2, attributes_str_names_2);
			if (attributes_nm_2.size() != nAttr || attributes_str_names_2.size() != nStrAttr) {
				MERR("Attrributes # inconsistency betweens samples %d and 0\n", i);
				return -1;
			}

			for (const string& attr : attributes_nm_2) {
				if (attr_names.find(attr) == attr_names.end()) {
					MERR("Attrributes names inconsistency betweens samples %d and 0 : extra attribute %s\n", i, attr.c_str());
					return -1;
				}
			}

			for (const string &attr : attributes_str_names_2) {
				if (str_attr_names.find(attr) == str_attr_names.end()) {
					MERR("Attrributes names inconsistency betweens samples %d and 0 : extra attribute %s\n", i, attr.c_str());
					return -1;
				}
			}

		}
	}

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
		for (int j = 0; j < n_preds; j++)
			out_f << ",pred_" << to_string(j);
	}

	// Attributes
	if (write_attributes) {
		for (string name : attr_names)
			out_f << ",attr_" << name;
		for (string name : str_attr_names)
			out_f << ",str_attr_" << name;
	}

	// names of features
	for (int j = 0; j < col_names.size(); j++)
		out_f << "," << col_names[j];
	out_f << "\n";

	// data
	write_csv_data(out_f, write_attributes, col_names, 0);

	out_f.close();
	MLOG("Wrote [%zu] rows with %zu features in %s\n", samples.size(), data.size(), csv_fname.c_str());


	return 0;

}

int MedFeatures::add_to_csv_mat(const string &csv_fname, bool write_attributes, int start_idx) const
{
	ofstream out_f;

	// No Sanity check - we assume fitting to previous data, and self-consistency.

	out_f.open(csv_fname, ofstream::out | ofstream::app);

	if (!out_f.is_open()) {
		MERR("ERROR: MedFeatures::write_as_csv_mat() :: Can't open file %s for writing\n", csv_fname.c_str());
		return -1;
	}

	vector<string> col_names;
	get_feature_names(col_names);

	// data
	write_csv_data(out_f, write_attributes, col_names, start_idx);

	out_f.close();
	MLOG("Added [%zu] rows with %zu features in %s\n", samples.size(), data.size(), csv_fname.c_str());


	return 0;
}


void add_to_map(map<string, int> &m, const vector<string> &names) {
	int curr_size = (int)m.size();
	for (int i = 0; i < names.size(); ++i) {
		m[names[i]] = curr_size;
		++curr_size;
	}
}

void op_map(const map<string, int> &m, vector<string> &names) {
	names.resize(m.size() + 1);
	for (auto it = m.begin(); it != m.end(); ++it) {
		if (it->second >= names.size())
			names.resize(it->second + 1);
		names[it->second] = it->first;
	}
}

int max_ind_map(const map<string, int> &m) {
	int max = -1;
	if (!m.empty())
		max = m.begin()->second;
	for (const auto &it : m)
		if (max < it.second)
			max = it.second;
	return max;
}

// Read features (samples + weights + data) from a csv file with a header line
// Return -1 upon failure to open file, 0 upon success
//.......................................................................................
int MedFeatures::read_from_csv_mat(const string &csv_fname, bool read_time_raw)
{
	vector<string> pre_fields = { "serial" }; //fields that appears in all modes: without_weights, with_weights
	vector<string> fields_order = { "id", "time", "outcome", "outcome_time", "split" }; //fields for MedSample
	string weight_field_name = "weight"; //weight field name if appear
	string pred_prefix = "pred_";
	string attr_prefix = "attr_";
	string str_attr_prefix = "str_attr_";

	if (!file_exists(csv_fname))
		MTHROW_AND_ERR("File %s doesn't exist\n", csv_fname.c_str());

	MLOG("reading data from %s\n", csv_fname.c_str());
	ifstream inf;
	inf.open(csv_fname, ios::in);
	if (!inf.good())
		MTHROW_AND_ERR("can not open file for reading %s\n", csv_fname.c_str());

	int ncols = -1;
	string curr_line;
	vector<string> names;

	map<string, int> pos_fields_no_weight;
	add_to_map(pos_fields_no_weight, pre_fields);
	map<string, int> pos_fields_weight = pos_fields_no_weight;
	int sz = (int)pos_fields_weight.size();
	pos_fields_weight[weight_field_name] = sz;
	add_to_map(pos_fields_weight, fields_order);
	add_to_map(pos_fields_no_weight, fields_order);
	vector<string> fields_order_no_weight, fields_order_weight;
	op_map(pos_fields_weight, fields_order_weight);
	op_map(pos_fields_no_weight, fields_order_no_weight);
	vector<string> *curr_fields_order = &fields_order_no_weight;
	map<string, int> *curr_pos_fields = &pos_fields_no_weight;
	map<string, int> pos_attr, pos_str_attr;
	unordered_map<string, int> feature_name_pos;
	vector<int> pos_preds;

	while (getline(inf, curr_line)) {
		boost::trim(curr_line);
		vector<string> fields;
		boost::split(fields, curr_line, boost::is_any_of(","));
		if (ncols == -1) { // Header line	
			int idx = 0;
			vector<int> skiped_input_columns;
			string curr_f;
			for (; idx < pre_fields.size(); ++idx) {
				curr_f = pre_fields[idx];
				if (fields[idx].compare(curr_f) != 0) {
					MLOG("header_line=%s\nIn field %s, idx=(%d / %zu), got_field_header=%s\n",
						curr_line.c_str(), curr_f.c_str(), idx, pre_fields.size(), fields[idx].c_str());
					assert(fields[idx].compare(curr_f) == 0);
				}
			}
			if (fields[idx].compare(fields_order_weight[idx]) == 0) {
				++idx;
				curr_fields_order = &fields_order_weight;
				curr_pos_fields = &pos_fields_weight;
			}
			for (; idx < curr_pos_fields->size(); ++idx)
			{
				curr_f = curr_fields_order->at(idx);
				if (fields[idx + skiped_input_columns.size()].compare(curr_f) != 0) {
					//try also "outcomeTime":
					if (curr_f == "outcome_time" && fields[idx + skiped_input_columns.size()].compare("outcomeTime") == 0)
						continue;

					//search for field in curr_fields_order from idx and above:
					int found_idx = -1;
					for (int s_id = idx + 1; s_id < curr_fields_order->size() && found_idx < 0; ++s_id)
						if (fields[idx + skiped_input_columns.size()].compare(curr_fields_order->at(s_id)) == 0 || (curr_fields_order->at(s_id) == "outcome_time" && fields[idx + skiped_input_columns.size()].compare("outcomeTime") == 0))
							found_idx = s_id;
					//recover and change order:
					if (found_idx >= 0) {
						string found_pos = curr_fields_order->at(found_idx); //what found in fields position (to switch with)
						string curr_pos = curr_fields_order->at(idx); //original expected at idx => will move to found_idx, will look for later

						curr_fields_order->at(found_idx) = curr_pos;
						curr_pos_fields->at(curr_pos) = found_idx;

						curr_fields_order->at(idx) = found_pos; //what fields has currently
						curr_pos_fields->at(found_pos) = idx;

						MLOG("MedFeatures CSV reader :: found %s(should be found in %d) instead %s(%d, input_idx=%d).\n",
							fields[idx + skiped_input_columns.size()].c_str(), found_idx, curr_f.c_str(), idx, idx + skiped_input_columns.size());
					}
					else {
						skiped_input_columns.push_back(idx);
						MWARN("WARN: skipped field %s(%d) in header - saved for later\n", fields[idx].c_str(), idx);
						--idx;
						if (idx + 1 + skiped_input_columns.size() < fields.size())
							continue;
						else
							MTHROW_AND_ERR("In field %s, idx=(%d / %zu), got_field_header=%s, expected=%s. expected_order=[%s]\nheader_line=%s\n",
								curr_f.c_str(), idx, curr_pos_fields->size() - 1, fields[idx].c_str(),
								curr_f.c_str(), medial::io::get_list(*curr_fields_order).c_str(), curr_line.c_str());
					}
				}
			}
			//fetch all skiped_input_columns:
			for (int skip_idx : skiped_input_columns)
			{
				if (boost::starts_with(fields[skip_idx], pred_prefix)) {
					pos_preds.push_back(skip_idx);
					MLOG("Added field %s(%d) into prediction fields\n", fields[skip_idx].c_str(), skip_idx);
					continue;
				}
				if (boost::starts_with(fields[skip_idx], attr_prefix)) {
					string name = fields[skip_idx].substr(attr_prefix.size());
					pos_attr[name] = skip_idx;
					MLOG("Added field %s(%d) into numeric attributes fields\n", fields[skip_idx].c_str(), skip_idx);
					continue;
				}
				if (boost::starts_with(fields[skip_idx], str_attr_prefix)) {
					string name = fields[skip_idx].substr(str_attr_prefix.size());
					pos_str_attr[name] = skip_idx;
					MLOG("Added field %s(%d) into string attributes fields\n", fields[skip_idx].c_str(), skip_idx);
					continue;
				}
				//features:
				data[fields[skip_idx]] = vector<float>();
				attributes[fields[skip_idx]].normalized = attributes[fields[skip_idx]].imputed = false;
				names.push_back(fields[skip_idx]);
				feature_name_pos[fields[skip_idx]] = skip_idx;
				MLOG("Added field %s(%d) into features\n", fields[skip_idx].c_str(), skip_idx);
			}
			idx += (int)skiped_input_columns.size();

			// Predictions
			while (idx < fields.size() && boost::starts_with(fields[idx], pred_prefix)) {
				pos_preds.push_back(idx);
				++idx;
			}

			// Attributes
			while (idx < fields.size() && boost::starts_with(fields[idx], attr_prefix)) {
				string name = fields[idx].substr(attr_prefix.size());
				pos_attr[name] = idx;
				++idx;
			}

			// Str-Attributes
			while (idx < fields.size() && boost::starts_with(fields[idx], str_attr_prefix)) {
				string name = fields[idx].substr(str_attr_prefix.size());
				pos_str_attr[name] = idx;
				++idx;
			}

			// Read Features
			for (int i = idx; i < fields.size(); i++) {
				data[fields[i]] = vector<float>();
				attributes[fields[i]].normalized = attributes[fields[i]].imputed = false;
				names.push_back(fields[i]);
				feature_name_pos[fields[i]] = i;
			}

			ncols = (int)fields.size();
		}
		else { // Data lines
			if (fields.size() != ncols)
				MTHROW_AND_ERR("Expected %d fields, got %d fields in line: \'%s\'\n", ncols, (int)fields.size(), curr_line.c_str());

			if (curr_pos_fields->find(weight_field_name) != curr_pos_fields->end())
				weights.push_back(med_stof(fields[curr_pos_fields->at(weight_field_name)]));

			MedSample newSample;
			newSample.parse_from_string(fields, *curr_pos_fields, pos_preds, pos_attr, pos_str_attr, time_unit, (int)read_time_raw, ",");

			samples.push_back(newSample);

			for (int i = 0; i < names.size(); i++) {
				try {
					data[names[i]].push_back(med_stof(fields[feature_name_pos[names[i]]]));
				}
				catch (...) {
					MERR("Error in line %zu, column %d(%s), value was \"%s\"\n",
						data[names[i]].size() + 2, i, names[i].c_str(),
						fields[feature_name_pos[names[i]]].c_str());
					throw;
				}
			}
		}
	}

	// Check if attribute 'train_weight' exists
	string attr_weight_name = "train_weight";
	if (attributes.find(attr_weight_name) != attributes.end()) {
		if (!weights.empty()) {
			for (size_t i = 0; i < weights.size(); i++)
				if (weights[i] != samples[i].attributes[attr_weight_name])
					MTHROW_AND_ERR("Both weights and attr_train_weight given and are inconsistent. Cannot choose\n")
		}
		else {
			weights.resize(samples.size());
			for (size_t i = 0; i < weights.size(); i++)
				weights[i] = samples[i].attributes[attr_weight_name];
		}
	}

	inf.close();
	return 0;
}


// Filter data (and attributes) to include only selected features
// Return -1 if any of the selected features is not present. 0 upon success.
//.......................................................................................
int MedFeatures::filter(unordered_set<string>& selectedFeatures) {

	// Sanity
	for (string feature : selectedFeatures) {
		if (data.find(feature) == data.end()) {
			MERR("Error in MedFeatures::filter - Cannot find feature %s in Matrix\n", feature.c_str());
			vector<string> all_names;
			get_feature_names(all_names);
			string all_opts = medial::io::get_list(all_names, "\n");
			MERR("All Feature Options(%zu):\n%s\n", all_names.size(), all_opts.c_str());
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
		tags.erase(feature);
	}

	return 0;
}

// Get the corresponding MedSamples object .  Assuming samples are ordered in features (all id's samples are consecutive)
//.......................................................................................
void MedFeatures::get_samples(MedSamples& outSamples) const {

	for (auto& sample : samples) {
		if (outSamples.idSamples.size() && outSamples.idSamples.back().id == sample.id)
			outSamples.idSamples.back().samples.push_back(sample);
		else {
			MedIdSamples newIdSample;
			newIdSample.id = sample.id;
			newIdSample.split = sample.split;
			newIdSample.samples.push_back(sample);
			outSamples.idSamples.push_back(newIdSample);
		}
	}

}

// Find the max serial_id_cnt
//.......................................................................................
int MedFeatures::get_max_serial_id_cnt() const {

	int max = 0;
	for (auto& rec : data) {
		string name = rec.first;
		if (name.substr(0, 4) == "FTR_") {
			int n = stoi(name.substr(4, name.length()));
			if (n > max)
				max = n;
		}
	}

	return max;
}

//................................................................................................
int MedFeatures::prep_selected_list(vector<string>& search_str, unordered_set<string> &selected)
{
	for (auto &f : data) {
		for (auto &s : search_str)
			if (s != "" && f.first.find(s) != string::npos)
				selected.insert(f.first);
	}

	return 0;
}

template<class T> void medial::process::commit_selection(vector<T> &vec, const vector<int> &idx) {
	vector<T> filt(idx.size());
	for (size_t i = 0; i < idx.size(); ++i)
		filt[i] = vec[idx[i]];
	vec.swap(filt);
}
template void medial::process::commit_selection<MedSample *>(vector<MedSample *> &vec, const vector<int> &idx);
template void medial::process::commit_selection<const MedSample *>(vector<const MedSample *> &vec, const vector<int> &idx);
template void medial::process::commit_selection<float>(vector<float> &vec, const vector<int> &idx);
template void medial::process::commit_selection<double>(vector<double> &vec, const vector<int> &idx);
template void medial::process::commit_selection<int>(vector<int> &vec, const vector<int> &idx);

//Assume selected_indexes is "sorted"
void medial::process::filter_row_indexes_safe(MedFeatures &dataMat, const vector<int> &selected_indexes, bool op_flag) {
	MedFeatures filtered;
	filtered.time_unit = dataMat.time_unit;
	filtered.attributes = dataMat.attributes;

	int curr_ind = 0;
	if (!op_flag) {
		for (auto iit = dataMat.data.begin(); iit != dataMat.data.end(); ++iit)
			filtered.data[iit->first].reserve(selected_indexes.size());
		if (!dataMat.weights.empty())
			filtered.weights.reserve(selected_indexes.size());
		if (!dataMat.masks.empty()) {
			for (auto iit = dataMat.masks.begin(); iit != dataMat.masks.end(); ++iit)
				filtered.masks[iit->first].reserve(selected_indexes.size());
		}

		filtered.samples.reserve(selected_indexes.size());
		for (int i : selected_indexes) //all selected indexes
		{
			filtered.samples.push_back(dataMat.samples[i]);
			for (auto iit = dataMat.data.begin(); iit != dataMat.data.end(); ++iit)
				filtered.data[iit->first].push_back(iit->second[i]);

			if (!dataMat.weights.empty())
				filtered.weights.push_back(dataMat.weights[i]);
			if (!dataMat.masks.empty()) {
				for (auto iit = dataMat.masks.begin(); iit != dataMat.masks.end(); ++iit)
					filtered.masks[iit->first].push_back(iit->second[i]);
			}

		}
	}
	else {
		for (auto iit = dataMat.data.begin(); iit != dataMat.data.end(); ++iit)
			filtered.data[iit->first].reserve((int)dataMat.samples.size() - (int)selected_indexes.size());
		if (!dataMat.weights.empty())
			filtered.weights.reserve((int)dataMat.samples.size() - (int)selected_indexes.size());
		if (!dataMat.masks.empty()) {
			for (auto iit = dataMat.masks.begin(); iit != dataMat.masks.end(); ++iit)
				filtered.masks[iit->first].reserve(selected_indexes.size());
		}
		filtered.samples.reserve((int)dataMat.samples.size() - (int)selected_indexes.size());
		for (int i = 0; i < dataMat.samples.size(); ++i)
		{
			//remove selected row when matched:
			if (curr_ind < selected_indexes.size() && i == selected_indexes[curr_ind]) {
				++curr_ind;
				continue;
			}
			filtered.samples.push_back(dataMat.samples[i]);
			for (auto iit = dataMat.data.begin(); iit != dataMat.data.end(); ++iit)
				filtered.data[iit->first].push_back(iit->second[i]);
			if (!dataMat.weights.empty())
				filtered.weights.push_back(dataMat.weights[i]);
			if (!dataMat.masks.empty()) {
				for (auto iit = dataMat.masks.begin(); iit != dataMat.masks.end(); ++iit)
					filtered.masks[iit->first].push_back(iit->second[i]);
			}
		}
	}
	filtered.init_pid_pos_len();

	//dataMat = filtered;
	dataMat.samples.swap(filtered.samples);
	dataMat.data.swap(filtered.data);
	dataMat.weights.swap(filtered.weights);
	dataMat.pid_pos_len.swap(filtered.pid_pos_len);
	dataMat.attributes.swap(filtered.attributes);
	dataMat.tags.swap(filtered.tags);
	dataMat.time_unit = filtered.time_unit;
	dataMat.masks.swap(filtered.masks);
}

void medial::process::filter_row_indexes(MedFeatures &dataMat, vector<int> &selected_indexes, bool op_flag) {
	sort(selected_indexes.begin(), selected_indexes.end());
	filter_row_indexes_safe(dataMat, selected_indexes, op_flag);
}

void medial::process::down_sample(MedFeatures &dataMat, double take_ratio, bool with_repeats,
	vector<int> *selected_indexes) {
	int final_cnt = int(take_ratio * dataMat.samples.size());
	if (take_ratio >= 1) {
		return;
	}
	vector<int> all_selected_indexes;
	if (selected_indexes == NULL)
		selected_indexes = &all_selected_indexes;
	selected_indexes->resize(final_cnt);

	vector<bool> seen_index((int)dataMat.samples.size());
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<> dist_gen(0, (int)dataMat.samples.size() - 1);
	for (size_t k = 0; k < final_cnt; ++k) //for 0 and 1:
	{
		int num_ind = dist_gen(gen);
		if (!with_repeats) {
			while (seen_index[num_ind])
				num_ind = dist_gen(gen);
			seen_index[num_ind] = true;
		}
		(*selected_indexes)[k] = num_ind;
	}
	filter_row_indexes(dataMat, *selected_indexes);
}

double medial::process::reweight_by_general(MedFeatures &data_records, const vector<string> &groups,
	vector<float> &weigths, bool print_verbose) {
	if (groups.size() != data_records.samples.size())
		MTHROW_AND_ERR("data_records and groups should hsve same size\n");

	vector<float> full_weights(groups.size(), 1);
	vector<unordered_map<string, vector<int>>> list_label_groups(2);
	vector<unordered_map<string, int>> count_label_groups(2);
	vector<string> all_groups;
	unordered_set<string> seen_pid_0;
	unordered_set<string> seen_group;
	unordered_map<string, unordered_set<int>> group_to_seen_pid; //of_predicition

	for (size_t i = 0; i < data_records.samples.size(); ++i)
	{
		//int year = int(year_bin_size * round((it->registry.date / 10000) / year_bin_size));
		int label = int(data_records.samples[i].outcome > 0);

		if ((label > 0) && seen_group.find(groups[i]) == seen_group.end()) {
			all_groups.push_back(groups[i]);
			seen_group.insert(groups[i]);
		}

		list_label_groups[label][groups[i]].push_back((int)i);
		++count_label_groups[label][groups[i]];
		if (label == 0) {
			group_to_seen_pid[groups[i]].insert(data_records.samples[i].id);
		}
	}

	unordered_map<string, int> year_total;
	unordered_map<string, float> year_ratio;
	int i = 0;
	sort(all_groups.begin(), all_groups.end());
	if (print_verbose) {
		MLOG("Before Matching Total samples is %d on %d groups\n",
			(int)data_records.samples.size(), (int)all_groups.size());
		MLOG("Group" "\t" "Count_0" "\t" "Count_1" "\t" "ratio\n");
	}
	for (const string &grp : all_groups)
	{
		if (count_label_groups[0][grp] == 0)
			continue;
		year_total[grp] = count_label_groups[0][grp] + count_label_groups[1][grp];
		year_ratio[grp] = count_label_groups[1][grp] / float(count_label_groups[0][grp] + count_label_groups[1][grp]);
		++i;

		if (print_verbose)
			MLOG("%s\t%d\t%d\t%2.5f\n", grp.c_str(), count_label_groups[0][grp], count_label_groups[1][grp],
				count_label_groups[1][grp] / float(count_label_groups[1][grp] + count_label_groups[0][grp]));
	}

	float r_target = 0.5;
	double max_factor = 0;

	//For each year_bin - balance to this ratio using price_ratio weight for removing 1's labels:
	seen_pid_0.clear();
	unordered_map<string, float> group_to_factor;
	for (int k = int(all_groups.size() - 1); k >= 0; --k) {
		string &grp = all_groups[k];

		float base_ratio = count_label_groups[1][grp] / float(count_label_groups[1][grp] + count_label_groups[0][grp]);
		float factor = 0;
		if (base_ratio > 0)
			factor = float((r_target / (1 - r_target)) *
			(double(count_label_groups[0][grp]) / count_label_groups[1][grp]));
		if (factor > max_factor)
			max_factor = factor;


		if (factor > 0) {
			if (print_verbose)
				MLOG("Weighting group %s base_ratio=%f factor= %f\n",
					grp.c_str(), base_ratio, factor);
		}
		else
			MLOG("Dropping group %s Num_controls=%d with %d cases\n",
				grp.c_str(), count_label_groups[0][grp], count_label_groups[1][grp]);

		if (factor <= 0) {
			//list_label_groups[0].erase(grp); //for zero it's different
			//list_label_groups[1].erase(grp); //for zero it's different
			count_label_groups[0][grp] = 0;
			count_label_groups[1][grp] = 0;
			//set weights 0 for all:
			for (int ind : list_label_groups[1][grp])
				full_weights[ind] = 0;
			for (int ind : list_label_groups[0][grp])
				full_weights[ind] = 0;
		}
		else {
			group_to_factor[grp] = factor;
			//Commit Change to weights:
			for (int ind : list_label_groups[1][grp])
				full_weights[ind] = factor;
		}
	}
	//let's divide all by 1/max_factor:
	double total_cnt_after = 0, init_sample_count = 0;
	for (auto it = group_to_factor.begin(); it != group_to_factor.end(); ++it) {
		total_cnt_after += count_label_groups[0][it->first] + (double)count_label_groups[1][it->first] * it->second;
		init_sample_count += (count_label_groups[0][it->first] + count_label_groups[1][it->first]);
	}
	max_factor = float(total_cnt_after / init_sample_count); //change to take not max - keep same sum
	MLOG_D("init_sample_count=%d, total_cnt_after=%2.1f, max_factor=%2.1f, orig_size=%d"
		",all_groups.size()=%d, group_to_factor.size()=%d\n",
		(int)init_sample_count, (float)total_cnt_after, (float)max_factor,
		(int)data_records.samples.size(), (int)all_groups.size(), (int)group_to_factor.size());

	if (max_factor > 0) {
		for (size_t i = 0; i < full_weights.size(); ++i)
			//if (data_records.samples[i].outcome > 0)
			full_weights[i] /= (float)max_factor;
		for (auto it = group_to_factor.begin(); it != group_to_factor.end(); ++it)
			group_to_factor[it->first] /= (float)max_factor;
	}

	//Commit on all records:
	MedFeatures filtered;
	filtered.time_unit = data_records.time_unit;
	filtered.attributes = data_records.attributes;

	vector<int> all_selected_indexes;
	for (size_t k = 0; k < list_label_groups.size(); ++k) //for 0 and 1:
	{
		for (auto it = list_label_groups[k].begin(); it != list_label_groups[k].end(); ++it) //for each group
		{
			vector<int> &ind_list = it->second;
			all_selected_indexes.insert(all_selected_indexes.end(), ind_list.begin(), ind_list.end());
		}
	}

	sort(all_selected_indexes.begin(), all_selected_indexes.end());
	filter_row_indexes(data_records, all_selected_indexes);
	weigths.clear();
	weigths.resize(all_selected_indexes.size(), 1);
	//filter full_weights to weights:
	for (size_t i = 0; i < all_selected_indexes.size(); ++i)
		weigths[i] = full_weights[all_selected_indexes[i]];

	if (print_verbose) {
		MLOG("After Matching Size=%d:\n", (int)data_records.samples.size());
		unordered_map<string, vector<int>> counts_stat;
		unordered_map<string, vector<double>> weight_stats;
		for (size_t i = 0; i < all_selected_indexes.size(); ++i)
		{
			if (counts_stat[groups[all_selected_indexes[i]]].empty()) {
				counts_stat[groups[all_selected_indexes[i]]].resize(2);
				weight_stats[groups[all_selected_indexes[i]]].resize(2);
			}
			++counts_stat[groups[all_selected_indexes[i]]][data_records.samples[i].outcome > 0];
			weight_stats[groups[all_selected_indexes[i]]][data_records.samples[i].outcome > 0] +=
				weigths[i];
		}
		MLOG("Group\tCount_0\tCount_1\tratio\tweight_cases\tweighted_ratio\n");
		for (const string &grp : all_groups)
			if (group_to_factor.find(grp) != group_to_factor.end())
				MLOG("%s\t%d\t%d\t%2.5f\t%2.5f\t%2.5f\n", grp.c_str(),
					counts_stat[grp][0], counts_stat[grp][1],
					counts_stat[grp][1] / double(counts_stat[grp][1] + counts_stat[grp][0]),
					group_to_factor.at(grp),
					weight_stats[grp][1] / double(weight_stats[grp][1] + weight_stats[grp][0]));
		//print_by_year(data_records.samples);
	}
	if (max_factor > 0)
		return 1.0 / max_factor;
	else
		return (double)0;
}

void  medial::process::match_by_general(MedFeatures &data_records, const vector<string> &groups, vector<int> &filtered_row_ids, float price_ratio, int min_grp_size, bool print_verbose) {
	medial::process::match_by_general(data_records, groups, filtered_row_ids, price_ratio, -1.0, min_grp_size, print_verbose);
}

void  medial::process::match_by_general(MedFeatures &data_records, const vector<string> &groups, vector<int> &filtered_row_ids, float price_ratio, float max_ratio, int min_grp_size, bool print_verbose) {
	if (groups.size() != data_records.samples.size())
		MTHROW_AND_ERR("data_records and groups should hsve same size\n");

	vector<unordered_map<string, vector<int>>> list_label_groups(2);
	vector<unordered_map<string, int>> count_label_groups(2);
	vector<string> all_groups;
	unordered_set<string> seen_pid_0;
	unordered_set<string> seen_group;
	unordered_map<string, unordered_set<int>> group_to_seen_pid; //of_predicition

	for (size_t i = 0; i < data_records.samples.size(); ++i)
	{
		//int year = int(year_bin_size * round((it->registry.date / 10000) / year_bin_size));
		int label = int(data_records.samples[i].outcome > 0);

		if ((label > 0) && seen_group.find(groups[i]) == seen_group.end()) {
			all_groups.push_back(groups[i]);
			seen_group.insert(groups[i]);
		}

		list_label_groups[label][groups[i]].push_back((int)i);
		++count_label_groups[label][groups[i]];
		if (label == 0) {
			group_to_seen_pid[groups[i]].insert(data_records.samples[i].id);
		}
	}
	//remove groups with only controls:
	for (auto it = list_label_groups[0].begin(); it != list_label_groups[0].end(); ++it)
		if (seen_group.find(it->first) == seen_group.end()) {
			MWARN("Warning: group %s has only %d controls with no cases- skipping\n",
				it->first.c_str(), (int)it->second.size());
			list_label_groups[0][it->first].clear();
		}

	unordered_map<string, int> year_total;
	unordered_map<string, float> year_ratio;
	vector<float> all_ratios((int)all_groups.size());
	int i = 0;
	sort(all_groups.begin(), all_groups.end());
	if (print_verbose) {
		MLOG("Before Matching Total samples is %d on %d groups\n",
			(int)data_records.samples.size(), (int)all_groups.size());
		MLOG("Group"  "\t"  "Count_0"  "\t"  "Count_1"  "\t"  "ratio" "\t" "required_price_ratio" "\n");
	}
	for (const string &grp : all_groups)
	{
		if (count_label_groups[0][grp] == 0)
			continue;
		year_total[grp] = count_label_groups[0][grp] + count_label_groups[1][grp];
		year_ratio[grp] = count_label_groups[1][grp] / float(count_label_groups[0][grp] + count_label_groups[1][grp]);
		all_ratios[i] = year_ratio[grp];
		++i;
	}
	//Choose ratio to balance all for:
	sort(all_ratios.begin(), all_ratios.end());
	if (print_verbose) {
		vector<float> controls_sum(all_groups.size()), cases_sum(all_groups.size());
		for (const string &grp : all_groups) {
			float grp_ratio = year_ratio[grp];
			int ratio_ind = medial::process::binary_search_index(all_ratios.data(),
				all_ratios.data() + all_ratios.size() - 1, grp_ratio);
			if (ratio_ind < 0) {
				MWARN("warning: bug in binary search - matching(effects just verbose printing)\n");
				break;
			}
			for (const string &grp_calc : all_groups) {
				float grp_ratio_clc = count_label_groups[1][grp_calc] / float(count_label_groups[1][grp_calc] + count_label_groups[0][grp_calc]);
				if (grp_ratio_clc < grp_ratio)
					controls_sum[ratio_ind] += (count_label_groups[1][grp_calc] + count_label_groups[0][grp_calc]) -
					count_label_groups[1][grp_calc] / grp_ratio;
				else
					cases_sum[ratio_ind] += (count_label_groups[1][grp_calc] -
					(count_label_groups[1][grp_calc] + count_label_groups[0][grp_calc])*grp_ratio) /
						(1 - grp_ratio);
			}
		}
		for (const string &grp : all_groups) {
			float grp_ratio = year_ratio[grp];
			int ratio_ind = medial::process::binary_search_index(all_ratios.data(),
				all_ratios.data() + all_ratios.size() - 1, grp_ratio);

			if (ratio_ind < 0)
				break;

			float factor_needed_down = 0, factor_needed_up = -1;
			if (ratio_ind < all_ratios.size() - 1) {
				float diff_cases = abs(cases_sum[ratio_ind] - cases_sum[ratio_ind + 1]);
				float diff_controls = abs(controls_sum[ratio_ind] - controls_sum[ratio_ind + 1]);
				if (diff_cases == 0)
					factor_needed_up = -1;
				else
					factor_needed_up = diff_controls / diff_cases;
			}

			if (ratio_ind > 0) {
				float diff_cases = abs(cases_sum[ratio_ind] - cases_sum[ratio_ind - 1]);
				float diff_controls = abs(controls_sum[ratio_ind] - controls_sum[ratio_ind - 1]);
				if (diff_cases == 0)
					factor_needed_down = -1;
				else
					factor_needed_down = diff_controls / diff_cases;
			}

			if (count_label_groups[0][grp] == 0)
				grp_ratio = 1; //just for correct printing
			if (factor_needed_up > 0) {
				if (factor_needed_up < factor_needed_down)
					MLOG("%s\t%d\t%d\t%f\t[NOT_AN_OPTION]\n", grp.c_str(), count_label_groups[0][grp], count_label_groups[1][grp]
						, grp_ratio);
				else
					MLOG("%s\t%d\t%d\t%f\t[%2.2f-%2.2f]\n", grp.c_str(), count_label_groups[0][grp], count_label_groups[1][grp]
						, grp_ratio, factor_needed_down, factor_needed_up);
			}
			else {
				MLOG("%s\t%d\t%d\t%f\t[%2.2f-]\n", grp.c_str(), count_label_groups[0][grp], count_label_groups[1][grp]
					, grp_ratio, factor_needed_down);
			}
		}
	}

	float r_target = 0;
	float best_cost = -1;
	int best_0_rem = 0, best_1_rem = 0;

	for (size_t k = 0; k < all_ratios.size(); ++k)
	{
		//evaluate if this was choosen:
		float curr_target = all_ratios[k];
		if (curr_target == 0)
			continue;
		float cost = 0;
		int tot_0_rem = 0, tot_1_rem = 0;
		for (auto it = count_label_groups[1].begin(); it != count_label_groups[1].end(); ++it) {
			float cost_val = 0;
			if (year_ratio[it->first] > curr_target) { //remove 1's (too much 1's):
				float shrink_factor = (curr_target * count_label_groups[0][it->first]) /
					(count_label_groups[1][it->first] - (count_label_groups[1][it->first] * curr_target)); //to multiply by 1's
				cost_val = (count_label_groups[1][it->first] - int(round(shrink_factor*count_label_groups[1][it->first]))) * price_ratio;
				tot_1_rem += count_label_groups[1][it->first] - int(round(shrink_factor*count_label_groups[1][it->first]));
			}
			else {
				float shrink_factor = 1;
				if (count_label_groups[0][it->first] > 0)
					shrink_factor = (1 - curr_target) * count_label_groups[1][it->first] /
					(curr_target * count_label_groups[0][it->first]); //to multiply by 0's
				cost_val = (float)(count_label_groups[0][it->first] - int(round(shrink_factor*count_label_groups[0][it->first])));
				tot_0_rem += count_label_groups[0][it->first] - int(round(shrink_factor*count_label_groups[0][it->first]));
			}
			cost += cost_val;
		}

		if (print_verbose)
			MLOG("Sampling ratio = %2.3f - Cost = %2.3f removing [%d,%d]\n", curr_target, cost, tot_0_rem, tot_1_rem);

		if (best_cost == -1 || cost < best_cost) {
			best_cost = cost;
			r_target = curr_target;
			best_0_rem = tot_0_rem;
			best_1_rem = tot_1_rem;
		}
	}
	//r_target = prctil(all_ratios, 0.5);
	if (!print_verbose)
		MLOG_D("Best Target is %2.3f so retargeting balance to it. cost=%2.3f remove [%d, %d]\n", r_target, best_cost, best_0_rem
			, best_1_rem);
	else
		MLOG("Best Target is %2.3f so retargeting balance to it. cost=%2.3f remove [%d, %d]\n", r_target, best_cost, best_0_rem
			, best_1_rem);

	if (max_ratio > 0) {
		float min_r = 1 / (1.0 + max_ratio);
		if (r_target < min_r)
			r_target = min_r;
		if (!print_verbose)
			MLOG_D("Correcting target to %2.3f due to max_ratio = %2.3f\n", r_target, max_ratio);
		else
			MLOG("Correcting target to %2.3f due to max_ratio = %2.3f\n", r_target, max_ratio);
	}

	//For each year_bin - balance to this ratio using price_ratio weight for removing 1's labels:
	seen_pid_0.clear();
	vector<int> skip_grp_indexs;
	for (int k = int(all_groups.size() - 1); k >= 0; --k) {
		string &grp = all_groups[k];
		int target_size = 0;
		int remove_size = 0;
		int ind = 0;
		float shrink_factor = 1;
		if (year_ratio[grp] > r_target) {
			shrink_factor = (r_target * count_label_groups[0][grp]) / (count_label_groups[1][grp] - (count_label_groups[1][grp] * r_target)); //to multiply by 1's
			ind = 1;
		}
		else {
			if (count_label_groups[0][grp] > 0)
				shrink_factor = (1 - r_target) * count_label_groups[1][grp] / (r_target * count_label_groups[0][grp]); //to multiply by 0's
			else
				shrink_factor = 1;
		}
		target_size = int(round(shrink_factor*count_label_groups[ind][grp]));
		remove_size = count_label_groups[ind][grp] - target_size;

		if (print_verbose)
			cout << "Doing group " << grp << " ind=" << ind << " target_size=" << target_size
			<< " removing= " << remove_size << endl;
		if (count_label_groups[0][grp] < min_grp_size || count_label_groups[1][grp] < min_grp_size) {
			MWARN("Warning: matching group has very small counts - skipping group=%s [%d, %d]\n",
				grp.c_str(), count_label_groups[0][grp], count_label_groups[1][grp]);
			list_label_groups[0][grp].clear();
			list_label_groups[1][grp].clear();
			skip_grp_indexs.push_back(k);
			continue;
		}
		unordered_set<int> seen_year_pid;
		shuffle(list_label_groups[ind][grp].begin(), list_label_groups[ind][grp].end(), globalRNG::get_engine());
		if (target_size > list_label_groups[ind][grp].size())
			MERR("ERROR/BUG: try to shrink %d into %d\n"
				, (int)list_label_groups[ind][grp].size(), target_size);
		else
			list_label_groups[ind][grp].resize(target_size); //for zero it's different
	}

	for (int k = 0; k < skip_grp_indexs.size(); ++k)
		all_groups.erase(all_groups.begin() + skip_grp_indexs[k]);

	//Commit on all records:
	for (size_t k = 0; k < list_label_groups.size(); ++k) //for 0 and 1:
		for (auto it = list_label_groups[k].begin(); it != list_label_groups[k].end(); ++it) //for each year
		{
			vector<int> &ind_list = it->second;
			filtered_row_ids.insert(filtered_row_ids.end(), ind_list.begin(), ind_list.end());
		}


	filter_row_indexes(data_records, filtered_row_ids);

	if (print_verbose) {
		MLOG("After Matching Size=%d:\n", (int)data_records.samples.size());
		print_stats(data_records, filtered_row_ids, groups);
	}
}

// Helper functions for multi-class matching
#define MATCHING_EPS 0.0001
void get_sampling_ratios(vector<float>& targetRatios, vector<vector<float>>& ratios, vector<vector<float>>& samplingRatios) {

	int nClasses = (int)targetRatios.size();
	int nGroups = (int)ratios.size();
	samplingRatios.resize(nGroups);

	for (int i = 0; i < nGroups; i++) {
		vector<float> ratioOfRatios(nClasses);
		for (int j = 0; j < nClasses; j++)
			ratioOfRatios[j] = ratios[i][j] / targetRatios[j];

		// Minimum, non-zerio r-of-r
		float minROR = 0;
		for (int j = 0; j < nClasses; j++) {
			if (ratioOfRatios[j] > 0 && (minROR == 0 || minROR > ratioOfRatios[j]))
				minROR = ratioOfRatios[j];
		}

		// Sampling ratios
		samplingRatios[i].assign(nClasses, 1.0);
		for (int j = 0; j < nClasses; j++) {
			if (ratioOfRatios[j] > 0)
				samplingRatios[i][j] = minROR / ratioOfRatios[j];
		}
	}
}
float get_multi_class_matching_loss(vector<float>& targetRatios, vector<vector<float>>& ratios, vector<vector<int>>& counts, vector<float>& price_ratios, int verbose) {

	vector<vector<float>> samplingRatios;
	get_sampling_ratios(targetRatios, ratios, samplingRatios);

	float loss = 0;
	int nClasses = (int)targetRatios.size();
	int nGroups = (int)ratios.size();

	for (int i = 0; i < nGroups; i++) {
		for (int j = 0; j < nClasses; j++)
			loss += price_ratios[j] * counts[i][j] * (1 - samplingRatios[i][j]);
	}

	return loss;
}

// Check step size for increasing targetRatios[i] and decrasing targetRatios[j] until minROR is changed.
float get_step_size(vector<vector<float>>& ratios, vector<float>& targetRatios, vector<float>& minROR, vector<int>& minROR_Idx, int i, int j) {


	float stepSize = 1.0 - MATCHING_EPS - targetRatios[i];
	if (targetRatios[j] - (0.0 + MATCHING_EPS) < stepSize)
		stepSize = targetRatios[j] - (0.0 + MATCHING_EPS);

	for (size_t g = 0; g < minROR.size(); g++) {
		float _minROR = 0.0;
		int _minROR_Idx = -1;
		if (minROR_Idx[g] == j) {
			// Non-i classes - we can decrease j until we get ROR[j] == ROR[c]
			for (size_t c = 0; c < targetRatios.size(); c++) {
				if (c != i && c != j) {
					float testStep = (ratios[g][c] * targetRatios[j] - ratios[g][j] * targetRatios[c]) / ratios[g][c];
					if (testStep < stepSize) {
						_minROR = ratios[g][c] / targetRatios[c];
						_minROR_Idx = (int)c;
						stepSize = testStep;
					}
				}
			}
			// i - we increase i and decrease j simultanously until  we get ROR[j] == ROR[j]
			float testStep = (ratios[g][i] * targetRatios[j] - ratios[g][j] * targetRatios[i]) / (ratios[g][i] + ratios[g][j]);
			if (testStep < stepSize) {
				_minROR = ratios[g][i] / (targetRatios[i] + testStep);
				_minROR_Idx = i;
				stepSize = testStep;
			}

			if (_minROR_Idx != -1) {
				minROR[g] = _minROR;
				minROR_Idx[g] = _minROR_Idx;
			}
		}
		else if (minROR_Idx[g] != i) {
			// We can increase i until we get ROR[i] = ROR[c]
			float testStep = ratios[g][i] / minROR[g] - targetRatios[i];
			if (testStep < stepSize) {
				minROR_Idx[g] = i;
				stepSize = testStep;
			}
		}
		// If minROR_Idx[j] = j, and we decrease targetRatios[j], there is not boundary.
	}

	return stepSize;
}
float do_greedy_search(vector<float>& targetRatios, vector<vector<float>>& ratios, vector<vector<int>>& counts, vector<float>& price_ratios, int verbose) {

	int nClasses = (int)targetRatios.size();
	int nGroups = (int)ratios.size();

	// minimal Ratio of Ratios per group 
	vector<float> minROR(nGroups);
	vector<int> minROR_Idx(nGroups);
	for (int i = 0; i < nGroups; i++) {
		minROR[i] = ratios[i][0] / targetRatios[0];
		minROR_Idx[i] = 0;
		for (int j = 1; j < nClasses; j++) {
			if (ratios[i][j] / targetRatios[j] < minROR[i]) {
				minROR[i] = ratios[i][j] / targetRatios[j];
				minROR_Idx[i] = j;
			}
		}
	}

	float loss = get_multi_class_matching_loss(targetRatios, ratios, counts, price_ratios, verbose);

	// Search
	bool keepGoing = true;
	int nSteps = 0;
	while (keepGoing) {
		keepGoing = false;

		// Verbosity
		if (verbose) {
			MLOG("minROR");
			for (int i = 0; i < nGroups; i++)
				MLOG("\t%.2f", minROR[i]);
			MLOG("\n");
			MLOG("Idx");
			for (int i = 0; i < nGroups; i++)
				MLOG("\t%d", minROR_Idx[i]);
			MLOG("\n");
			MLOG("TargetR");
			for (int i = 0; i < nClasses; i++)
				MLOG("\t%.4f", targetRatios[i]);
			MLOG("\n");
		}

		// Check all neighbours
		for (int i = 0; i < nClasses; i++) {
			for (int j = 0; j < nClasses; j++) {
				if (i != j) {
					nSteps++;

					vector<float> origMinROR = minROR;
					vector<int> origMinROR_Idx = minROR_Idx;
					float stepSize = get_step_size(ratios, targetRatios, minROR, minROR_Idx, i, j);

					if (stepSize > 0.001) {
						targetRatios[i] += stepSize;
						targetRatios[j] -= stepSize;
						float newLoss = get_multi_class_matching_loss(targetRatios, ratios, counts, price_ratios, verbose);

						if (verbose)
							MLOG("Multi-Class Matching : Loss for %d/%d +- %f => (%f/%f) = %f", i, j, stepSize, targetRatios[i], targetRatios[j], newLoss);

						if (newLoss < loss) {
							keepGoing = true;
							loss = newLoss;
							if (verbose)
								MLOG("  -- Found new Min Loss\n");
							break;
						}
						if (verbose)
							MLOG("\n");

						targetRatios[i] -= stepSize;
						targetRatios[j] += stepSize;
					}
					minROR = origMinROR;
					minROR_Idx = origMinROR_Idx;
				}
			}
			if (keepGoing)
				break;
		}
	}

	if (verbose)
		MLOG("Multi-Class Matching : Number of steps to (local) minimum = %d. Loss = %f\n", nSteps, loss);
	return loss;
}

int prepare_for_matching(vector<MedSample>& samples, const vector<string>& groups, vector<int>& class_idx, vector<string>& groups_v, vector<vector<float>>& ratios, vector<vector<int>>& counts, int verbose) {

	// Collect classes
	set<int> classes_set;
	for (MedSample& sample : samples)
		classes_set.insert((int)sample.outcome);
	int nClasses = (int)classes_set.size();


	int maxClass = 0;
	for (int _class : classes_set) {
		if (_class > maxClass)
			maxClass = _class;
	}

	class_idx.resize(maxClass + 1);
	int idx = 0;
	for (int _class : classes_set)
		class_idx[_class] = idx++;

	// Collect group counts
	map<string, vector<int>> group_class_counts;
	map<string, int> group_tot_counts;
	set<string> all_groups;
	for (size_t i = 0; i < samples.size(); ++i)
	{
		int iClass = class_idx[(int)samples[i].outcome];
		if (all_groups.find(groups[i]) == all_groups.end())
			group_class_counts[groups[i]].assign(nClasses, 0);
		group_tot_counts[groups[i]] ++;
		all_groups.insert(groups[i]);
		group_class_counts[groups[i]][iClass]++;
	}

	// Get ratios and counts
	int nGroups = (int)all_groups.size();
	groups_v.clear();
	groups_v.insert(groups_v.end(), all_groups.begin(), all_groups.end());
	ratios.resize(nGroups, vector<float>(nClasses));
	counts.resize(nGroups, vector<int>(nClasses));
	for (int i = 0; i < nGroups; i++) {
		for (int j = 0; j < nClasses; j++) {
			float ratio = (group_class_counts[groups_v[i]][j] + 0.0) / group_tot_counts[groups_v[i]];
			if (ratio < MATCHING_EPS)
				ratio = (float)MATCHING_EPS;
			ratios[i][j] = ratio;

			counts[i][j] = group_class_counts[groups_v[i]][j];
		}
	}

	// Verbosity
	if (verbose) {
		MLOG("GROUP");
		for (int i = 0; i < nClasses; i++)
			MLOG("\tCLS_%d", i);
		MLOG("\n");
		for (int i = 0; i < nGroups; i++) {
			MLOG("%s", groups_v[i].c_str());
			for (int j = 0; j < nClasses; j++)
				MLOG("\t%d", counts[i][j]);
			MLOG("\n");
		}
	}

	return nClasses;
}

float get_matching_dist(vector<MedSample>& samples, const vector<string> &groups, vector<float>& targetRatios, vector<float>& price_ratios, int nRand, int verbose,
	vector<vector<float>>& ratios, vector<vector<int>>& counts, vector<string>& groups_v, vector<int>& class_idx) {

	// Sanity
	if (groups.size() != samples.size())
		MTHROW_AND_ERR("data samples and groups should hsve same size\n");

	int nClasses = prepare_for_matching(samples, groups, class_idx, groups_v, ratios, counts, verbose);

	if (nClasses != price_ratios.size())
		MTHROW_AND_ERR("price_ratio and number of classes are not compatible\n");
	int nGroups = (int)ratios.size();

	// verbosity
	if (verbose) {
		MLOG("P.Ratio");
		for (int j = 0; j < nClasses; j++)
			MLOG("\t%.3f", price_ratios[j]);
		MLOG("\n");
	}

	// Find optimal ratios - search on the discrete steps
	// Start with random ratio
	targetRatios.resize(nClasses);
	float sum, loss;

	vector<pair<float, float>> rRanges(nClasses);
	for (int i = 0; i < nClasses; i++) {
		rRanges[i] = { ratios[0][i],ratios[0][i] };
		for (int j = 1; j < nGroups; j++) {
			if (ratios[j][i] < rRanges[i].first)
				rRanges[i].first = ratios[j][i];
			if (ratios[j][i] > rRanges[i].second)
				rRanges[i].second = ratios[j][i];
		}
	}

	// Do multiple greedy searches from random starting points

	// Generate random intial vectors to have consistency when threading
	vector<vector<float>> randRatios(nRand, vector<float>(nClasses));
	vector<float> randLosses(nRand);

	for (int i = 0; i < nRand; i++) {
		sum = 0;
		for (int j = 0; j < nClasses; j++) {
			randRatios[i][j] = rRanges[j].first + (rRanges[j].second - rRanges[j].first)*(globalRNG::rand() / (globalRNG::max() + 1.0));
			sum += randRatios[i][j];
		}
		for (int j = 0; j < nClasses; j++)
			randRatios[i][j] /= sum;
	}


#pragma omp parallel for
	for (int i = 0; i < nRand; i++)
		randLosses[i] = do_greedy_search(randRatios[i], ratios, counts, price_ratios, 0);

	int idx = 0;
	for (int i = 1; i < nRand; i++) {
		if (randLosses[i] < randLosses[idx])
			idx = i;
	}

	loss = randLosses[idx];
	targetRatios = randRatios[idx];

	if (verbose) {
		MLOG("Minimal loss for %d random points = %f\n", nRand, loss);

		MLOG("Ratio");
		for (int i = 0; i < nClasses; i++)
			MLOG("\t%.4f", targetRatios[i]);
		MLOG("\n");
	}

	return loss;
}

void get_filtered_row_ids(vector<MedSample>& samples, const vector<string>& groups, vector<float>& targetRatios, vector<vector<float>>& ratios, vector<string>& groups_v, vector<int>& class_idx,
	vector<int>& filtered_row_ids) {

	int nClasses = (int)targetRatios.size();
	int nGroups = (int)ratios.size();

	vector<vector<float>> samplingRatios(nGroups, vector<float>(nClasses));
	get_sampling_ratios(targetRatios, ratios, samplingRatios);

	map<string, vector<vector<int>>> indices;
	for (string& group : groups_v)
		indices[group].resize(nClasses);

	for (int i = 0; i < samples.size(); i++)
		indices[groups[i]][class_idx[(int)samples[i].outcome]].push_back(i);

	for (size_t i = 0; i < nGroups; i++) {
		for (int j = 0; j < nClasses; j++) {
			vector<int>& vec = indices[groups_v[i]][j];
			shuffle(vec.begin(), vec.end(), globalRNG::get_engine());
			filtered_row_ids.insert(filtered_row_ids.end(), vec.begin(), vec.begin() + (int)(0.5 + samplingRatios[i][j] * vec.size()));
		}
	}
	sort(filtered_row_ids.begin(), filtered_row_ids.end());
}

float medial::process::match_multi_class(MedFeatures& data, const vector<string> &groups, vector<int> &filtered_row_ids, vector<float>& price_ratios, int nRand, int verbose) {

	vector<float> targetRatios;
	vector <vector<float>> ratios;
	vector<string> groups_v;
	vector<int> class_idx;
	vector<vector<int>> counts;
	float loss = get_matching_dist(data.samples, groups, targetRatios, price_ratios, nRand, verbose, ratios, counts, groups_v, class_idx);

	// sample according to optimal ratio
	get_filtered_row_ids(data.samples, groups, targetRatios, ratios, groups_v, class_idx, filtered_row_ids);

	filter_row_indexes(data, filtered_row_ids);
	return loss;
}

float medial::process::match_multi_class(vector<MedSample>& samples, const vector<string> &groups, vector<int> &filtered_row_ids, vector<float>& price_ratios, int nRand, int verbose) {

	vector<float> targetRatios;
	vector <vector<float>> ratios;
	vector<string> groups_v;
	vector<int> class_idx;
	vector<vector<int>> counts;
	float loss = get_matching_dist(samples, groups, targetRatios, price_ratios, nRand, verbose, ratios, counts, groups_v, class_idx);

	// sample according to optimal ratio
	get_filtered_row_ids(samples, groups, targetRatios, ratios, groups_v, class_idx, filtered_row_ids);
	if (verbose)
		MLOG("filtered_row_ids size = %d\n", (int)filtered_row_ids.size());

	for (size_t i = 0; i < filtered_row_ids.size(); i++)
		samples[i] = samples[filtered_row_ids[i]];
	samples.resize(filtered_row_ids.size());

	return loss;
}

void medial::process::match_multi_class_to_dist(MedFeatures& data, const vector<string> &groups, vector<int> &filtered_row_ids, vector<float> probs) {

	vector<vector<float>> ratios;
	vector<string> groups_v;
	vector<int> class_idx;
	vector<vector<int>> counts;

	prepare_for_matching(data.samples, groups, class_idx, groups_v, ratios, counts, 0);
	get_filtered_row_ids(data.samples, groups, probs, ratios, groups_v, class_idx, filtered_row_ids);
	filter_row_indexes(data, filtered_row_ids);

	return;
}

void medial::process::match_multi_class_to_dist(vector<MedSample>& samples, const vector<string> &groups, vector<int> &filtered_row_ids, vector<float> probs) {

	vector<vector<float>> ratios;
	vector<string> groups_v;
	vector<int> class_idx;
	vector<vector<int>> counts;

	prepare_for_matching(samples, groups, class_idx, groups_v, ratios, counts, 0);
	get_filtered_row_ids(samples, groups, probs, ratios, groups_v, class_idx, filtered_row_ids);

	for (size_t i = 0; i < filtered_row_ids.size(); i++)
		samples[i] = samples[filtered_row_ids[i]];
	samples.resize(filtered_row_ids.size());

	return;
}

void medial::process::split_matrix(const MedFeatures& matrix, vector<int>& folds, int iFold,
	MedFeatures& trainMatrix, MedFeatures& testMatrix, const vector<string> *selected_features) {
	MedFeatures *matrixPtrs[2] = { &trainMatrix, &testMatrix };

	// Prepare
	vector<string> features;
	if (selected_features == NULL || selected_features->empty())
		matrix.get_feature_names(features);
	else
		features = *selected_features;

	for (const string& ftr : features) {
		trainMatrix.attributes[ftr] = matrix.attributes.at(ftr);
		testMatrix.attributes[ftr] = matrix.attributes.at(ftr);
	}
	//realloc memory:
	int selection_cnt = 0;
	for (int i = 0; i < matrix.samples.size(); i++)
		selection_cnt += int(folds[i] == iFold);
	trainMatrix.samples.reserve((int)matrix.samples.size() - selection_cnt);
	testMatrix.samples.reserve(selection_cnt);
	for (string& ftr : features) {
		trainMatrix.data[ftr].reserve((int)matrix.samples.size() - selection_cnt);
		testMatrix.data[ftr].reserve(selection_cnt);
	}
	if (!matrix.weights.empty()) {
		trainMatrix.weights.reserve((int)matrix.samples.size() - selection_cnt);
		testMatrix.weights.reserve(selection_cnt);
	}

	// Fill
	for (int i = 0; i < matrix.samples.size(); i++) {
		int ptrIdx = (folds[i] == iFold) ? 1 : 0;
		matrixPtrs[ptrIdx]->samples.push_back(matrix.samples[i]);
		for (string& ftr : features)
			matrixPtrs[ptrIdx]->data[ftr].push_back(matrix.data.at(ftr)[i]);

		if (!matrix.weights.empty())
			matrixPtrs[ptrIdx]->weights.push_back(matrix.weights[i]);
	}
}

void medial::process::split_matrix(const MedFeatures& matrix, unordered_map<int, int>& folds, int iFold,
	MedFeatures& trainMatrix, MedFeatures& testMatrix, const vector<string> *selected_features) {
	MedFeatures *matrixPtrs[2] = { &trainMatrix, &testMatrix };

	// Prepare
	vector<string> features;
	if (selected_features == NULL || selected_features->empty())
		matrix.get_feature_names(features);
	else
		features = *selected_features;

	for (const string& ftr : features) {
		trainMatrix.attributes[ftr] = matrix.attributes.at(ftr);
		testMatrix.attributes[ftr] = matrix.attributes.at(ftr);
	}
	//realloc memory:
	int selection_cnt = 0;
	for (int i = 0; i < matrix.samples.size(); i++)
		selection_cnt += int(folds[matrix.samples[i].id] == iFold);
	trainMatrix.samples.reserve((int)matrix.samples.size() - selection_cnt);
	testMatrix.samples.reserve(selection_cnt);
	for (string& ftr : features) {
		trainMatrix.data[ftr].reserve((int)matrix.samples.size() - selection_cnt);
		testMatrix.data[ftr].reserve(selection_cnt);
	}
	if (!matrix.weights.empty()) {
		trainMatrix.weights.reserve((int)matrix.samples.size() - selection_cnt);
		testMatrix.weights.reserve(selection_cnt);
	}

	// Fill
	for (int i = 0; i < matrix.samples.size(); i++) {
		int ptrIdx = (folds[matrix.samples[i].id] == iFold) ? 1 : 0;
		matrixPtrs[ptrIdx]->samples.push_back(matrix.samples[i]);
		for (string& ftr : features)
			matrixPtrs[ptrIdx]->data[ftr].push_back(matrix.data.at(ftr)[i]);

		if (!matrix.weights.empty())
			matrixPtrs[ptrIdx]->weights.push_back(matrix.weights[i]);
	}
}

void medial::process::convert_prctile(vector<float> &features_prctiles) {
	unordered_map<float, vector<int>> val_to_inds;
	vector<float> sorted_uniqu_vals;
	for (int k = 0; k < features_prctiles.size(); ++k) {
		float val = features_prctiles[k];
		if (val_to_inds.find(val) == val_to_inds.end())
			sorted_uniqu_vals.push_back(val);
		val_to_inds[val].push_back(k);
	}
	sort(sorted_uniqu_vals.begin(), sorted_uniqu_vals.end());
	int cum_sum_size = 0;
	for (size_t k = 0; k < sorted_uniqu_vals.size(); ++k)
	{
		float prctile_val = float(cum_sum_size) / features_prctiles.size();
		//set this prctile val in all indexes:
		for (int ind : val_to_inds[sorted_uniqu_vals[k]])
			features_prctiles[ind] = prctile_val;
		cum_sum_size += (int)val_to_inds[sorted_uniqu_vals[k]].size();
	}
}

void medial::process::match_to_prior(const vector<float> &outcome,
	const vector<float> &group_values, float target_prior, vector<int> &sel_idx) {
	unordered_map<float, vector<int>> val_to_inds;
	for (size_t i = 0; i < group_values.size(); ++i)
		val_to_inds[group_values[i]].push_back((int)i);

	//sub sample each group to match this prior:
	for (auto it = val_to_inds.begin(); it != val_to_inds.end(); ++it)
	{
		double grp_prior = 0;
		vector<vector<int>> grp_inds(2);
		for (size_t i = 0; i < it->second.size(); ++i)
			grp_inds[outcome[it->second[i]] > 0].push_back(it->second[i]);
		grp_prior = double(grp_inds[1].size()) / it->second.size();
		int grp_sel = int(grp_prior > target_prior);
		vector<int> *inds = &grp_inds[grp_sel];
		int sub_sample_count;
		if (grp_prior > target_prior)
			sub_sample_count = target_prior * grp_inds[1 - grp_sel].size() / (1 - target_prior);
		else
			sub_sample_count = (1 - target_prior) * grp_inds[1 - grp_sel].size() / target_prior;
		if (sub_sample_count > inds->size())
			sub_sample_count = (int)inds->size();
		shuffle(inds->begin(), inds->end(), globalRNG::get_engine());
		inds->resize(sub_sample_count); //subsample in inds

										//add fully groups
		sel_idx.insert(sel_idx.end(), grp_inds[0].begin(), grp_inds[0].end());
		sel_idx.insert(sel_idx.end(), grp_inds[1].begin(), grp_inds[1].end());
	}
}

double medial::process::match_to_prior(MedSamples &samples, float target_prior, vector<int> &sel_idx) {
	int size_f = samples.nSamples();
	if (size_f == 0)
		MTHROW_AND_ERR("Error : sampels is empty\n");
	vector<float> fetched_labels; fetched_labels.reserve(size_f);
	vector<float> all_in_same(size_f);
	vector<MedSample *> pointers_to_smps; pointers_to_smps.reserve(size_f);
	double pr = 0;
	for (size_t i = 0; i < samples.idSamples.size(); ++i)
		for (size_t j = 0; j < samples.idSamples[i].samples.size(); ++j)
		{
			fetched_labels.push_back(samples.idSamples[i].samples[j].outcome);
			pr += int(samples.idSamples[i].samples[j].outcome > 0);
			pointers_to_smps.push_back(&samples.idSamples[i].samples[j]);
		}
	pr /= size_f;
	if (pr < target_prior) {
		vector<MedIdSamples> to_change;
		medial::process::match_to_prior(fetched_labels, all_in_same, target_prior, sel_idx);
		medial::process::commit_selection(pointers_to_smps, sel_idx);
		//sort by pid:
		sort(pointers_to_smps.begin(), pointers_to_smps.end(), [](const MedSample *a, const MedSample *b) {
			if (a->id == b->id)
				return b->time > a->time;
			return b->id > a->id;
		});
		//aggregate pointers_to_smps into to_change:
		for (size_t i = 0; i < pointers_to_smps.size(); ++i)
		{
			if (to_change.empty() || to_change.back().id != pointers_to_smps[i]->id) {
				MedIdSamples smp(pointers_to_smps[i]->id);
				smp.split = pointers_to_smps[i]->split;
				to_change.push_back(smp);
			}
			to_change.back().samples.push_back(*pointers_to_smps[i]);
		}
		MLOG("Changing prior: was %2.3f%% and changed to %2.3f%%\n", 100 * pr, 100 * target_prior);
		samples.idSamples = move(to_change);
		samples.sort_by_id_date();
		medial::print::print_samples_stats(samples);
	}
	return pr;
}

double medial::process::match_to_prior(MedFeatures &features, float target_prior, vector<int> &sel_idx) {
	int size_f = (int)features.samples.size();
	if (size_f == 0)
		MTHROW_AND_ERR("Error : sampels is empty\n");
	vector<float> fetched_labels; fetched_labels.reserve(size_f);
	vector<float> all_in_same(size_f);
	double pr = 0;
	for (size_t i = 0; i < features.samples.size(); ++i) {
		fetched_labels.push_back(features.samples[i].outcome);
		pr += int(features.samples[i].outcome > 0);
	}
	pr /= size_f;
	if (pr < target_prior) {
		vector<MedIdSamples> to_change;
		medial::process::match_to_prior(fetched_labels, all_in_same, target_prior, sel_idx);
		medial::process::filter_row_indexes(features, sel_idx);

		MLOG("Changing prior: was %2.3f%% and changed to %2.3f%%\n", 100 * pr, 100 * target_prior);
		medial::print::print_samples_stats(features.samples);
	}
	return pr;
}

void medial::process::match_to_prior(MedFeatures &features,
	const vector<string> &group_values, float target_prior, vector<int> &sel_idx, bool print_verbose) {
	if (target_prior <= 0 || target_prior >= 1)
		MTHROW_AND_ERR("Error - medial::process::match_to_prior - bad target_prior (%f), shoulf by between 0 -1\n",
			target_prior);
	//rewrote it again to accept vector<string> and not vector<float> as group_values

	sel_idx.clear();
	unordered_map<string, vector<int>> val_to_inds;
	for (size_t i = 0; i < group_values.size(); ++i)
		val_to_inds[group_values[i]].push_back((int)i);

	if (print_verbose) {
		MLOG("Before Matching Size=%d:\n", (int)features.samples.size());
		print_stats(features, sel_idx, group_values);
	}

	//sub sample each group to match this prior:
	for (auto it = val_to_inds.begin(); it != val_to_inds.end(); ++it)
	{
		double grp_prior = 0;
		vector<vector<int>> grp_inds(2);
		for (size_t i = 0; i < it->second.size(); ++i)
			grp_inds[features.samples[it->second[i]].outcome > 0].push_back(it->second[i]);
		grp_prior = double(grp_inds[1].size()) / it->second.size();
		int grp_sel = int(grp_prior > target_prior);
		vector<int> *inds = &grp_inds[grp_sel];
		int sub_sample_count;
		if (grp_prior > target_prior)
			sub_sample_count = target_prior * grp_inds[1 - grp_sel].size() / (1 - target_prior);
		else
			sub_sample_count = (1 - target_prior) * grp_inds[1 - grp_sel].size() / target_prior;
		if (sub_sample_count > inds->size())
			sub_sample_count = (int)inds->size();
		shuffle(inds->begin(), inds->end(), globalRNG::get_engine());
		inds->resize(sub_sample_count); //subsample in inds

										//add fully groups
		sel_idx.insert(sel_idx.end(), grp_inds[0].begin(), grp_inds[0].end());
		sel_idx.insert(sel_idx.end(), grp_inds[1].begin(), grp_inds[1].end());
	}
	filter_row_indexes(features, sel_idx);

	if (print_verbose) {
		MLOG("After Matching Size=%d:\n", (int)features.samples.size());
		print_stats(features, sel_idx, group_values);
	}
}

/// Return number of splits, also check mismatches between idSample and internal MedSamples and set idSamples.split if missing
int medial::process::nSplits(vector<MedSample>& samples) {

	int maxSplit = -2;
	unordered_set<int> splits;
	for (auto& sample : samples) {
		if (sample.split > maxSplit)
			maxSplit = sample.split;
		splits.insert(sample.split);
	}

	if (maxSplit == -1)
		return 1;
	else {
		if (splits.find(-1) != splits.end())
			MTHROW_AND_ERR("Cannot handle -1 and splits\n");

		for (int split = 0; split < maxSplit; split++) {
			if (splits.find(split) == splits.end())
				MTHROW_AND_ERR("Missing split %d\n", split);
		}
		return maxSplit + 1;
	}

}


//-------------------------------------------------------------------------------------------------------
// assumed data is already initialized
int MedFeatures::init_masks()
{
	masks.clear();
	if (data.empty()) return 0;

	for (auto &e : data) {
		masks[e.first] = vector<unsigned char>();
		masks[e.first].resize(data[e.first].size(), 0);
	}
	return 0;
}

//-------------------------------------------------------------------------------------------------------
int MedFeatures::get_masks_as_mat(MedMat<unsigned char> &masks_mat)
{
	masks_mat.clear();
	if (masks.empty()) return 0;

	int nrows = (int)masks.begin()->second.size();
	int ncols = (int)masks.size();
	masks_mat.resize(nrows, ncols);

	vector<string> names;
	get_feature_names(names);
	masks_mat.set_signals(names);
	//Set metadata from samples:
	masks_mat.recordsMetadata.resize(samples.size());
	for (size_t i = 0; i < samples.size(); ++i)
	{
		RecordData &d = masks_mat.recordsMetadata[i];
		const MedSample &s = samples[i];
		d.id = s.id;
		d.date = s.time;
		d.outcomeTime = s.outcomeTime;
		d.split = s.split;
		d.label = s.outcome;
		if (!s.prediction.empty())
			d.pred = s.prediction[0];
		if (!weights.empty())
			d.weight = weights[i];
	}

#pragma omp parallel for
	for (int i = 0; i < names.size(); i++) {
		unsigned char *p_mask = masks[names[i]].data();
		for (int j = 0; j < nrows; j++)
			masks_mat(j, i) = p_mask[j];
	}

	return 0;
}

//-------------------------------------------------------------------------------------------------------
int MedFeatures::mark_imputed_in_masks(float _missing_val)
{
	if (masks.empty()) init_masks();
	if (masks.empty()) return 0;

	int nrows = (int)masks.begin()->second.size();
	//int ncols = (int)masks.size();
	vector<string> names;
	get_feature_names(names);


#pragma omp parallel for
	for (int i = 0; i < names.size(); i++) {
		unsigned char *p_mask = masks[names[i]].data();
		float *p_data = data[names[i]].data();
		for (int j = 0; j < nrows; j++)
			if (p_data[j] == _missing_val)
				p_mask[j] |= MedFeatures::imputed_mask;
	}

	return 0;
}


//-------------------------------------------------------------------------------------------------------
void MedFeatures::round_data(float r)
{
	if (r <= 0) return;
	vector<string> names;
	get_feature_names(names);
	int n_feat = (int)names.size();

#pragma omp parallel for
	for (int j = 0; j < n_feat; j++) {
		float *p_data = data[names[j]].data();
		int len = (int)data[names[j]].size();
		for (int i = 0; i < len; i++)
			p_data[i] = roundf(p_data[i] / r)*r;
	}

}


//-------------------------------------------------------------------------------------------------------
void MedFeatures::noise_data(float r)
{
	if (r <= 0) return;
	vector<string> names;
	get_feature_names(names);
	int n_feat = (int)names.size();

#pragma omp parallel for
	for (int j = 0; j < n_feat; j++) {
		float *p_data = data[names[j]].data();
		int len = (int)data[names[j]].size();
		for (int i = 0; i < len; i++) {
			float noise = r * (rand_1()*2.0f - 1.0f);
			p_data[i] = p_data[i] + noise;
		}
	}

}

// Sort Features by id + time
//-------------------------------------------------------------------------------------------------------
void MedFeatures::samples_sort() {

	int nSamples = (int)samples.size();
	vector<pair<int, MedSample>> sorted_inds(nSamples);
	for (int i = 0; i < nSamples; i++) {
		sorted_inds[i].first = i;
		sorted_inds[i].second = samples[i];
	}

	sort(sorted_inds.begin(), sorted_inds.end(),
		[](const pair<int, MedSample> &c1, const pair<int, MedSample> &c2) {return ((c1.second.id < c2.second.id) || (c1.second.id == c2.second.id && c1.second.time < c2.second.time)); });

	vector<MedSample> origSamples = samples;
	for (int i = 0; i < nSamples; i++)
		samples[i] = origSamples[sorted_inds[i].first];

	if (weights.size()) {
		vector<float> origWeights = weights;
		for (int i = 0; i < nSamples; i++)
			weights[i] = origWeights[sorted_inds[i].first];
	}

	for (auto& rec : data) {
		string name = rec.first;
		vector<float> origData = data[name];
		for (int i = 0; i < nSamples; i++)
			data[name][i] = origData[sorted_inds[i].first];

		if (masks.find(name) != masks.end()) {
			vector<unsigned char> origMask = masks[name];
			for (int i = 0; i < nSamples; i++)
				masks[name][i] = origMask[sorted_inds[i].first];
		}
	}

}

// Get feature name that matches a substring
//-------------------------------------------------------------------------------------------------------
string MedFeatures::resolve_name(string& substr) const {
	// Exact name ?
	if (data.find(substr) != data.end())
		return substr;
	else {
		vector<string> names;
		get_feature_names(names);
		return names[find_in_feature_names(names, substr)];
	}
}
