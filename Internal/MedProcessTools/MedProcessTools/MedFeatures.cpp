#define _CRT_SECURE_NO_WARNINGS

#include "MedFeatures.h"
#include <random>

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
	for (int i = 0; i < (int)datap.size(); i++) {
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
	int index = 0;
	for (auto& ss : samples) {
		RecordData rd;
		rd.time = (long)ss.outcomeTime; rd.label = ss.outcome; ; rd.split = ss.split; rd.weight = 0.0;
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
		mat.normalized_flag &= (int)attributes[name].normalized;
	for (string& name : namesToTake)
		mat.normalized_flag &= (int)attributes[name].imputed;

	int index = 0;
	mat.signals.insert(mat.signals.end(), namesToTake.begin(), namesToTake.end());
	for (int i = 0; i < nrows; i++) {
		//for (auto& ss : samples) {
		auto &ss = samples[idx[i]];
		RecordData rd;
		rd.time = (long)ss.outcomeTime;  rd.weight = 0.0;  rd.label = ss.outcome; ; rd.split = ss.split;
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
		smp.outcomeTime = rd.time;
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
// Return -1 upon failure to open file, 0 upon success
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
		for (int j = 0; j < n_preds; j++)
			out_f << ",pred_" + to_string(j);
	}

	// names of features
	for (int j = 0; j < col_names.size(); j++)
		out_f << "," + col_names[j];
	out_f << "\n";

	// data
	for (int i = 0; i < samples.size(); i++) {

		out_f << to_string(i); // serial
		if (weights.size()) out_f << "," + to_string(weights[i]); // Weights

		// sample
		out_f << "," + to_string(samples[i].id);
		out_f << "," + to_string(samples[i].time);
		out_f << "," + to_string(samples[i].outcome);
		out_f << "," + to_string(samples[i].outcomeTime);
		out_f << "," + to_string(samples[i].split);

		// predictions
		for (int j = 0; j < n_preds; j++)
			out_f << "," + to_string(samples[i].prediction[j]);

		// features
		for (int j = 0; j < col_names.size(); j++)
			out_f << "," + to_string(data[col_names[j]][i]);
		out_f << "\n";
	}

	out_f.close();
	return 0;
}

// Read features (samples + weights + data) from a csv file with a header line
// Return -1 upon failure to open file, 0 upon success
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
		return -1;
	}

	int ncols = -1;
	string curr_line;
	vector<string> names;
	int weighted = 0;
	int max_pred = -1;
	vector<int> field_ind_to_pred_index;
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
				if (!boost::starts_with(fields[i], "pred_")) {
					data[fields[i]] = vector<float>();
					attributes[fields[i]].normalized = attributes[fields[i]].imputed = false;
					names.push_back(fields[i]);
				}
				else {
					int candi = stoi(boost::replace_all_copy(fields[i], "pred_", ""));
					if (candi > max_pred) {
						max_pred = candi;
						field_ind_to_pred_index.resize(max_pred + 1, -1);
						field_ind_to_pred_index[max_pred] = i;
					}
				}
			}
			ncols = (int)fields.size();
		}
		else { // Data lines
			if (fields.size() != ncols) {
				string msg = "expected " + to_string(ncols) + " fields, got " + to_string((int)fields.size()) + "fields in line: " + curr_line.c_str() + "\n";
				throw runtime_error(msg.c_str());
			}

			idx++;

			if (weighted)
				weights.push_back(stof(fields[idx++]));

			MedSample newSample;
			newSample.id = stoi(fields[idx++]);
			newSample.time = stoi(fields[idx++]);
			newSample.outcome = stof(fields[idx++]);
			newSample.outcomeTime = stoi(fields[idx++]);
			newSample.split = stoi(fields[idx++]);

			newSample.prediction.resize(max_pred + 1);

			for (int i = 0; i <= max_pred; ++i)
				if (field_ind_to_pred_index[i] >= 0)
					newSample.prediction[i] = stof(fields[field_ind_to_pred_index[i]]);
			samples.push_back(newSample);

			for (int i = 0; i < names.size(); i++)
				data[names[i]].push_back(stof(fields[idx + i]));
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

// Get the corresponding MedSamples object .  Assuming samples are ordered in features (all id's samples are consecutive)
//.......................................................................................
void MedFeatures::get_samples(MedSamples& outSamples) {

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
int MedFeatures::get_max_serial_id_cnt() {

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

void medial::process::filter_row_indexes(MedFeatures &dataMat, vector<int> &selected_indexes, bool op_flag) {
	MedFeatures filtered;
	filtered.time_unit = dataMat.time_unit;
	filtered.attributes = dataMat.attributes;

	sort(selected_indexes.begin(), selected_indexes.end());

	int curr_ind = 0;
	if (!op_flag) {
		for (auto iit = dataMat.data.begin(); iit != dataMat.data.end(); ++iit)
			filtered.data[iit->first].reserve(selected_indexes.size());
		if (!dataMat.weights.empty())
			filtered.weights.reserve(selected_indexes.size());
		filtered.samples.reserve(selected_indexes.size());
		for (int i : selected_indexes) //all selected indexes
		{
			filtered.samples.push_back(dataMat.samples[i]);
			for (auto iit = dataMat.data.begin(); iit != dataMat.data.end(); ++iit)
				filtered.data[iit->first].push_back(iit->second[i]);

			if (!dataMat.weights.empty())
				filtered.weights.push_back(dataMat.weights[i]);
		}
	}
	else {
		for (auto iit = dataMat.data.begin(); iit != dataMat.data.end(); ++iit)
			filtered.data[iit->first].reserve((int)dataMat.samples.size() - (int)selected_indexes.size());
		if (!dataMat.weights.empty())
			filtered.weights.reserve((int)dataMat.samples.size() - (int)selected_indexes.size());
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
}

void medial::process::down_sample(MedFeatures &dataMat, double take_ratio) {
	int final_cnt = int(take_ratio * dataMat.samples.size());
	if (take_ratio >= 1) {
		return;
	}

	vector<int> all_selected_indexes(final_cnt);
	vector<bool> seen_index((int)dataMat.samples.size());
	random_device rd;
	mt19937 gen;
	uniform_int_distribution<> dist_gen(0, (int)dataMat.samples.size() - 1);
	for (size_t k = 0; k < final_cnt; ++k) //for 0 and 1:
	{
		int num_ind = dist_gen(gen);
		while (seen_index[num_ind])
			num_ind = dist_gen(gen);
		seen_index[num_ind] = true;

		all_selected_indexes[k] = num_ind;
	}
	filter_row_indexes(dataMat, all_selected_indexes);
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
	int init_sample_count = (int)data_records.samples.size();

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
		MLOG("Before Mathcing Total samples is %d on %d groups\n",
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
		if (print_verbose) {
			cout << grp << "\t" << count_label_groups[0][grp] << "\t" << count_label_groups[1][grp] << "\t"
				<< count_label_groups[1][grp] / float(count_label_groups[1][grp] + count_label_groups[0][grp]) << endl;
		}

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

		if (print_verbose)
			if (factor > 0)
				MLOG("Weighting group %s base_ratio=%f factor= %f\n",
					grp.c_str(), base_ratio, factor);
			else
				MLOG("Dropping group %s Num_controls=%d with zero cases\n",
					grp.c_str(), count_label_groups[0][grp]);

		if (factor == 0) {
			list_label_groups[0].erase(grp); //for zero it's different
			list_label_groups[1].erase(grp); //for zero it's different
		}
		else {
			group_to_factor[grp] = factor;
			//Commit Change to weights:
			for (int ind : list_label_groups[1][grp])
				full_weights[ind] = factor;
		}
	}
	//let's divide all by 1/max_factor:
	double total_cnt_after = 0;
	for (auto it = group_to_factor.begin(); it != group_to_factor.end(); ++it)
		total_cnt_after += (count_label_groups[0][it->first] + count_label_groups[1][it->first]) * it->second;
	max_factor = float(total_cnt_after / (double)init_sample_count); //change to take not max - keep same sum

	if (max_factor > 0) {
		for (size_t i = 0; i < full_weights.size(); ++i)
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
		MLOG("Group\tCount_0\tCount_1\tratio\tweighted_ratio\n");
		for (const string &grp : all_groups)
			if (group_to_factor.find(grp) != group_to_factor.end())
				MLOG("%s\t%d\t%d\t%2.5f\t%2.5f\n", grp.c_str(),
					counts_stat[grp][0], counts_stat[grp][1],
					counts_stat[grp][1] / double(counts_stat[grp][1] + counts_stat[grp][0]),
					weight_stats[grp][1] / double(weight_stats[grp][1] + weight_stats[grp][0]));
		//print_by_year(data_records.samples);
	}
	if (max_factor > 0)
		return 1.0 / max_factor;
	else
		return (double)0;
}

void  medial::process::match_by_general(MedFeatures &data_records, const vector<string> &groups,
	vector<int> &filtered_row_ids, float price_ratio, bool print_verbose) {
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

	unordered_map<string, int> year_total;
	unordered_map<string, float> year_ratio;
	vector<float> all_ratios((int)all_groups.size());
	int i = 0;
	sort(all_groups.begin(), all_groups.end());
	if (print_verbose) {
		MLOG("Before Mathcing Total samples is %d on %d groups\n",
			(int)data_records.samples.size(), (int)all_groups.size());
		MLOG("Group"  "\t"  "Count_0"  "\t"  "Count_1"  "\t"  "ratio\n");
	}
	for (const string &grp : all_groups)
	{
		if (count_label_groups[0][grp] == 0)
			continue;
		year_total[grp] = count_label_groups[0][grp] + count_label_groups[1][grp];
		year_ratio[grp] = count_label_groups[1][grp] / float(count_label_groups[0][grp] + count_label_groups[1][grp]);
		all_ratios[i] = year_ratio[grp];
		++i;
		if (print_verbose)
			MLOG("%s\t%d\t%d\t%f\n", grp.c_str(), count_label_groups[0][grp], count_label_groups[1][grp]
				, count_label_groups[1][grp] / float(count_label_groups[1][grp] + count_label_groups[0][grp]));

	}


	//Choose ratio to balance all for:
	sort(all_ratios.begin(), all_ratios.end());
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

	//For each year_bin - balance to this ratio using price_ratio weight for removing 1's labels:
	seen_pid_0.clear();
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

		unordered_set<int> seen_year_pid;
		random_shuffle(list_label_groups[ind][grp].begin(), list_label_groups[ind][grp].end());
		if (target_size > list_label_groups[ind][grp].size())
			MERR("ERROR/BUG: try to shrink %d into %d\n"
				, (int)list_label_groups[ind][grp].size(), target_size);
		else
			list_label_groups[ind][grp].resize(target_size); //for zero it's different
	}

	//Commit on all records:
	MedFeatures filtered;
	filtered.time_unit = data_records.time_unit;
	filtered.attributes = data_records.attributes;

	for (size_t k = 0; k < list_label_groups.size(); ++k) //for 0 and 1:
		for (auto it = list_label_groups[k].begin(); it != list_label_groups[k].end(); ++it) //for each year
		{
			vector<int> ind_list = it->second;
			filtered_row_ids.insert(filtered_row_ids.end(), ind_list.begin(), ind_list.end());
		}


	filter_row_indexes(data_records, filtered_row_ids);

	if (print_verbose) {
		MLOG("After Matching Size=%d:\n", (int)data_records.samples.size());
		unordered_map<string, vector<int>> counts_stat;
		for (size_t i = 0; i < filtered_row_ids.size(); ++i)
		{
			if (counts_stat[groups[filtered_row_ids[i]]].empty())
				counts_stat[groups[filtered_row_ids[i]]].resize(2);
			++counts_stat[groups[filtered_row_ids[i]]][data_records.samples[i].outcome > 0];
		}
		MLOG("Group\tCount_0\tCount_1\tratio\n");
		for (const string &grp : all_groups)
			MLOG("%s\t%d\t%d\t%2.5f\n", grp.c_str(),
				counts_stat[grp][0], counts_stat[grp][1], counts_stat[grp][1] / double(counts_stat[grp][1] + counts_stat[grp][0]));
		//print_by_year(data_records.samples);
	}
}
