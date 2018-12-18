#define _CRT_SECURE_NO_WARNINGS

#include <MedProcessTools/MedProcessTools/MedFeatures.h>
#include <MedUtils/MedUtils/MedUtils.h>
#include <random>
#include <omp.h>
#include <MedIO/MedIO/MedIO.h>

#include <Logger/Logger/Logger.h>
#define LOCAL_SECTION LOG_MEDFEAT
#define LOCAL_LEVEL	LOG_DEF_LEVEL

int MedFeatures::global_serial_id_cnt = 0;

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
		mat.normalized_flag &= (int)attributes.at(name).normalized;
	for (string& name : namesToTake)
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
int MedFeatures::write_as_csv_mat(const string &csv_fname) const
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
			out_f << ",pred_" << to_string(j);
	}

	// names of features
	for (int j = 0; j < col_names.size(); j++)
		out_f << "," << col_names[j];
	out_f << "\n";

	// data
	for (int i = 0; i < samples.size(); i++) {

		out_f << to_string(i); // serial
		if (weights.size()) out_f << "," + to_string(weights[i]); // Weights

																  // sample
		out_f << "," << samples[i].id;
		out_f << "," << samples[i].time;
		out_f << "," << samples[i].outcome;
		out_f << "," << samples[i].outcomeTime;
		out_f << "," << samples[i].split;

		// predictions
		for (int j = 0; j < n_preds; j++)
			out_f << "," << samples[i].prediction[j];

		// features
		for (int j = 0; j < col_names.size(); j++)
			out_f << "," << data.at(col_names[j])[i];
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

			idx += (max_pred + 1);
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
		MLOG("Before Mathcing Total samples is %d on %d groups\n",
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
				factor_needed_up = diff_controls / diff_cases;
			}
			if (ratio_ind > 0) {
				float diff_cases = abs(cases_sum[ratio_ind] - cases_sum[ratio_ind - 1]);
				float diff_controls = abs(controls_sum[ratio_ind] - controls_sum[ratio_ind - 1]);
				factor_needed_down = diff_controls / diff_cases;
			}
			if (count_label_groups[0][grp] == 0)
				grp_ratio = 1; //just for correct printing
			if (factor_needed_up > 0)
				MLOG("%s\t%d\t%d\t%f\t[%2.2f-%2.2f]\n", grp.c_str(), count_label_groups[0][grp], count_label_groups[1][grp]
					, grp_ratio, factor_needed_down, factor_needed_up);
			else
				MLOG("%s\t%d\t%d\t%f\t[%2.2f-]\n", grp.c_str(), count_label_groups[0][grp], count_label_groups[1][grp]
					, grp_ratio, factor_needed_down);
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
	int min_grp_size = 5;
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
		random_shuffle(list_label_groups[ind][grp].begin(), list_label_groups[ind][grp].end());
		if (target_size > list_label_groups[ind][grp].size())
			MERR("ERROR/BUG: try to shrink %d into %d\n"
				, (int)list_label_groups[ind][grp].size(), target_size);
		else
			list_label_groups[ind][grp].resize(target_size); //for zero it's different
	}

	for (int k = 0; k < skip_grp_indexs.size(); ++k)
		all_groups.erase(all_groups.begin() + skip_grp_indexs[k]);

	//Commit on all records:
	MedFeatures filtered;
	filtered.time_unit = data_records.time_unit;
	filtered.attributes = data_records.attributes;

	for (size_t k = 0; k < list_label_groups.size(); ++k) //for 0 and 1:
		for (auto it = list_label_groups[k].begin(); it != list_label_groups[k].end(); ++it) //for each year
		{
			vector<int> &ind_list = it->second;
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
		random_shuffle(inds->begin(), inds->end());
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
		sort(pointers_to_smps.begin(), pointers_to_smps.end(), [](MedSample *a, MedSample *b) {
			int pos = 0;
			if (a->id == b->id)
				return b->time > a->time;
			return b->id > a->id;
		});
		//aggregate pointers_to_smps into to_change:
		for (size_t i = 0; i < pointers_to_smps.size(); ++i)
		{
			if (to_change.empty() || to_change.back().id != pointers_to_smps[i]->id) {
				MedIdSamples smp(pointers_to_smps[i]->id);
				to_change.push_back(smp);
			}
			to_change.back().samples.push_back(*pointers_to_smps[i]);
		}
		MLOG("Changing prior: was %2.3f%% and changed to %2.3f%%\n", 100 * pr, 100 * target_prior);
		samples.idSamples.swap(to_change);
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


template<class T> float medial::stats::get_preds_auc_q(const vector<T> &preds, const vector<float> &y,
	const vector<float> *weights) {
	vector<T> pred_threshold;
	unordered_map<T, vector<int>> pred_indexes;
	double tot_true_labels = 0, tot_false_labels = 0;
	bool has_weights = weights != NULL && !weights->empty();
	if (has_weights)
		for (size_t i = 0; i < preds.size(); ++i)
		{
			pred_indexes[preds[i]].push_back((int)i);
			tot_true_labels += int(y[i] > 0) * (*weights)[i];
			tot_false_labels += int(y[i] <= 0) * (*weights)[i];
		}
	else
		for (size_t i = 0; i < preds.size(); ++i)
		{
			pred_indexes[preds[i]].push_back((int)i);
			tot_true_labels += int(y[i] > 0);
			tot_false_labels += int(y[i] <= 0);
		}
	pred_threshold.resize((int)pred_indexes.size());
	auto it = pred_indexes.begin();
	for (size_t i = 0; i < pred_threshold.size(); ++i)
	{
		pred_threshold[i] = it->first;
		++it;
	}
	sort(pred_threshold.begin(), pred_threshold.end());


	//From up to down sort:
	double t_cnt = 0;
	double f_cnt = 0;
	vector<float> true_rate = vector<float>((int)pred_indexes.size());
	vector<float> false_rate = vector<float>((int)pred_indexes.size());
	int st_size = (int)pred_threshold.size() - 1;
	if (has_weights)
		for (int i = st_size; i >= 0; --i)
		{
			vector<int> &indexes = pred_indexes[pred_threshold[i]];
			//calc AUC status for this step:
			for (int ind : indexes)
			{
				bool true_label = y[ind] > 0;
				t_cnt += int(true_label) * (*weights)[ind];
				f_cnt += int(!true_label) * (*weights)[ind];
			}
			true_rate[st_size - i] = float(t_cnt / tot_true_labels);
			false_rate[st_size - i] = float(f_cnt / tot_false_labels);
		}
	else
		for (int i = st_size; i >= 0; --i)
		{
			vector<int> &indexes = pred_indexes[pred_threshold[i]];
			//calc AUC status for this step:
			for (int ind : indexes)
			{
				bool true_label = y[ind] > 0;
				t_cnt += int(true_label);
				f_cnt += int(!true_label);
			}
			true_rate[st_size - i] = float(t_cnt / tot_true_labels);
			false_rate[st_size - i] = float(f_cnt / tot_false_labels);
		}

	float auc = false_rate[0] * true_rate[0] / 2;
	for (size_t i = 1; i < true_rate.size(); ++i)
		auc += (false_rate[i] - false_rate[i - 1]) * (true_rate[i - 1] + true_rate[i]) / 2;
	return auc;
}
template float medial::stats::get_preds_auc_q<float>(const vector<float> &preds, const vector<float> &y, const vector<float> *weights);
template float medial::stats::get_preds_auc_q<double>(const vector<double> &preds, const vector<float> &y, const vector<float> *weights);

template<class T> double medial::stats::mean_vec(const vector<T> &v, const vector<float> *weights) {
	bool has_weights = weights != NULL && !weights->empty();

	double s = 0, c = 0;
	if (has_weights)
		for (size_t i = 0; i < v.size(); ++i) {
			s += v[i] != MED_MAT_MISSING_VALUE ? v[i] * (*weights)[i] : 0;
			c += (*weights)[i] * int(v[i] != MED_MAT_MISSING_VALUE);
		}
	else
		for (size_t i = 0; i < v.size(); ++i) {
			s += v[i] != MED_MAT_MISSING_VALUE ? v[i] : 0;
			c += int(v[i] != MED_MAT_MISSING_VALUE);
		}

	if (c == 0)
		return MED_MAT_MISSING_VALUE;
	return s / c;
}
template double medial::stats::mean_vec<float>(const vector<float> &v, const vector<float> *weights);
template double medial::stats::mean_vec<double>(const vector<double> &v, const vector<float> *weights);

template<class T> double medial::stats::std_vec(const vector<T> &v, T mean, const vector<float> *weights) {
	bool has_weights = weights != NULL && !weights->empty();

	double s = 0, c = 0;
	if (has_weights)
		for (size_t i = 0; i < v.size(); ++i) {
			s += v[i] != MED_MAT_MISSING_VALUE ? (*weights)[i] * (v[i] - mean) * (v[i] - mean) : 0;
			c += (*weights)[i] * int(v[i] != MED_MAT_MISSING_VALUE);
		}
	else
		for (size_t i = 0; i < v.size(); ++i) {
			s += v[i] != MED_MAT_MISSING_VALUE ? (v[i] - mean) * (v[i] - mean) : 0;
			c += int(v[i] != MED_MAT_MISSING_VALUE);
		}

	if (c == 0)
		return MED_MAT_MISSING_VALUE;
	return sqrt(s / c);
}
template double medial::stats::std_vec<float>(const vector<float> &v, float mean, const vector<float> *weights);
template double medial::stats::std_vec<double>(const vector<double> &v, double mean, const vector<float> *weights);

template <typename T, typename S>
double medial::stats::get_kendall_tau(const vector<T>& preds, const vector<S>& y, const vector<float> *weights) {
	//return kendallTau(preds, y);
	double tau = 0, cnt = 0;
	if (weights == NULL || weights->empty()) {
		unordered_map<S, vector<T>> label_to_scores;
		for (size_t i = 0; i < y.size(); ++i)
			label_to_scores[y[i]].push_back(preds[i]);
		for (auto it = label_to_scores.begin(); it != label_to_scores.end(); ++it)
			sort(it->second.begin(), it->second.end());
		for (auto it = label_to_scores.begin(); it != label_to_scores.end(); ++it)
		{
			auto bg = it;
			++bg;
			vector<T> *preds = &it->second;
			int pred_i_bigger;
			double pred_i_smaller;
			for (auto jt = bg; jt != label_to_scores.end(); ++jt)
			{
				vector<T> *preds_comp = &jt->second;
				double p_size = (double)preds_comp->size();
				for (T pred : *preds)
				{
					pred_i_bigger = medial::process::binary_search_position(preds_comp->data(), preds_comp->data() + preds_comp->size() - 1, pred);
					pred_i_smaller = p_size - medial::process::binary_search_position_last(preds_comp->data(), preds_comp->data() + preds_comp->size() - 1, pred);
					if (it->first > jt->first)
						//tau += pred_i_bigger;
						tau += pred_i_bigger - pred_i_smaller;
					else
						//tau += pred_i_smaller;
						tau += pred_i_smaller - pred_i_bigger;
				}
				cnt += p_size * preds->size();
			}
		}

		if (cnt > 1)
			tau /= cnt;

		return (float)tau;
	}
	unordered_map<S, vector<T>> label_to_scores;
	unordered_map<S, vector<pair<int, T>>> label_to_scores_w;
	for (size_t i = 0; i < y.size(); ++i)
		label_to_scores[y[i]].push_back(preds[i]);
	for (size_t i = 0; i < y.size(); ++i)
		label_to_scores_w[y[i]].push_back(pair<int, T>((int)i, preds[i]));
	for (auto it = label_to_scores_w.begin(); it != label_to_scores_w.end(); ++it)
		sort(it->second.begin(), it->second.end(), ComparePairBySecond<int, T>());
	for (auto it = label_to_scores.begin(); it != label_to_scores.end(); ++it)
		sort(it->second.begin(), it->second.end());

	vector<double> group_weights(label_to_scores.size());
	vector<vector<double>> group_weights_cumsum(label_to_scores.size());
	int iter = 0;
	for (auto it = label_to_scores_w.begin(); it != label_to_scores_w.end(); ++it) {
		for (size_t i = 0; i < it->second.size(); ++i) {
			group_weights[iter] += (*weights)[it->second[i].first];
			group_weights_cumsum[iter].push_back((*weights)[it->second[i].first]);
		}
		++iter;
	}
	//make cumsum:
	for (size_t i = 0; i < group_weights_cumsum.size(); ++i)
		for (size_t j = 1; j < group_weights_cumsum[i].size(); ++j)
			group_weights_cumsum[i][j] += group_weights_cumsum[i][j - 1];

	iter = 0;
	for (auto it = label_to_scores.begin(); it != label_to_scores.end(); ++it)
	{
		auto bg = it;
		++bg;
		vector<T> *preds = &it->second;
		//double i_size = preds->size();
		double i_size = group_weights[iter];
		double pred_i_bigger;
		double pred_i_smaller;
		int pred_i_bigger_i;
		int pred_i_smaller_i;
		int inside_group_idx = iter + 1;
		for (auto jt = bg; jt != label_to_scores.end(); ++jt)
		{
			vector<T> *preds_comp = &jt->second;
			//double p_size = (double)preds_comp->size();
			double p_size = group_weights[inside_group_idx];
			for (T pred : *preds)
			{
				pred_i_bigger_i = medial::process::binary_search_position(preds_comp->data(), preds_comp->data() + preds_comp->size() - 1, pred);
				pred_i_smaller_i = medial::process::binary_search_position_last(preds_comp->data(), preds_comp->data() + preds_comp->size() - 1, pred);
				if (pred_i_bigger_i < group_weights_cumsum[inside_group_idx].size())
					pred_i_bigger = group_weights_cumsum[inside_group_idx][pred_i_bigger_i];
				else
					pred_i_bigger = p_size;
				if (pred_i_smaller_i < group_weights_cumsum[inside_group_idx].size())
					pred_i_smaller = group_weights_cumsum[inside_group_idx][pred_i_smaller_i];
				else
					pred_i_smaller = p_size;
				if (it->first > jt->first)
					//tau += pred_i_bigger;
					tau += pred_i_bigger - (p_size - pred_i_smaller);
				else
					//tau += pred_i_smaller;
					tau += (p_size - pred_i_smaller) - pred_i_bigger;
			}
			cnt += p_size * i_size;
			++inside_group_idx;
		}
		++iter;
	}

	if (cnt > 0)
		tau /= cnt;

	return tau;
}
template double medial::stats::get_kendall_tau<double, double>(const vector<double>& preds, const vector<double>& y, const vector<float> *weights);
template double medial::stats::get_kendall_tau<float, float>(const vector<float>& preds, const vector<float>& y, const vector<float> *weights);
template double medial::stats::get_kendall_tau<double, float>(const vector<double>& preds, const vector<float>& y, const vector<float> *weights);
template double medial::stats::get_kendall_tau<float, double>(const vector<float>& preds, const vector<double>& y, const vector<float> *weights);

float medial::stats::get_rmse(const vector<float> &preds, const vector<float> &y, const vector<float> *weights) {
	double res = 0;
	if (weights == NULL || weights->empty()) {
		for (size_t i = 0; i < y.size(); ++i)
			res += (y[i] - preds[i]) * (y[i] - preds[i]);
		res /= y.size();
		res = sqrt(res);
		return (float)res;
	}

	double cnt = 0;
	for (size_t i = 0; i < y.size(); ++i) {
		res += (*weights)[i] * (y[i] - preds[i]) * (y[i] - preds[i]);
		cnt += (*weights)[i];
	}
	res /= cnt;
	res = sqrt(res);
	return (float)res;
}

float medial::stats::get_accuracy(const vector<float> &preds, const vector<float> &y, const vector<float> *weights) {
	double res = 0;
	if (weights == NULL || weights->empty()) {
		for (size_t i = 0; i < y.size(); ++i)
			res += y[i] == preds[i];
		return float(res / y.size());
	}
	double cnt = 0;
	for (size_t i = 0; i < y.size(); ++i) {
		res += (*weights)[i] * (y[i] == preds[i]);
		cnt += (*weights)[i];
	}
	return float(res / cnt);
}
