#include "DuplicateProcessor.h"

#define LOCAL_SECTION LOG_FEATCLEANER
#define LOCAL_LEVEL	LOG_DEF_LEVEL

void DuplicateProcessor::init_defaults() {
	processor_type = FTR_PROCESS_DUPLICATE;
	resample_cnt = 50;
}

int DuplicateProcessor::init(map<string, string>& mapper) {
	for (const auto &it : mapper) {
		//! [DuplicateProcessor::init]
		if (it.first == "resample_cnt")
			resample_cnt = med_stoi(it.second);
		else
			MTHROW_AND_ERR("Error DuplicateProcessor::init - unknown param %s\n", it.first.c_str());
		//! [DuplicateProcessor::init]
	}
	return 0;
}

int DuplicateProcessor::_apply(MedFeatures& features, unordered_set<int>& ids) {
	MedFeatures batch;
	batch.attributes = features.attributes;
	batch.tags = features.tags;
	batch.time_unit = features.time_unit;
	batch.medf_missing_value = features.medf_missing_value;
	batch.samples.resize(features.samples.size() * resample_cnt);
	for (int i = 0; i < batch.samples.size(); ++i) {
		//Id index  int(i / resample_cnt), inside sample row:=  i % resample_cnt
		int r_id = int(i / resample_cnt);
		batch.samples[i].id = features.samples[r_id].id;
		batch.samples[i].time = features.samples[r_id].time;
	}
	batch.init_pid_pos_len();

	for (auto &it : features.data) {
		vector<float> &v = batch.data[it.first];
		v.resize(batch.samples.size());
		for (int i = 0; i < batch.samples.size(); ++i) {
			int r_id = int(i / resample_cnt);
			v[i] = it.second[r_id];
		}
	}

	features = move(batch);

	return 0;
}