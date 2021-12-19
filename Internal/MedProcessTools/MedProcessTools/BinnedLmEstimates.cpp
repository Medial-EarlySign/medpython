#define _CRT_SECURE_NO_WARNINGS

#include "MedProcessTools/MedProcessTools/FeatureGenerator.h"
#include "MedProcessTools/MedProcessTools/MedProcessUtils.h"
#include "Logger/Logger/Logger.h"
#include <omp.h>

#define LOCAL_SECTION LOG_FTRGNRTR
#define LOCAL_LEVEL	LOG_DEF_LEVEL

#define BINNED_LM_MAX_AGE 110
#define BINNED_LM_MINIMAL_NUM_PER_AGE 100

//.......................................................................................
//.......................................................................................
// BinnedLinearModels : Apply a set of liner models to generate features
//.......................................................................................
//.......................................................................................

// Default Values
int def_bin_bounds[] = { -(30 * 9),-(30 * 4),-30,30,30 * 4,30 * 9,30 * 18,30 * 360 };
int def_nbin_bounds = sizeof(def_bin_bounds) / sizeof(int);

int def_estimation_points[] = { 90,180 };
int def_nestimation_points = sizeof(def_estimation_points) / sizeof(int);

int def_min_period = -(30 * 18);
int def_max_period = 30 * 9999;

double def_rfactor = 0.99;

// Initialization
//.......................................................................................
void BinnedLmEstimates::set_names() {

	if (names.empty()) {
		string base_name = "FTR_" + int_to_string_digits(serial_id, 6) + "." + signalName + ".Estimate.";
		for (int point : params.estimation_points) {
			string name = base_name + std::to_string(point);
			if (time_channel != 0 || val_channel != 0)
				name += ".t" + std::to_string(time_channel) + "v" + std::to_string(val_channel);
			names.push_back(name);
		}
	}

}

//.......................................................................................
void BinnedLmEstimates::set(string& _signalName) {

	signalName = _signalName;

	init_defaults();
	names.clear();
	set_names();

	req_signals.resize(3);
	req_signals[0] = "GENDER";
	req_signals[1] = "BDATE";
	req_signals[2] = signalName;
}

//.......................................................................................
void BinnedLmEstimates::init_defaults() {

	generator_type = FTR_GEN_BINNED_LM;

	params.bin_bounds.resize(def_nbin_bounds);
	for (int i = 0; i < def_nbin_bounds; i++)
		params.bin_bounds[i] = def_bin_bounds[i];

	params.max_period = def_max_period;
	params.min_period = def_min_period;
	params.rfactor = (float)def_rfactor;

	params.estimation_points.resize(def_nestimation_points);
	for (int i = 0; i < def_nestimation_points; i++)
		params.estimation_points[i] = def_estimation_points[i];

	req_signals.resize(2);
	req_signals[0] = "GENDER";
	req_signals[1] = "BDATE";

	signalId = -1;
	bdateId = -1;
	genderId = -1;

	time_unit_periods = global_default_windows_time_unit;
	sampling_strategy = BINNED_LM_TAKE_ALL;
}

//.......................................................................................
void BinnedLmEstimates::set(string& _signalName, BinnedLmEstimatesParams* _params) {

	signalName = _signalName;

	params.bin_bounds = _params->bin_bounds;
	params.max_period = _params->max_period;
	params.min_period = _params->min_period;
	params.rfactor = _params->rfactor;
	params.estimation_points = _params->estimation_points;

	set_names();

	req_signals.resize(3);
	req_signals[0] = "GENDER";
	req_signals[1] = "BDATE";
	req_signals[2] = signalName;
}

//..............................................................................
int BinnedLmEstimates::init(map<string, string>& mapper) {
	for (auto entry : mapper) {
		string field = entry.first;
		//! [BinnedLmEstimates::init]
		if (field == "bin_bounds") {
			if (init_dvec(entry.second, params.bin_bounds) == -1) {
				fprintf(stderr, "Cannot initialize bin_bounds for LM\n");
				return -1;
			}
		}
		else if (field == "max_period") params.max_period = stoi(entry.second);
		else if (field == "min_period") params.min_period = stoi(entry.second);
		else if (field == "rfactor") params.rfactor = stof(entry.second);
		else if (field == "signalName" || field == "signal") signalName = entry.second;
		else if (field == "estimation_points") {
			if (init_dvec(entry.second, params.estimation_points) == -1) {
				fprintf(stderr, "Cannot initialize estimation_points for LM\n");
				return -1;
			}
		}
		else if (field == "time_unit") time_unit_periods = med_time_converter.string_to_type(entry.second);
		else if (field == "time_channel") time_channel = stoi(entry.second);
		else if (field == "val_channel") val_channel = stoi(entry.second);
		else if (field == "tags") boost::split(tags, entry.second, boost::is_any_of(","));
		else if (field == "sampling") set_sampling_strategy(entry.second);
		else if (field == "weights_generator") iGenerateWeights = stoi(entry.second);
		else if (field != "fg_type")
			MLOG("Unknonw parameter \'%s\' for BinnedLmEstimates\n", field.c_str());
		//! [BinnedLmEstimates::init]
	}

	names.clear();
	set_names();

	req_signals.resize(3);
	req_signals[0] = "GENDER";
	req_signals[1] = "BDATE";
	req_signals[2] = signalName;


	return 0;
}

//..............................................................................
void BinnedLmEstimates::set_sampling_strategy(string& strategy) {

	boost::to_lower(strategy);

	if (strategy == "all" || strategy == "take_all")
		sampling_strategy = BINNED_LM_TAKE_ALL;
	else if (strategy == "first" || strategy == "stop_at_first")
		sampling_strategy = BINNED_LM_STOP_AT_FIRST;
	else if (strategy == "last" || strategy == "stop_at_last")
		sampling_strategy = BINNED_LM_STOP_AT_LAST;
	else
		MTHROW_AND_ERR("Unknonwn sampling strategy \'%s\' for BinnedLM\n", strategy.c_str());
}

//..............................................................................
void BinnedLmEstimates::set_signal_ids(MedSignals& sigs) {

	signalId = sigs.sid(signalName);
	genderId = sigs.sid("GENDER");
	bdateId = sigs.sid(req_signals[1]);
}

// Learn a generator
//.......................................................................................
int BinnedLmEstimates::_learn(MedPidRepository& rep, const MedSamples& samples, vector<RepProcessor *> processors) {

	// Sanity check
	if (signalId == -1 || genderId == -1 || bdateId == -1) {
		MERR("BinnedLmEstimates::_learn - Uninitialized signalId (%s)\n", signalName.c_str());
		return -1;
	}

	if (time_unit_sig == MedTime::Undefined)	time_unit_sig = rep.sigs.Sid2Info[signalId].time_unit;

	size_t nperiods = params.bin_bounds.size();
	size_t nmodels = 1ll << nperiods;
	size_t nfeatures = nperiods * (INT64_C(1) << nperiods);
	size_t nids = samples.idSamples.size();

	// Required signals
	vector<int> all_req_signal_ids_v;
	vector<unordered_set<int> > current_required_signal_ids(processors.size());
	vector<FeatureGenerator *> generators = { this };
	unordered_set<int> extra_req_signal_ids;
	handle_required_signals(processors, generators, extra_req_signal_ids, all_req_signal_ids_v, current_required_signal_ids);

	// Collect Data
	int nthreads = omp_get_max_threads();
	vector<PidDynamicRec> recs(nthreads);

	vector<float> values;
	vector<int> ages, times, genders;
	vector<int> id_firsts(nids), id_lasts(nids);

#pragma omp parallel for
	for (int i = 0; i < nids; i++) {
		UniversalSigVec usv, ageUsv;
		int byear, gender, age;
		int n_th = omp_get_thread_num();
		PidDynamicRec &rec = recs[n_th];
		int id = samples.idSamples[i].id;

		// Gender
		gender = medial::repository::get_value(rep, id, genderId);

		assert(gender != -1);

		// Get signal (if not virtual)
#pragma omp critical
		id_firsts[i] = (int)ages.size();
		if (rep.sigs.Sid2Info[signalId].virtual_sig == 0) {
			rep.uget(id, signalId, usv);

			if (usv.len == 0) {
#pragma omp critical
				id_lasts[i] = id_firsts[i];
				continue;
			}
		}

		// Get last time-point
		vector<int> time_points;
		int last_time_point = -1;
		if (sampling_strategy == BINNED_LM_STOP_AT_FIRST) {
			last_time_point = usv.Time(0, time_channel);
			time_points.push_back(last_time_point);
		}
		else if (sampling_strategy == BINNED_LM_STOP_AT_LAST) {
			last_time_point = usv.Time(usv.len - 1, time_channel);
			time_points.push_back(last_time_point);
		}

		int nvalues = 0;
		if (processors.size()) {
			rec.init_from_rep(std::addressof(rep), id, all_req_signal_ids_v, 1);

			// BYear/Age
			prepare_for_age(rec, ageUsv, age, byear);

			// Apply Processors
			vector<vector<float>> dummy_attributes_mat;
			for (unsigned int j = 0; j < processors.size(); ++j)
				processors[j]->conditional_apply(rec, time_points, current_required_signal_ids[j], dummy_attributes_mat);

			// Collect values and ages
			rec.uget(signalId, 0, usv);
#pragma omp critical
			for (int j = 0; j < usv.len; ++j) {
				if (sampling_strategy != BINNED_LM_TAKE_ALL && usv.Time(j, time_channel) > last_time_point)
					break;
				values.push_back(usv.Val(j, val_channel));
				get_age(usv.Time(j, time_channel), time_unit_sig, age, byear);
				ages.push_back(age);
				genders.push_back(gender);
				times.push_back(med_time_converter.convert_times(time_unit_sig, time_unit_periods, usv.Time(j, time_channel)));
				nvalues++;
			}
		}
		else {
			// BYear/Age
			prepare_for_age(rep, id, ageUsv, age, byear);
#pragma omp critical
			for (int j = 0; j < usv.len; j++) {
				if (sampling_strategy != BINNED_LM_TAKE_ALL && usv.Time(i, time_channel) > last_time_point)
					break;

				values.push_back(usv.Val(j, val_channel));
				get_age(usv.Time(j, time_channel), time_unit_sig, age, byear);
				ages.push_back(age);
				genders.push_back(gender);
				times.push_back(med_time_converter.convert_times(time_unit_sig, time_unit_periods, usv.Time(j, time_channel)));
				nvalues++;
			}
		}
#pragma omp critical
		id_lasts[i] = id_firsts[i] + nvalues - 1;
	}

	// Allocate
	int num = (int)values.size();
	if (num == 0) {
		MERR("No Data Collected for %s\n", signalName.c_str());
		return -1;
	}

	models.resize(nmodels);
	xmeans.resize(nfeatures, 0);
	xsdvs.resize(nfeatures, 0);
	ymeans.resize(nmodels, 0);
	ysdvs.resize(nmodels, 0);

	vector<double> sums[2];
	vector<int> nums[2];
	for (int igender = 0; igender < 2; igender++) {
		sums[igender].resize(BINNED_LM_MAX_AGE + 1, 0);
		nums[igender].resize(BINNED_LM_MAX_AGE + 1, 0);
		means[igender].resize(BINNED_LM_MAX_AGE + 1, 0);
	}

	// Gender/Age - Collect Data for means and standard deviations
	for (int i = 0; i < num; i++) {
		if (ages[i] <= BINNED_LM_MAX_AGE && ages[i] >= 1) {
			nums[genders[i] - 1][ages[i]] ++;
			sums[genders[i] - 1][ages[i]] += values[i];
		}
	}

	// Gender/Age - correct for missing data
	for (int igender = 0; igender < 2; igender++) {
		int most_common_age = 0;
		for (int iage = 1; iage <= BINNED_LM_MAX_AGE; iage++) {
			if (nums[igender][iage] > nums[igender][most_common_age])
				most_common_age = iage;
		}

		if (nums[igender][most_common_age] == 0) {
			MDBG(DEBUG_LOG_LEVEL, "No %s found for gender %d. Are we in a single gender mode ?\n", signalName.c_str(), igender + 1);
			continue;
		}
		else if (nums[igender][most_common_age] < BINNED_LM_MINIMAL_NUM_PER_AGE) {
			MERR("Not enough tests for %s (most common age = %d has only %d samples)\n", signalName.c_str(), most_common_age, nums[igender][most_common_age]);
			return -1;
		}

		for (int iage = most_common_age; iage <= BINNED_LM_MAX_AGE; iage++) {
			if (nums[igender][iage] < BINNED_LM_MINIMAL_NUM_PER_AGE) {
				int missing_num = BINNED_LM_MINIMAL_NUM_PER_AGE - nums[igender][iage];
				nums[igender][iage] = BINNED_LM_MINIMAL_NUM_PER_AGE;
				sums[igender][iage] += sums[igender][iage - 1] * ((0.0 + missing_num) / nums[igender][iage - 1]);
			}
			means[igender][iage] = (float)(sums[igender][iage] / nums[igender][iage]);
		}

		for (int iage = most_common_age - 1; iage >= 0; iage--) {
			if (nums[igender][iage] < BINNED_LM_MINIMAL_NUM_PER_AGE) {
				int missing_num = BINNED_LM_MINIMAL_NUM_PER_AGE - nums[igender][iage];
				nums[igender][iage] = BINNED_LM_MINIMAL_NUM_PER_AGE;
				sums[igender][iage] += sums[igender][iage + 1] * ((0.0 + missing_num) / nums[igender][iage + 1]);
			}

			means[igender][iage] = (float)(sums[igender][iage] / nums[igender][iage]);
		}
	}
	
	// Collect data
	MedMat<float> x((int)values.size(), (int)nperiods);
	vector<float> y(values.size());
	vector<int> types(values.size());
	int irow = 0;

	for (unsigned int i = 0; i < nids; i++) {

		if (id_lasts[i] <= id_firsts[i])
			continue;
		int gender = genders[id_firsts[i]];

		for (int idx1 = id_firsts[i]; idx1 <= id_lasts[i]; idx1++) {
			if (ages[idx1] > BINNED_LM_MAX_AGE || ages[idx1] < 1)
				continue;

			// Add line + type to data matrix
			int type = 0;
			int iperiod = (int)nperiods;
			int jperiod = (int)nperiods;

			float type_sum = 0;
			int type_num = 0;
			for (int idx2 = id_firsts[i]; idx2 <= id_lasts[i]; idx2++) {
				if (ages[idx2] > BINNED_LM_MAX_AGE || ages[idx2] < 1)
					continue;
				if (idx1 == idx2)
					continue;

				int gap = times[idx1] - times[idx2];
				if (gap < params.min_period || gap > params.max_period)
					continue;

				while (iperiod > 0 && gap <= params.bin_bounds[iperiod - 1])
					iperiod--;

				if (iperiod != jperiod) {
					if (type_num) {
						type += 1 << jperiod;
						x(irow, jperiod) = type_sum / type_num;
					}
					type_sum = 0;
					type_num = 0;
					jperiod = iperiod;
				}

				if (iperiod != nperiods) {
					type_sum += values[idx2] - means[gender - 1][ages[idx2]];
					type_num++;
				}
			}

			if (type_num) {
				type += 1 << jperiod;
				x(irow, jperiod) = type_sum / type_num;
			}

			y[irow] = values[idx1] - means[gender - 1][ages[idx1]];
			types[irow++] = type;
		}
	}
	
	// Build model for each class 
	// Collect Data
	int inrows = irow;

	models[0].n_ftrs = 0;
	for (int type = 1; type < nmodels; type++) {
		vector<int> cols(nperiods, 0);

		int ncols = 0;
		for (int iperiod = 0; iperiod < nperiods; iperiod++) {
			if (type & (1 << iperiod))
				cols[ncols++] = iperiod;
		}

		int jnrows = 0;
		for (int i = 0; i < inrows; i++) {
			if ((types[i] & type) == type)
				jnrows++;
		}

		if (jnrows < ncols) {
			fprintf(stderr, "Not enough samples of type %d (%d, required - %d) %s\n", type, jnrows, ncols, signalName.c_str());
			return -1;
		}

		MedMat<float> tx(ncols, jnrows);
		MedMat<float> ty(jnrows, 1);
		tx.transposed_flag = 1;

		int jrow = 0;
		for (int i = 0; i < inrows; i++) {
			if ((types[i] & type) == type) {
				ty(jrow, 0) = y[i];
				for (int j = 0; j < ncols; j++)
					tx(j, jrow) = x(i, cols[j]);
				jrow++;
			}
		}

		// Normalize
		tx.normalize(2);
		for (int j = 0; j < ncols; j++) {
			xmeans[nperiods*type + j] = tx.avg[j];
			xsdvs[nperiods*type + j] = tx.std[j];
		}

		ty.normalize();
		ymeans[type] = ty.avg[0];
		ysdvs[type] = ty.std[0];

		models[type].params.rfactor = params.rfactor;
		models[type].learn(tx, ty);
		//		models[type].print(stdout, "model." + signalName + "." + to_string(type));

	}

	return 0;
}

// generate new feature(s)
//.......................................................................................
int BinnedLmEstimates::_generate(PidDynamicRec& rec, MedFeatures& features, int index, int num, vector<float *> &_p_data) {

	// Sanity check
	if (signalId == -1 || genderId == -1 || bdateId == -1) {
		MERR("BinnedLmEstimates::_generate - Uninitialized signalId(%)\n", signalName.c_str());
		return -1;
	}

	if (time_unit_sig == MedTime::Undefined)	time_unit_sig = rec.my_base_rep->sigs.Sid2Info[signalId].time_unit;


	size_t nperiods = params.bin_bounds.size();
	size_t nfeatures = nperiods * (INT64_C(1) << nperiods);

	int iperiod = (int)nperiods;
	int jperiod = (int)nperiods;

	int byear, gender, age;

	// BYear/Age
	UniversalSigVec usv, ageUsv;
	prepare_for_age(rec, ageUsv, age, byear);

	// Gender
	gender = medial::repository::get_value(rec, genderId);

	if (means[gender - 1][0] == 0) {
		MERR("No age-dependent mean found for %s for gender %d\n", signalName.c_str(), gender);
		return -1;
	}

	// Features
	MedMat<float> x(1, (int)nfeatures);
	for (int i = 0; i < num; i++) {
		rec.uget(signalId, i);
		int last_time = med_time_converter.convert_times(features.time_unit, time_unit_periods, features.samples[index + i].time);
		int last_sig_time = med_time_converter.convert_times(features.time_unit, time_unit_sig, features.samples[index + i].time);

		for (unsigned int ipoint = 0; ipoint < params.estimation_points.size(); ipoint++) {
			int type = 0;
			float type_sum = 0;
			int type_num = 0;

			float *p_feat = _p_data[ipoint] + index;
			int days_diff = last_time - params.estimation_points[ipoint];
			if (days_diff < 0) {
				MWARN("WARN :: BinnedLmEstimates::_generate - got date for pid %d and signal %s with time window before 1900\n",
					rec.pid, signalName.c_str());
				days_diff = 0;
			}
			int target_time = med_time_converter.convert_times(time_unit_periods, time_unit_sig, days_diff);

			for (int j = 0; j < rec.usv.len; j++) {
				int time = rec.usv.Time(j, time_channel);
				if (time > last_sig_time)
					break;

				int gap = target_time - time;
				if (gap < params.min_period || gap > params.max_period)
					continue;

				while (iperiod > 0 && gap <= params.bin_bounds[iperiod - 1])
					iperiod--;

				if (iperiod != jperiod) {
					if (type_num) {
						type += 1 << jperiod;
						x(0, jperiod) = type_sum / type_num;
					}
					type_sum = 0;
					type_num = 0;
					jperiod = iperiod;
				}

				if (iperiod != nperiods) {
					get_age(rec.usv.Time(j, time_channel), time_unit_sig, age, byear);
					if (age > BINNED_LM_MAX_AGE) age = BINNED_LM_MAX_AGE;
					if (age < 1) age = 1;
					type_sum += rec.usv.Val(j, val_channel) - means[gender - 1][age];
					type_num++;
				}
			}

			if (type_num) {
				type += 1 << jperiod;
				x(0, jperiod) = type_sum / type_num;
			}

			get_age(last_time, time_unit_periods, age, byear);
			if (age > BINNED_LM_MAX_AGE) age = BINNED_LM_MAX_AGE;
			if (age < 1) age = 1;

			// Predict
			if (type) {
				float pred = 0;
				int j = 0;
				for (iperiod = 0; iperiod < nperiods; iperiod++) {
					if (type & (1 << iperiod)) {
						pred += models[type].b[j] * (x(0, iperiod) - xmeans[nperiods*type + j]) / xsdvs[nperiods*type + j];
						j++;
					}
				}

				p_feat[i] = (pred + ymeans[type] + means[gender - 1][age]);
			}
			else {
				p_feat[i] = missing_val;
			}
		}
	}

	return 0;
}

// Get pointers to data vectors
//.......................................................................................
void BinnedLmEstimates::get_p_data(MedFeatures &features, vector<float *> &_p_data) {

	p_data.clear();

	if (iGenerateWeights) {
		if (names.size() != 1)
			MTHROW_AND_ERR("Cannot generate weights using a multi-feature generator (type %d generates %d features)\n", generator_type, (int)names.size())
		else
			_p_data[0] = &(features.weights[0]);
	}
	else {
		for (unsigned int ipoint = 0; ipoint < params.estimation_points.size(); ipoint++)
			_p_data.push_back(&(features.data[names[ipoint]][0]));
	}

	return;
}

//.......................................................................................
// Filter generated features according to a set. return number of valid features (does not affect single-feature genertors, just returns 1/0 if feature name in set)
int BinnedLmEstimates::filter_features(unordered_set<string>& validFeatures) {

	vector<int> e_points;
	vector<string> names_new;
	for (int i = 0; i < names.size(); i++) {
		if (validFeatures.find(names[i]) != validFeatures.end()) {
			e_points.push_back(params.estimation_points[i]);
			names_new.push_back(names[i]);
		}
	}

	names = move(names_new);
	params.estimation_points = move(e_points);

	return ((int)names.size());
}

// Age Related functions
//.......................................................................................
void BinnedLmEstimates::prepare_for_age(PidDynamicRec& rec, UniversalSigVec& ageUsv, int &age, int &byear) {
	int bdate = medial::repository::get_value(rec, bdateId);
	assert(bdate != -1);
	if (req_signals[1] == "BDATE")
		byear = int(bdate / 10000);
	else
		byear = bdate;
}

void BinnedLmEstimates::prepare_for_age(MedPidRepository& rep, int id, UniversalSigVec& ageUsv, int &age, int &byear) {
	int bdate = medial::repository::get_value(rep, id, bdateId);
	assert(bdate != -1);
	if (req_signals[1] == "BDATE")
		byear = int(bdate / 10000);
	else
		byear = bdate;
}

inline void BinnedLmEstimates::get_age(int time, int time_unit_from, int& age, int byear) {
	age = med_time_converter.convert_times(time_unit_from, MedTime::Date, time) / 10000 - byear;
}

void BinnedLmEstimates::dprint(const string &pref, int fg_flag) {
	if (fg_flag > 0) {
		MLOG("%s :: BinnedLmEstimates : serial_id %d : ", pref.c_str(), serial_id);
		MLOG("names(%d) : ", names.size());
		if (fg_flag > 1) for (auto &name : names) MLOG("%s,", name.c_str());
		MLOG(" tags(%d) : ", tags.size());
		if (fg_flag > 1) for (auto &t : tags) MLOG("%s,", t.c_str());
		MLOG(" req_signals(%d) : ", req_signals.size());
		if (fg_flag > 1) for (auto &rsig : req_signals) MLOG("%s,", rsig.c_str());
		string s = "";
		for (size_t i = 0; i < params.estimation_points.size(); ++i)
			s += to_string(params.estimation_points[i]) + ", ";
		MLOG(" params.period=[%d, %d], params.estimation_points(%zu):[%s]",
			params.min_period, params.max_period, params.estimation_points.size(), s.c_str());
		MLOG("\n");
	}
}

// Print predictor
//.......................................................................................
void BinnedLmEstimates::print()
{
	string prefix = "BinnedLmEstimates(" + signalName + ") :: ";

	string sout = "";

	sout += prefix + "nmodels:" + to_string(models.size()) + "\n";

	for (int i = 0; i < models.size(); i++) {
		string mprefix = prefix + "model " + to_string(i) + " :: ";
		sout += mprefix + "len: " + to_string(models[i].b.size()) + " b0: " + to_string(models[i].b0) + " b:: ";
		for (int j = 0; j < models[i].b.size(); j++)
			sout += to_string(j) + ":" + to_string(models[i].b[j]) + " ";
		sout += "\n";
	}

	sout += prefix + "xmeans (" + to_string(xmeans.size()) + ")" + " ::";
	for (int i = 0; i < xmeans.size(); i++) {
		sout += to_string(i) + ":" + to_string(xmeans[i]) + " ";
	}
	sout += "\n";

	sout += prefix + "ymeans (" + to_string(ymeans.size()) + ")" + " ::";
	for (int i = 0; i < ymeans.size(); i++) {
		sout += to_string(i) + ":" + to_string(ymeans[i]) + " ";
	}
	sout += "\n";

	sout += prefix + "xsdvs (" + to_string(xsdvs.size()) + ")" + " ::";
	for (int i = 0; i < xsdvs.size(); i++) {
		sout += to_string(i) + ":" + to_string(xsdvs[i]) + " ";
	}
	sout += "\n";

	sout += prefix + "ysdvs (" + to_string(ysdvs.size()) + ")" + " ::";
	for (int i = 0; i < ysdvs.size(); i++) {
		sout += to_string(i) + ":" + to_string(ysdvs[i]) + " ";
	}
	sout += "\n";

	sout += prefix + "means[0] (" + to_string(means[0].size()) + ")" + " ::";
	for (int i = 0; i < means[0].size(); i++) {
		sout += to_string(i) + ":" + to_string(means[0][i]) + " ";
	}
	sout += "\n";

	sout += prefix + "means[1] (" + to_string(means[1].size()) + ")" + " ::";
	for (int i = 0; i < means[1].size(); i++) {
		sout += to_string(i) + ":" + to_string(means[1][i]) + " ";
	}
	sout += "\n";

	MLOG("%s", sout.c_str());
}

