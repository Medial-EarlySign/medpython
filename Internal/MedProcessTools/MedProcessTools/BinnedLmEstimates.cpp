#define _CRT_SECURE_NO_WARNINGS

#include "MedProcessTools/MedProcessTools/FeatureGenerator.h"
#include "MedProcessTools/MedProcessTools/MedProcessUtils.h"
#include "Logger/Logger/Logger.h"

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

//.......................................................................................
void BinnedLmEstimates::set_names() {

	if (names.empty()) {
		string base_name = signalName + ".Estimate.";
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
	req_signals[1] = "BYEAR";
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

	req_signals ={ "GENDER", "BYEAR" };

	signalId = -1; 
	byearId = -1;
	genderId = -1;
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

	req_signals ={ "GENDER", "BYEAR", signalName };

}

//..............................................................................
int BinnedLmEstimates::init(map<string, string>& mapper) {
	init_defaults();

	for (auto entry : mapper) {
		string field = entry.first;

		if (field == "bin_bounds") {
			if (init_dvec(entry.second, params.bin_bounds) == -1) {
				fprintf(stderr, "Cannot initialize bin_bounds for LM\n");
				return -1;
			}
		}
		else if (field == "max_period") params.max_period = stoi(entry.second);
		else if (field == "min_period") params.min_period = stoi(entry.second);
		else if (field == "rfactor") params.rfactor = stof(entry.second);
		else if (field == "signalName") signalName = entry.second;
		else if (field == "estimation_points") {
			if (init_dvec(entry.second, params.estimation_points) == -1) {
				fprintf(stderr, "Cannot initialize estimation_points for LM\n");
				return -1;
			}
		}
		else MLOG("Unknonw parameter \'%s\' for BinnedLmEstimates\n", field.c_str());
	}

	names.clear();
	set_names();

	req_signals = { "GENDER", "BYEAR", signalName };

	return 0;
}

// Learn a generator
//.......................................................................................
int BinnedLmEstimates::_learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *> processors) {

	// Sanity check
	if (signalId == -1) {
		MERR("Uninitialized signalId\n");
		return -1;
	}
	assert(rep.sigs.type(signalId) == T_DateVal);

	if (byearId == -1)
		byearId = rep.dict.id("BYEAR");

	if (genderId == -1)
		genderId = rep.dict.id("GENDER");

	size_t nperiods = params.bin_bounds.size();
	size_t nmodels = 1 << nperiods;
	size_t nfeatures = nperiods * (INT64_C(1) << nperiods);

	// Collect Data
	int len, byear, gender;
	PidDynamicRec rec;

	vector<float> values;
	vector<int> ages, days, genders;
	vector<int> id_firsts(ids.size()), id_lasts(ids.size());

	for (unsigned int i = 0; i < ids.size(); i++) {
		int id = ids[i];

		// BYear
		SVal *bYearSignal = (SVal *)rep.get(id, byearId, len);
		assert(len == 1);
		byear = (int)(bYearSignal[0].val);

		// Gender
		SVal *genderSignal = (SVal *)rep.get(id, genderId, len);
		assert(len == 1);
		gender = (int)(genderSignal[0].val);

		// Get signal
		SDateVal *signal = (SDateVal *)rep.get(id, signalId, len);
		id_firsts[i] = (int) ages.size();

		if (len == 0) {
			id_lasts[i] = id_firsts[i];
			continue;
		}

		if (processors.size()) {

			// Apply processing at last time point only
			vector<int> time_points(1, signal[len-1].date + 1);

			// Init Dynamic Rec
			rec.init_from_rep(std::addressof(rep), id, req_signal_ids, 1);

			// Apply
			for (auto& processor : processors)
				processor->apply(rec, time_points, req_signal_ids);

			// Collect values and ages
			SDateVal *clnSignal = (SDateVal *)rec.get(signalId, 0, len);
			for (int i = 0; i < len; i++) {
				values.push_back(clnSignal[i].val);
				ages.push_back(clnSignal[i].date / 10000 - byear);
				genders.push_back(gender);
				days.push_back(get_day(clnSignal[i].date));
			}
		}
		else {
			for (int i = 0; i < len; i++) {
				values.push_back(signal[i].val);
				ages.push_back(signal[i].date / 10000 - byear);
				genders.push_back(gender);
				days.push_back(get_day(signal[i].date));
			}
		}
		id_lasts[i] = id_firsts[i] + len - 1;
	}

	// Allocate
	int num = (int) values.size();
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
		if (ages[i] <= BINNED_LM_MAX_AGE) {
			nums[genders[i]-1][ages[i]] ++;
			sums[genders[i]-1][ages[i]] += values[i];
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
			MDBG(DEBUG_LOG_LEVEL,"No %s found for gender %d. Are we in a single gender mode ?\n", signalName.c_str(), igender + 1);
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
	MedMat<float> x((int) values.size(), (int) nperiods);
	vector<float> y(values.size());
	vector<int> types(values.size());
	int irow = 0;

	for (unsigned int i = 0; i < ids.size(); i++) {
		if (id_lasts[i] <= id_firsts[i])
			continue;
		int gender = genders[id_firsts[i]];

		for (int idx1 = id_firsts[i]; idx1 <= id_lasts[i]; idx1++) {
			if (ages[idx1] > BINNED_LM_MAX_AGE)
				continue;

			// Add line + type to data matrix
			int type = 0;
			int iperiod = (int) nperiods;
			int jperiod = (int) nperiods;

			float type_sum = 0;
			int type_num = 0;
			for (int idx2 = id_firsts[i]; idx2 <= id_lasts[i]; idx2++) {
				if (idx1 == idx2)
					continue;

				int gap = days[idx1] - days[idx2];
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
					type_sum += values[idx2] - means[gender-1][ages[idx2]];
					type_num++;
				}
			}

			if (type_num) {
				type += 1 << jperiod;
				x(irow, jperiod) = type_sum / type_num;
			}

			y[irow] = values[idx1] - means[gender-1][ages[idx1]];
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
			fprintf(stderr, "Not enough samples of type %d (%d, required - %d)\n", type, jnrows, ncols);
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
//		models[type].print(stderr, "model" + to_string(type));

	}


	return 0;
}

// generate new feature(s)
//.......................................................................................
int BinnedLmEstimates::Generate(PidDynamicRec& rec, MedFeatures& features, int index, int num) {

	// Sanity check
	if (signalId == -1) {
		MERR("Uninitialized signalId\n");
		return -1;
	}
	assert(rec.my_base_rep->sigs.type(signalId) == T_DateVal);

	if (byearId == -1)
		byearId = rec.my_base_rep->dict.id("BYEAR");

	if (genderId == -1)
		genderId = rec.my_base_rep->dict.id("GENDER");

	size_t nperiods = params.bin_bounds.size();
	size_t nmodels = 1 << nperiods;
	size_t nfeatures = nperiods * (INT64_C(1) << nperiods);

	int iperiod = (int) nperiods;
	int jperiod = (int) nperiods;
	
	int len, byear, gender;

	vector<float> values;
	vector<int> ages, days;

	// BYear
	SVal *bYearSignal = (SVal *)rec.get(byearId, len);
	assert(len == 1);
	byear = (int)(bYearSignal[0].val);

	// Gender
	SVal *genderSignal = (SVal *)rec.get(genderId, len);
	assert(len == 1);
	gender = (int)(genderSignal[0].val);

	if (means[gender - 1][0] == 0) {
		MERR("No age-dependent mean found for %s for gender %d\n", signalName.c_str(), gender);
		return -1;
	}

	// Features
	MedMat<float> x(1, (int) nfeatures);
	for (int i = 0; i < num; i++) {
		SDateVal *signal = (SDateVal *)rec.get(signalId, i, len);
		int last_date = features.samples[index+i].date;
		int last_day = get_day(last_date);

		for (unsigned int ipoint = 0; ipoint < params.estimation_points.size(); ipoint++) {
			int type = 0;
			float type_sum = 0;
			int type_num = 0;

			float *p_feat = &(features.data[names[0]][index]);
			int target_day = last_day - params.estimation_points[ipoint];

			for (int j = 0; j < len; j++) {
				if (signal[j].date > last_date)
					break;

				int gap = target_day - get_day(signal[j].date);
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
					int age = (int)signal[j].date / 10000 - byear;
					if (age > BINNED_LM_MAX_AGE)
						age = BINNED_LM_MAX_AGE;
					type_sum += signal[j].val - means[gender-1][age];
					type_num++;
				}
			}

			if (type_num) {
				type += 1 << jperiod;
				x(0, jperiod) = type_sum / type_num;
			}

			int age = last_date / 10000 - byear - params.estimation_points[ipoint] / 365;
			if (age > BINNED_LM_MAX_AGE)
				age = BINNED_LM_MAX_AGE;

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

				features.data[names[ipoint]][index + i] = (pred + ymeans[type] + means[gender-1][age]);
			}
			else {
				features.data[names[ipoint]][index + i] = missing_val;
			}
		}
	}

	return 0;
}

// (De)Serialization
//.......................................................................................
size_t BinnedLmEstimates::get_size() {

	size_t size = 0;

	// signalName
	size += sizeof(size_t);
	size += signalName.length() + 1;

	// Params
	size += sizeof(size_t); size += params.bin_bounds.size() * sizeof(int);
	size += 2 * sizeof(int);
	size += sizeof(float);
	size += sizeof(size_t); 
	size += params.estimation_points.size() * sizeof(int);

	size_t nperiods = params.bin_bounds.size();
	size_t nmodels = 1 << nperiods;
	size_t nfeatures = nperiods * (INT64_C(1) << nperiods);

	// Means and Sdvs
	size += (2 * nfeatures + nmodels) * sizeof(float);
	size += 2 * (BINNED_LM_MAX_AGE+1) * sizeof(float);

	// Models
	for (auto& model : models)
		size += model.get_size();

	return size;

}

extern char signalName_c[MAX_NAME_LEN + 1];

//.......................................................................................
size_t BinnedLmEstimates::serialize(unsigned char *blob) {

	size_t ptr = 0;

	// SignalName
	size_t nameLen = signalName.length();
	assert(nameLen < MAX_NAME_LEN);

	strcpy(signalName_c, signalName.c_str());

	memcpy(blob + ptr, &nameLen, sizeof(size_t)); ptr += sizeof(size_t);
	memcpy(blob + ptr, signalName_c, nameLen + 1); ptr += nameLen + 1;

	size_t npoints = params.estimation_points.size();
	size_t nperiods = params.bin_bounds.size();
	size_t nmodels = 1 << nperiods;
	size_t nfeatures = nperiods * nmodels;

	// Params
	memcpy(blob + ptr, &nperiods, sizeof(size_t));  ptr += sizeof(size_t);
	memcpy(blob + ptr, &(params.bin_bounds[0]), nperiods * sizeof(int)); ptr += nperiods * sizeof(int);
	memcpy(blob + ptr, &params.min_period, sizeof(int)); ptr += sizeof(int);
	memcpy(blob + ptr, &params.max_period, sizeof(int)); ptr += sizeof(int);
	memcpy(blob + ptr, &params.rfactor, sizeof(float)); ptr += sizeof(float);
	memcpy(blob + ptr, &npoints, sizeof(size_t));  ptr += sizeof(size_t);
	memcpy(blob + ptr, &(params.estimation_points[0]), npoints * sizeof(int)); ptr += npoints * sizeof(int);

	// Means and Sdvs
	memcpy(blob + ptr, &(xmeans[0]), nfeatures*sizeof(float)); ptr += nfeatures*sizeof(float);
	memcpy(blob + ptr, &(xsdvs[0]), nfeatures*sizeof(float)); ptr += nfeatures*sizeof(float);
	memcpy(blob + ptr, &(ymeans[0]), nmodels*sizeof(float)); ptr += nmodels*sizeof(float);
	memcpy(blob + ptr, &(means[0][0]), (BINNED_LM_MAX_AGE + 1) * sizeof(float)); ptr += (BINNED_LM_MAX_AGE + 1)*sizeof(float);
	memcpy(blob + ptr, &(means[1][0]), (BINNED_LM_MAX_AGE + 1) * sizeof(float)); ptr += (BINNED_LM_MAX_AGE + 1)*sizeof(float);

	// Models
	for (auto& model : models) 
		ptr += model.serialize(blob + ptr);

	return ptr;
}

//.......................................................................................
size_t BinnedLmEstimates::deserialize(unsigned char *blob) {

	size_t ptr = 0;

	size_t npoints, nperiods;
	
	// SignalName
	size_t nameLen;
	memcpy(&nameLen, blob + ptr, sizeof(size_t)); ptr += sizeof(size_t);
	assert(nameLen < MAX_NAME_LEN);

	memcpy(signalName_c, blob + ptr, nameLen + 1); ptr += nameLen + 1;
	signalName = signalName_c;

	// Params
	memcpy(&nperiods, blob + ptr, sizeof(size_t));  ptr += sizeof(size_t);
	params.bin_bounds.resize(nperiods);
	memcpy(&(params.bin_bounds[0]), blob + ptr, nperiods * sizeof(int)); ptr += nperiods * sizeof(int);
	memcpy(&params.min_period, blob + ptr, sizeof(int)); ptr += sizeof(int);
	memcpy(&params.max_period, blob + ptr, sizeof(int)); ptr += sizeof(int);
	memcpy(&params.rfactor, blob + ptr, sizeof(float)); ptr += sizeof(float);
	memcpy(&npoints, blob + ptr, sizeof(size_t));  ptr += sizeof(size_t);
	params.estimation_points.resize(npoints);
	memcpy(&(params.estimation_points[0]), blob + ptr, npoints * sizeof(int)); ptr += npoints * sizeof(int);

	size_t nmodels = 1 << nperiods;
	size_t nfeatures = nperiods *nmodels;

	// Means and Sdvs
	xmeans.resize(nfeatures*sizeof(float));
	memcpy(&(xmeans[0]), blob + ptr, nfeatures*sizeof(float)); ptr += nfeatures*sizeof(float);
	xsdvs.resize(nfeatures*sizeof(float));
	memcpy(&(xsdvs[0]), blob + ptr, nfeatures*sizeof(float)); ptr += nfeatures*sizeof(float);
	ymeans.resize(nmodels*sizeof(float));
	memcpy(&(ymeans[0]), blob + ptr, nmodels*sizeof(float)); ptr += nmodels*sizeof(float);

	means[0].resize(BINNED_LM_MAX_AGE + 1);
	memcpy(&(means[0][0]), blob + ptr, (BINNED_LM_MAX_AGE + 1)*sizeof(float)); ptr += (BINNED_LM_MAX_AGE + 1)*sizeof(float);

	means[1].resize(BINNED_LM_MAX_AGE + 1);
	memcpy(&(means[1][0]), blob + ptr, (BINNED_LM_MAX_AGE + 1)*sizeof(float)); ptr += (BINNED_LM_MAX_AGE + 1)*sizeof(float);

	// Models
	models.resize(nmodels);
	for (auto& model : models) 
		ptr += model.deserialize(blob + ptr);

	set_names();

	req_signals.resize(3);
	req_signals[0] = "GENDER";
	req_signals[1] = "BYEAR";
	req_signals[2] = signalName;

	return ptr;
}