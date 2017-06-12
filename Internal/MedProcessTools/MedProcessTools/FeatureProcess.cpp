#define _CRT_SECURE_NO_WARNINGS

#define LOCAL_SECTION LOG_FEATCLEANER
#define LOCAL_LEVEL	LOG_DEF_LEVEL

#include "FeatureProcess.h"
#include "DoCalcFeatProcessor.h"
#include <omp.h>

//=======================================================================================
// Feature Processors
//=======================================================================================
// Processor types
FeatureProcessorTypes feature_processor_name_to_type(const string& processor_name) 
{
	if (processor_name == "multi_processor" || processor_name == "multi")
		return FTR_PROCESS_MULTI;
	else if (processor_name == "basic_outlier_cleaner" || processor_name == "basic_cleaner" || processor_name == "basic_cln")
		return FTR_PROCESS_BASIC_OUTLIER_CLEANER;
	else if (processor_name == "normalizer")
		return FTR_PROCESS_NORMALIZER;
	else if (processor_name == "imputer")
		return FTR_PROCESS_IMPUTER;
	else if (processor_name == "do_calc")
		return FTR_PROCESS_DO_CALC;
	else if (processor_name == "univariate_selector")
		return FTR_PROCESS_UNIVARIATE_SELECTOR;
	else
		return FTR_PROCESS_LAST;
}

// Initialization
//.......................................................................................
FeatureProcessor* FeatureProcessor::make_processor(string processor_name) {

	return make_processor(feature_processor_name_to_type(processor_name));
}

//.......................................................................................
FeatureProcessor * FeatureProcessor::make_processor(string processor_name, string init_string) {

	FeatureProcessorTypes type = feature_processor_name_to_type(processor_name);
	return make_processor(type , init_string);
}

//.......................................................................................
FeatureProcessor * FeatureProcessor::make_processor(FeatureProcessorTypes processor_type) {

	if (processor_type == FTR_PROCESS_MULTI)
		return new MultiFeatureProcessor;
	else if (processor_type == FTR_PROCESS_BASIC_OUTLIER_CLEANER)
		return new FeatureBasicOutlierCleaner;
	else if (processor_type == FTR_PROCESS_NORMALIZER)
		return new FeatureNormalizer;
	else if (processor_type == FTR_PROCESS_IMPUTER)
		return new FeatureImputer;
	else if (processor_type == FTR_PROCESS_DO_CALC)
		return new DoCalcFeatProcessor;
	else if (processor_type == FTR_PROCESS_UNIVARIATE_SELECTOR)
		return new UnivariateFeatureSelector;
	else
		return NULL;

}

//.......................................................................................
FeatureProcessor * FeatureProcessor::make_processor(FeatureProcessorTypes processor_type, string init_string) {

	FeatureProcessor *newProcessor = make_processor(processor_type);
	newProcessor->init_from_string(init_string);
	return newProcessor;
}

//.......................................................................................
int FeatureProcessor::learn(MedFeatures& features) {

	// All Ids - mark as an empty set
	unordered_set<int> temp;
	return Learn(features, temp);

}

//.......................................................................................
int FeatureProcessor::apply(MedFeatures& features) {

	// All Ids - mark as an empty set
	unordered_set<int> temp;
	return Apply(features, temp);
}

// (De)Serialize
//.......................................................................................
size_t FeatureProcessor::get_processor_size() {
	return sizeof(processor_type) + get_size();
}

//.......................................................................................
size_t FeatureProcessor::processor_serialize(unsigned char *blob) {

	size_t ptr = 0;
	memcpy(blob + ptr, &processor_type, sizeof(FeatureProcessorTypes)); ptr += sizeof(FeatureProcessorTypes);
	ptr += serialize(blob + ptr);

	return ptr;
}

//=======================================================================================
// MultiFeatureProcessor
//=======================================================================================
int MultiFeatureProcessor::Learn(MedFeatures& features, unordered_set<int>& ids) {

	// Create processors
	if (processors.size() == 0) {
		vector<string> features_to_process;
		for (auto& rec : features.data)
			features_to_process.push_back(rec.first);
		add_processors_set(members_type, features_to_process, init_string);
	}

	//	for (auto& cleaner : cleaners) {
	int RC = 0;
#pragma omp parallel for schedule(dynamic)
	for (int j = 0; j<processors.size(); j++) {
		int rc = processors[j]->Learn(features, ids);
//#pragma omp critical
		if (rc < 0) RC = -1;
	}

	return RC;
}

//.......................................................................................
int MultiFeatureProcessor::Apply(MedFeatures& features, unordered_set<int>& ids) {

	int RC = 0;
#pragma omp parallel for schedule(dynamic)
	for (int j = 0; j<processors.size(); j++) {
		int rc = processors[j]->Apply(features, ids);
//#pragma omp critical
		if (rc < 0) RC = -1;
	}

	return RC;
}
void MultiFeatureProcessor::get_feature_names(vector<string>& all_feature_names) {
	all_feature_names.clear();
	for (auto p : processors) {
		vector<string> my_feature_names;
		p->get_feature_names(my_feature_names);
		all_feature_names.insert(all_feature_names.end(), my_feature_names.begin(), my_feature_names.end());
	}
}

// Add processors
//.......................................................................................
void  MultiFeatureProcessor::add_processors_set(FeatureProcessorTypes type, vector<string>& features) {
	
	for (string& feature : features) {
		FeatureProcessor *processor = FeatureProcessor::make_processor(type);
		processor->set_feature_name(feature);
		processors.push_back(processor);
	}
}

void  MultiFeatureProcessor::add_processors_set(FeatureProcessorTypes type, vector<string>& features, string init_string) {
	
	for (string& feature : features) {
		FeatureProcessor *processor = FeatureProcessor::make_processor(type, init_string);
		processor->set_feature_name(feature);
		processors.push_back(processor);
	}

}

// Filter according to a subset of features
//.......................................................................................
int  MultiFeatureProcessor::filter(unordered_set<string>& features) {

	int idx = 0;
	for (int i = 0; i < processors.size(); i++) {
		if (features.find(processors[i]->feature_name) != features.end())
			processors[idx++] = processors[i];
	}

	processors.resize(idx);
	return (int)processors.size();

}

// Copy
//.......................................................................................
void MultiFeatureProcessor::copy(FeatureProcessor *processor) {

	MultiFeatureProcessor *tempProcessor = dynamic_cast<MultiFeatureProcessor*>(processor);
	assert(tempProcessor != 0);

	*this = *tempProcessor;

	processors.resize(tempProcessor->processors.size());
	for (int i = 0; i < processors.size(); i++) {
		processors[i] = make_processor(tempProcessor->processors[i]->processor_type);
		processors[i]->copy(tempProcessor->processors[i]);
	}
}

// Serialization
//.......................................................................................
size_t MultiFeatureProcessor::get_size() {

	// Number of processors
	size_t size = sizeof(int);

	for (auto& processor : processors)
		size += processor->get_processor_size();

	return size;
}

//.......................................................................................
size_t MultiFeatureProcessor::serialize(unsigned char *blob) {

	size_t ptr = 0;

	int nProcessors = (int)processors.size();
	memcpy(blob + ptr, &nProcessors, sizeof(int)); ptr += sizeof(int);

	for (auto& processor : processors)
		ptr += processor->processor_serialize(blob + ptr);

	return ptr;
}

//.......................................................................................
size_t MultiFeatureProcessor::deserialize(unsigned char *blob) {

	size_t ptr = 0;

	// number of cleaners
	int nProcessors;
	memcpy(&nProcessors, blob + ptr, sizeof(int)); ptr += sizeof(int);
	processors.resize(nProcessors);

	processors.resize(nProcessors);
	for (int i = 0; i < nProcessors; i++) {
		FeatureProcessorTypes type;
		memcpy(&type, blob + ptr, sizeof(FeatureProcessorTypes)); ptr += sizeof(FeatureProcessorTypes);
		processors[i] = FeatureProcessor::make_processor(type);
		ptr += processors[i]->deserialize(blob + ptr);
	}

	return ptr;
}

//=======================================================================================
// FeatureBasicOutlierCleaner
//=======================================================================================
//.......................................................................................
int FeatureBasicOutlierCleaner::Learn(MedFeatures& features, unordered_set<int>& ids) {

	if (params.type == VAL_CLNR_ITERATIVE)
		return iterativeLearn(features, ids);
	else if (params.type == VAL_CLNR_QUANTILE)
		return quantileLearn(features, ids);
	else {
		MERR("Unknown cleaning type %d\n", params.type);
		return -1;
	}
}

//.......................................................................................
int FeatureBasicOutlierCleaner::iterativeLearn(MedFeatures& features, unordered_set<int>& ids) {

	// Sanity
	if (features.data.find(feature_name) == features.data.end()) {
		MERR("Cannot find signal %s in features data\n", feature_name.c_str());
		return -1;
	}

	// Get all values
	vector<float> values;
	get_all_values(features, feature_name, ids, values);

	// Get bounds
	if (values.size() == 0)
		MWARN("EMPTY_VECTOR:: feature [%s] has 0 values\n", feature_name.c_str());

	int rc = get_iterative_min_max(values);
	return rc;
}

//.......................................................................................
int FeatureBasicOutlierCleaner::quantileLearn(MedFeatures& features, unordered_set<int>& ids) {

	// Sanity
	if (features.data.find(feature_name) == features.data.end()) {
		MERR("Cannot find signal %s in features data\n", feature_name.c_str());
		return -1;
	}

	// Get all values
	vector<float> values;
	get_all_values(features, feature_name, ids, values);

	// Get bounds
	return get_quantile_min_max(values);
}

// Clean
//.......................................................................................
int FeatureBasicOutlierCleaner::Apply(MedFeatures& features, unordered_set<int>& ids) {

	// Sanity
	if (features.data.find(feature_name) == features.data.end()) {
		MERR("Cannot find signal %s in features data\n", feature_name.c_str());
		return -1;
	}

	// Clean
	bool empty = ids.empty();
	vector<float>& data = features.data[feature_name];
	for (unsigned int i = 0; i < features.samples.size(); i++) {
		if ((empty || ids.find(features.samples[i].id) != ids.end()) && data[i] != params.missing_value) {
			if (params.doRemove && (data[i] < removeMin - NUMERICAL_CORRECTION_EPS || data[i] > removeMax + NUMERICAL_CORRECTION_EPS))
				data[i] = params.missing_value;
			else if (params.doTrim) {
				if (data[i] < trimMin)
					data[i] = trimMin;
				else if (data[i] > trimMax)
					data[i] = trimMax;
			}
		}
	}
	return 0;
}

// (De)Serialization
//.......................................................................................
size_t FeatureBasicOutlierCleaner::get_size() {
	return MedSerialize::get_size(processor_type, feature_name, params.doTrim, params.doRemove, trimMax, trimMin, removeMax, removeMin);
}

//.......................................................................................
size_t FeatureBasicOutlierCleaner::serialize(unsigned char *blob) {
	return MedSerialize::serialize(blob, processor_type, feature_name, params.doTrim, params.doRemove, trimMax, trimMin, removeMax, removeMin);
}

//.......................................................................................
size_t FeatureBasicOutlierCleaner::deserialize(unsigned char *blob) {
	return MedSerialize::deserialize(blob, processor_type, feature_name, params.doTrim, params.doRemove, trimMax, trimMin, removeMax, removeMin);
}

//=======================================================================================
// FeatureNormalizer
//=======================================================================================
//.......................................................................................
int FeatureNormalizer::Learn(MedFeatures& features, unordered_set<int>& ids) {

	// Get all values
	vector<float> values;
	get_all_values(features, feature_name, ids, values);

	vector<float> wgts(values.size(), 1.0);
	int rc = get_moments(values, wgts, missing_value, mean, sd);
	if (sd == 1) {
		MLOG("got sd=1.0 in feature %s....\n", feature_name.c_str());
	}

	if (sd == 0)
		MTHROW_AND_ERR("FeatureNormalizer learn sd: %f mean: %f size: %d", sd, mean, (int)values.size());
	return rc;
}

// Apply
//.......................................................................................
int FeatureNormalizer::Apply(MedFeatures& features, unordered_set<int>& ids) {

	// Sanity
	if (features.data.find(feature_name) == features.data.end()) {
		MERR("Cannot find signal %s in features data\n", feature_name.c_str());
		return -1;
	}

	// Attribute
	features.attributes[feature_name].normalized = true;
	if (fillMissing)
		features.attributes[feature_name].imputed = true;

	// Clean
	bool empty = ids.empty();
	vector<float>& data = features.data[feature_name];
	for (unsigned int i = 0; i < features.samples.size(); i++) {
		if ((empty || ids.find(features.samples[i].id) != ids.end())) {
			if (data[i] != missing_value) {
				data[i] -= mean;
				if (normalizeSd)
					data[i] /= sd;
			} else if (fillMissing)
					data[i] = 0;
		}
		if (!isfinite(data[i]))
			MTHROW_AND_ERR("FeatureNormalizer sd: %f mean: %f", sd, mean);

	}
	return 0;
}

// Init
//.......................................................................................
int FeatureNormalizer::init(map<string, string>& mapper) {

	init_defaults();

	for (auto entry : mapper) {
		string field = entry.first;

		if (field == "missing_value") missing_value = stof(entry.second);
		else if (field == "normalizeSd") normalizeSd = (stoi(entry.second) != 0);
		else if (field == "fillMissing") fillMissing = (stoi(entry.second) != 0);
		else if (field != "names" && field != "fp_type")
				MLOG("Unknonw parameter \'%s\' for FeatureNormalizer\n", field.c_str());
	}

	return 0;
}

// (De)Serialization
//.......................................................................................
size_t FeatureNormalizer::get_size() {
	return MedSerialize::get_size(processor_type, feature_name, mean, sd, normalizeSd, fillMissing);
}

//extern char signalName_c[MAX_NAME_LEN + 1];

//.......................................................................................
size_t FeatureNormalizer::serialize(unsigned char *blob) {
	return MedSerialize::serialize(blob, processor_type, feature_name, mean, sd, normalizeSd, fillMissing);
}

//.......................................................................................
size_t FeatureNormalizer::deserialize(unsigned char *blob) {
	return MedSerialize::deserialize(blob, processor_type, feature_name, mean, sd, normalizeSd, fillMissing);
}


//=======================================================================================
// FeatureImputer
//=======================================================================================
void FeatureImputer::print()
{
	MLOG("Imputer: Feat: %s nMoments: %d :: ", feature_name.c_str(), moments.size());
	for (auto moment : moments)
		MLOG("%f ", moment);
	MLOG("\n");
}

// Learn
//.......................................................................................
int FeatureImputer::Learn(MedFeatures& features, unordered_set<int>& ids) {

	// Sanity
	if (features.data.find(feature_name) == features.data.end()) {
		MERR("Cannot find signal %s in features data\n", feature_name.c_str());
		return -1;
	}

	for (int i = 0; i < imputerStrata.nStratas(); i++) {
		if (features.data.find(imputerStrata.stratas[i].name) == features.data.end()) {
			MERR("Cannot find signal %s in features data\n", imputerStrata.stratas[i].name.c_str());
			return -1;
		}
	}

	// Get all values
	vector<float> values;
	get_all_values(features, feature_name, ids, values);

	// Get all strata values
	vector<vector<float> > strataValues(imputerStrata.nStratas());
	for (int i = 0; i < imputerStrata.nStratas(); i++)
		get_all_values(features, imputerStrata.stratas[i].name, ids, strataValues[i]);

	// Collect
	imputerStrata.getFactors();

	vector<vector<float> > stratifiedValues(imputerStrata.nValues());
	for (int i = 0; i < values.size(); i++) {
		if (values[i] != missing_value) {
			int index = 0;
			for (int j = 0; j < imputerStrata.nStratas(); j++)
				index += imputerStrata.factors[j] * imputerStrata.stratas[j].getIndex(strataValues[j][i], missing_value);
			stratifiedValues[index].push_back(values[i]);
		}
	}

	// Get moments
	moments.resize(stratifiedValues.size());
	for (unsigned int i = 0; i < stratifiedValues.size(); i++) {
		if (moment_type == IMPUTE_MMNT_MEAN)
			get_mean(stratifiedValues[i], moments[i]);
		else if (moment_type == IMPUTE_MMNT_MEDIAN) {
			if (stratifiedValues[i].size() > 0)
				sort_and_get_median(stratifiedValues[i], moments[i]);
			else
				moments[i] = missing_value;
		}
		else if (moment_type == IMPUTE_MMNT_COMMON)
			get_common(stratifiedValues[i], moments[i]);
		else {
			MERR("Unknown moment type %d for imputing %s\n", moment_type, feature_name.c_str());
		}

	}

//#pragma omp critical
//	print();
	return 0;
}

// Apply
//.......................................................................................
int FeatureImputer::Apply(MedFeatures& features, unordered_set<int>& ids) {

	// Sanity
	if (features.data.find(feature_name) == features.data.end()) {
		MERR("Cannot find signal %s in features data\n", feature_name.c_str());
		return -1;
	}

	for (int i = 0; i < imputerStrata.nStratas(); i++) {
		if (features.data.find(imputerStrata.stratas[i].name) == features.data.end()) {
			MERR("Cannot find signal %s in features data\n", imputerStrata.stratas[i].name.c_str());
			return -1;
		}
	}

	// Attribute
	features.attributes[feature_name].imputed = true;


	// Impute
	imputerStrata.getFactors();
	vector<float>& data = features.data[feature_name];
	vector<vector<float> *> strataData(imputerStrata.nStratas());
	for (int j = 0; j < imputerStrata.nStratas(); j++)
		strataData[j] = &(features.data[imputerStrata.stratas[j].name]);

	for (unsigned int i = 0; i < features.samples.size(); i++) {
		if (data[i] == missing_value) {

			int index = 0;
			for (int j = 0; j < imputerStrata.nStratas(); j++)
				index += imputerStrata.factors[j] * imputerStrata.stratas[j].getIndex((*strataData[j])[i], missing_value);
			data[i] = moments[index];
		}
	}

	return 0;
}

// Init : starta can be a vector, separated by ":"
//.......................................................................................
int FeatureImputer::init(map<string, string>& mapper) {

	init_defaults();
	vector<string> strata;

	for (auto entry : mapper) {
		string field = entry.first;

		if (field == "moment_type") moment_type = (imputeMomentTypes)stoi(entry.second);
		else if (field == "strata") {
			boost::split(strata, entry.second, boost::is_any_of(":"));
			for (string& stratum : strata) addStrata(stratum);
		}
		else if (field != "names" && field != "fp_type")
				MLOG("Unknown parameter \'%s\' for FeatureImputer\n", field.c_str());
	}

	return 0;
}

//.......................................................................................
void FeatureImputer::addStrata(string& init_string) {

	vector<string> fields;
	boost::split(fields, init_string, boost::is_any_of(","));

	if (fields.size() != 4)
		MLOG("Cannot initialize strata from \'%s\'. Ignoring\n", init_string.c_str());
	else
		addStrata(fields[0], stof(fields[3]), stof(fields[1]), stof(fields[2]));

}

// (De)Serialization
//.......................................................................................
size_t FeatureImputer::get_size() {
	return MedSerialize::get_size(processor_type, feature_name, missing_value, imputerStrata, moment_type, moments);
}

//.......................................................................................
size_t FeatureImputer::serialize(unsigned char *blob) {
	return MedSerialize::serialize(blob, processor_type, feature_name, missing_value, imputerStrata, moment_type, moments);
}

//.......................................................................................
size_t FeatureImputer::deserialize(unsigned char *blob) {
	return MedSerialize::deserialize(blob, processor_type, feature_name, missing_value, imputerStrata, moment_type, moments);
}


//.......................................................................................
size_t featureSetStrata::get_size() {
	return MedSerialize::get_size(stratas, factors);
}

//.......................................................................................
size_t featureSetStrata::serialize(unsigned char *blob) {
	return MedSerialize::serialize(blob, stratas, factors);
}

//.......................................................................................
size_t featureSetStrata::deserialize(unsigned char *blob) {
	return MedSerialize::deserialize(blob, stratas, factors);
}

//.......................................................................................
size_t featureStrata::get_size() {
	return MedSerialize::get_size(name, resolution, min, max, nValues);
}

//.......................................................................................
size_t featureStrata::serialize(unsigned char *blob) {
	return MedSerialize::serialize(blob, name, resolution, min, max, nValues);
}


//.......................................................................................
size_t featureStrata::deserialize(unsigned char *blob) {
	return MedSerialize::deserialize(blob, name, resolution, min, max, nValues);
}

//=======================================================================================
// FeatureSelector
//=======================================================================================
// Learn : Add required to feature selected by inheriting classes
//.......................................................................................
int FeatureSelector::Learn(MedFeatures& features, unordered_set<int>& ids) {

	// select, ignoring requirments
	if (_learn(features, ids) < 0)
		return -1;

	// Add required signals
	unordered_set<string> selectedFeatures;
	for (string& feature : selected)
		selectedFeatures.insert(feature);
	
	unordered_set<string> missingRequired;
	for (string feature : required) {
		if (selectedFeatures.find(feature) == selectedFeatures.end())
			missingRequired.insert(feature);
	}

	// Keep maximum numToSelect ...
	if (numToSelect > 0) {
		int nMissing = (int)missingRequired.size();
		int nSelected = (int)selected.size();
		if (nSelected + nMissing > numToSelect)
			selected.resize(numToSelect - nMissing);
	}

	for (string feature : missingRequired)
		selected.push_back(feature);

	return 0;
}

// Apply selection : Ignore set of ids
//.......................................................................................
int FeatureSelector::Apply(MedFeatures& features, unordered_set<int>& ids) {

	unordered_set<string> selectedFeatures;
	for (string& feature : selected)
		selectedFeatures.insert(feature);

	return features.filter(selectedFeatures);
}

// (De)Serialization
//.......................................................................................
size_t FeatureSelector::get_size() {
	return MedSerialize::get_size(processor_type, selected);
}

//.......................................................................................
size_t FeatureSelector::serialize(unsigned char *blob) {
	return MedSerialize::serialize(blob, processor_type, processor_type, selected);
}

//.......................................................................................
size_t FeatureSelector::deserialize(unsigned char *blob) {
	return MedSerialize::deserialize(blob, processor_type, processor_type, selected);
}

//=======================================================================================
// UnivariateFeatureSelector
//=======================================================================================
// Learn 
//.......................................................................................
int UnivariateFeatureSelector::_learn(MedFeatures& features, unordered_set<int>& ids) {

	// Get Stats
	vector<pair<string, float> > stats;
	if (params.method == UNIV_SLCT_PRSN) {
		getAbsPearsonCorrs(features, ids, stats);
	}
	else {
		MERR("Unknown method %d for univariate feature selection\n", params.method);
		return -1;
	}

	// Select
	sort(stats.begin(), stats.end(), [](const pair<string, float> &v1, const pair<string, float> &v2) {return (v1.second > v2.second); });

	if (numToSelect == 0) {
		// Select according to minimum value of stat
		for (auto& rec : stats) {
			if (rec.second < params.minStat)
				break;
			selected.push_back(rec.first);
		}
	}
	else {
		// Select according to number
		int n = (stats.size() > numToSelect) ? numToSelect : (int) stats.size();
		selected.resize(n);
		for (int i = 0; i < n; i++)
			selected[i] = stats[i].first;
	}

	return 0;
}

// Utility : Caluclate pearson correlations
//.......................................................................................
int UnivariateFeatureSelector::getAbsPearsonCorrs(MedFeatures& features, unordered_set<int>& ids, vector<pair<string,float> >& stats) {

	// Get outcome
	vector<float> label;
	get_all_outcomes(features, ids, label);

	stats.resize(features.data.size());
	vector<string> names(stats.size());
	features.get_feature_names(names);
	vector<vector<float>>values(stats.size());

#pragma omp parallel for 
	for (int i = 0; i < names.size(); i++) {
		int n;
		get_all_values(features, names[i], ids, values[i]);
		stats[i].first = names[i];
		stats[i].second = fabs(get_pearson_corr(values[i], label, n, missing_value));
		if (n == 0) stats[i].second = 0.0;
	}
}

// Init
//.......................................................................................
int UnivariateFeatureSelector::init(map<string, string>& mapper) {
	
	init_defaults();

	for (auto entry : mapper) {
		string field = entry.first;

		if (field == "missing_value") missing_value = stof(entry.second);
		else if (field == "numToSelect") numToSelect = stoi(entry.second);
		else if (field == "method") params.method = params.get_method(entry.second);
		else if (field == "minStat") params.minStat = stof(entry.second);
		else if (field == "required") boost::split(required, entry.second, boost::is_any_of(","));
		else if (field != "fp_type")
			MLOG("Unknonw parameter \'%s\' for FeatureSelector\n", field.c_str());
	}

	return 0;

}
//=======================================================================================
// Utilities
//=======================================================================================
//.......................................................................................
void get_all_values(MedFeatures& features, string& signalName, unordered_set<int>& ids, vector<float>& values) {

	if (ids.empty()) {

		int max_sample = 10000;
		int jump = 1;
		//MLOG("taking all values\n");
		int size = (int)features.data[signalName].size();
		if (size > max_sample)
			jump = size/max_sample;
		for (int i=0; i<size; i+=jump)
			values.push_back(features.data[signalName][i]);
		//values = features.data[signalName];

	} else {
		for (unsigned int i = 0; i < features.samples.size(); i++) {
			if (ids.find(features.samples[i].id) != ids.end())
				values.push_back(features.data[signalName][i]);
		}
	}
}

//.......................................................................................
void get_all_outcomes(MedFeatures& features, unordered_set<int>& ids, vector<float>& values) {

	for (unsigned int i = 0; i < features.samples.size(); i++) {
		if (ids.empty() || ids.find(features.samples[i].id) != ids.end())
			values.push_back(features.samples[i].outcome);
	}
}