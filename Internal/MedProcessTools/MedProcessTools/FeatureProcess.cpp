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
	else if (processor_name == "mrmr" || processor_name == "mrmr_selector")
		return FTR_PROCESSOR_MRMR_SELECTOR;
	else if (processor_name == "lasso")
		return FTR_PROCESSOR_LASSO_SELECTOR;
	else if (processor_name == "remove_deg")
		return FTR_PROCESS_REMOVE_DGNRT_FTRS;
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
	else if (processor_type == FTR_PROCESSOR_MRMR_SELECTOR)
		return new MRMRFeatureSelector;
	else if (processor_type == FTR_PROCESSOR_LASSO_SELECTOR)
		return new LassoSelector;
	else if (processor_type == FTR_PROCESS_REMOVE_DGNRT_FTRS)
		return new DgnrtFeatureRemvoer;
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

//.......................................................................................
string FeatureProcessor::resolve_feature_name(MedFeatures& features, string substr) {

	string real_feature_name = "";

	// Exact name ?
	if (features.data.find(substr) != features.data.end())
		return substr;

	// Or ...
	for (auto candidate : features.attributes)
		if (candidate.first.find(substr) != string::npos) {
			if (real_feature_name != "")
				throw runtime_error(string("source_feature_name [") + substr + "] matches both [" + real_feature_name + "] and [" + candidate.first + "]");
			real_feature_name = candidate.first;
		}
	if (real_feature_name == "") {
		string err = string("source_feature_name [") + substr + "] does not match any feature. Tried matching to these features:\n";
		for (auto candidate : features.attributes)
			err += candidate.first + "\n";
		throw runtime_error(err);
	}

	return real_feature_name;
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
	if (processors.size() == 0 && duplicate) {
		vector<string> features_to_process;
		for (auto& rec : features.data) {
			string name = rec.first;
			if (tag == "" || features.tags[name].find(tag) != features.tags[name].end())
				features_to_process.push_back(name);
		}
		add_processors_set(members_type, features_to_process, init_string);
	}

	int RC = 0;

	// Allow nested parallelism if one processor
	if (processors.size() == 1) omp_set_nested(1);

#pragma omp parallel for schedule(dynamic)
	for (int j = 0; j<processors.size(); j++) {
		int rc = processors[j]->Learn(features, ids);
//#pragma omp critical
		if (rc < 0) RC = -1;
	}

	omp_set_nested(0);
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

//.......................................................................................
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

// Init 
//.......................................................................................
int MultiFeatureProcessor::init(map<string, string>& mapper) {

	for (auto entry : mapper) {
		string field = entry.first;

		if (field == "tag") tag = entry.second;
	}

	return 0;
}

// Serialization
//.......................................................................................
size_t MultiFeatureProcessor::get_size() {

	size_t size = MedSerialize::get_size(members_type, init_string, duplicate, tag);
	size += sizeof(int);

	for (auto& processor : processors)
		size += processor->get_processor_size();

	return size;
}

//.......................................................................................
size_t MultiFeatureProcessor::serialize(unsigned char *blob) {

	size_t ptr = MedSerialize::serialize(blob, members_type, init_string, duplicate, tag);

	int nProcessors = (int)processors.size();
	memcpy(blob + ptr, &nProcessors, sizeof(int)); ptr += sizeof(int);

	for (auto& processor : processors)
		ptr += processor->processor_serialize(blob + ptr);

	return ptr;
}

//.......................................................................................
size_t MultiFeatureProcessor::deserialize(unsigned char *blob) {

	size_t ptr = MedSerialize::deserialize(blob, members_type, init_string, duplicate, tag);

	// number of processors
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

	// Resolve
	resolved_feature_name = resolve_feature_name(features, feature_name);

	// Get all values
	vector<float> values;
	get_all_values(features, resolved_feature_name, ids, values, params.max_samples);

	// Get bounds
	if (values.size() == 0)
		MWARN("EMPTY_VECTOR:: feature [%s] has 0 values\n", resolved_feature_name.c_str());

	int rc = get_iterative_min_max(values);
	return rc;
}

//.......................................................................................
int FeatureBasicOutlierCleaner::quantileLearn(MedFeatures& features, unordered_set<int>& ids) {

	// Resolve
	resolved_feature_name = resolve_feature_name(features, feature_name);

	// Get all values
	vector<float> values;
	get_all_values(features, resolved_feature_name, ids, values, params.max_samples);

	// Get bounds
	return get_quantile_min_max(values);
}

// Clean
//.......................................................................................
int FeatureBasicOutlierCleaner::Apply(MedFeatures& features, unordered_set<int>& ids) {

	// Resolve
	resolved_feature_name = resolve_feature_name(features, feature_name);

	// Clean
	bool empty = ids.empty();
	vector<float>& data = features.data[resolved_feature_name];
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
	return MedSerialize::get_size(processor_type, feature_name, resolved_feature_name, params.doTrim, params.doRemove, trimMax, trimMin, removeMax, removeMin);
}

//.......................................................................................
size_t FeatureBasicOutlierCleaner::serialize(unsigned char *blob) {
	return MedSerialize::serialize(blob, processor_type, feature_name, resolved_feature_name, params.doTrim, params.doRemove, trimMax, trimMin, removeMax, removeMin);
}

//.......................................................................................
size_t FeatureBasicOutlierCleaner::deserialize(unsigned char *blob) {
	return MedSerialize::deserialize(blob, processor_type, feature_name, resolved_feature_name, params.doTrim, params.doRemove, trimMax, trimMin, removeMax, removeMin);
}

//=======================================================================================
// FeatureNormalizer
//=======================================================================================
//.......................................................................................
int FeatureNormalizer::Learn(MedFeatures& features, unordered_set<int>& ids) {

	// Resolve
	resolved_feature_name = resolve_feature_name(features, feature_name);

	// Get all values
	vector<float> values;
	get_all_values(features, resolved_feature_name, ids, values, max_samples);

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

	// Resolve
	resolved_feature_name = resolve_feature_name(features, feature_name);

	// Attribute
	features.attributes[resolved_feature_name].normalized = true;
	if (fillMissing)
		features.attributes[resolved_feature_name].imputed = true;

	// Clean
	bool empty = ids.empty();
	vector<float>& data = features.data[resolved_feature_name];
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
		else if (field == "max_samples") max_samples = stoi(entry.second);
		else if (field != "names" && field != "fp_type" && field != "tag")
				MLOG("Unknonw parameter \'%s\' for FeatureNormalizer\n", field.c_str());
	}

	return 0;
}

// (De)Serialization
//.......................................................................................
size_t FeatureNormalizer::get_size() {
	return MedSerialize::get_size(processor_type, feature_name, resolved_feature_name, mean, sd, normalizeSd, fillMissing);
}

//extern char signalName_c[MAX_NAME_LEN + 1];

//.......................................................................................
size_t FeatureNormalizer::serialize(unsigned char *blob) {
	return MedSerialize::serialize(blob, processor_type, feature_name, resolved_feature_name, mean, sd, normalizeSd, fillMissing);
}

//.......................................................................................
size_t FeatureNormalizer::deserialize(unsigned char *blob) {
	return MedSerialize::deserialize(blob, processor_type, resolved_feature_name, feature_name, mean, sd, normalizeSd, fillMissing);
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

	// Resolve
	resolved_feature_name = resolve_feature_name(features, feature_name);

	for (int i = 0; i < imputerStrata.nStratas(); i++) {
		if (features.data.find(imputerStrata.stratas[i].name) == features.data.end()) {
			MERR("Cannot find signal %s in features data\n", imputerStrata.stratas[i].name.c_str());
			return -1;
		}
	}

	// Get all values
	vector<float> values;
	get_all_values(features, resolved_feature_name, ids, values, max_samples);

	// Get all strata values
	vector<vector<float> > strataValues(imputerStrata.nStratas());
	for (int i = 0; i < imputerStrata.nStratas(); i++)
		get_all_values(features, imputerStrata.stratas[i].name, ids, strataValues[i], max_samples);

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

	// Resolve
	resolved_feature_name = resolve_feature_name(features, feature_name);

	for (int i = 0; i < imputerStrata.nStratas(); i++) {
		if (features.data.find(imputerStrata.stratas[i].name) == features.data.end()) {
			MERR("Cannot find signal %s in features data\n", imputerStrata.stratas[i].name.c_str());
			return -1;
		}
	}

	// Attribute
	features.attributes[resolved_feature_name].imputed = true;

	// Impute
	imputerStrata.getFactors();
	vector<float>& data = features.data[resolved_feature_name];
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

		if (field == "moment_type") moment_type = getMomentType(entry.second); 
		else if (field == "max_samples") max_samples = stoi(entry.second);
		else if (field == "strata") {
			boost::split(strata, entry.second, boost::is_any_of(":"));
			for (string& stratum : strata) addStrata(stratum);
		}
		else if (field != "names" && field != "fp_type" && field != "tag")
				MLOG("Unknown parameter \'%s\' for FeatureImputer\n", field.c_str());
	}

	return 0;
}

//.......................................................................................
imputeMomentTypes FeatureImputer::getMomentType(string& entry) {

	boost::to_lower(entry);
	if (entry == "0" || entry == "mean")
		return IMPUTE_MMNT_MEAN;
	else if (entry == "1" || entry == "median")
		return IMPUTE_MMNT_MEDIAN;
	else if (entry == "2" || entry == "common")
		return IMPUTE_MMNT_COMMON;
	else
		return IMPUTE_MMNT_LAST;
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
	return MedSerialize::get_size(processor_type, feature_name, resolved_feature_name, missing_value, imputerStrata, moment_type, moments);
}

//.......................................................................................
size_t FeatureImputer::serialize(unsigned char *blob) {
	return MedSerialize::serialize(blob, processor_type, feature_name, resolved_feature_name, missing_value, imputerStrata, moment_type, moments);
}

//.......................................................................................
size_t FeatureImputer::deserialize(unsigned char *blob) {
	return MedSerialize::deserialize(blob, processor_type, feature_name, resolved_feature_name, missing_value, imputerStrata, moment_type, moments);
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
// Utilities
//=======================================================================================
//.......................................................................................
void get_all_values(MedFeatures& features, string& signalName, unordered_set<int>& ids, vector<float>& values, int max_sample) {

	values.clear();
	if (ids.empty()) {

		int jump = 1;
		int size = (int)features.data[signalName].size();
		if (max_sample > 0 && max_sample < size)
			jump = size/max_sample;
		
		for (int i=0; i<size; i+=jump)
			values.push_back(features.data[signalName][i]);

	} else {
		for (unsigned int i = 0; i < features.samples.size(); i++) {
			if (ids.find(features.samples[i].id) != ids.end())
				values.push_back(features.data[signalName][i]);
		}
	}
}

//.......................................................................................
void get_all_outcomes(MedFeatures& features, unordered_set<int>& ids, vector<float>& values, int max_sample) {

	values.clear();
	if (ids.empty()) {

		int jump = 1;
		int size = (int)features.samples.size();
		if (max_sample > 0 && max_sample < size)
			jump = size / max_sample;

		for (int i = 0; i<size; i += jump)
			values.push_back(features.samples[i].outcome);
		//values = features.data[signalName];

	}
	else {
		for (unsigned int i = 0; i < features.samples.size(); i++) {
			if (ids.find(features.samples[i].id) != ids.end())
				values.push_back(features.samples[i].outcome);
		}
	}
}

//.......................................................................................
void smearBins(vector<int>& bins, int nBins, int reqNbins) {

	float f = (float)nBins/(float)reqNbins;
	vector<vector<int> > newBins(nBins);
	for (int iBin = 0; iBin < reqNbins; iBin++) {
		int OrigBin = (int)(iBin * f);
		newBins[OrigBin].push_back(iBin);
	}

	for (int i = 0; i < bins.size(); i++) {
		int origBin = bins[i];
		int nNewBins = (int)newBins[origBin].size();
		bins[i] = newBins[origBin][nNewBins*(rand() / ((int)RAND_MAX))];
	}
}
