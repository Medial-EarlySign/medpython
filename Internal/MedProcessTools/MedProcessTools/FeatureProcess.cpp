#define _CRT_SECURE_NO_WARNINGS

#define LOCAL_SECTION LOG_FEATCLEANER
#define LOCAL_LEVEL	LOG_DEF_LEVEL

#include "FeatureProcess.h"
#include "DoCalcFeatProcessor.h"
#include "PredictorImputer.h"
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
	else if (processor_name == "iterative_imputer")
		return FTR_PROCESS_ITERATIVE_IMPUTER;
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
	else if (processor_name == "tags_selector")
		return FTR_PROCESSOR_TAGS_SELECTOR;
	else if (processor_name == "importance_selector")
		return FTR_PROCESSOR_IMPORTANCE_SELECTOR;
	else if (processor_name == "iterative_selector")
		return FTR_PROCESSOR_ITERATIVE_SELECTOR;
	else if (processor_name == "pca")
		return FTR_PROCESS_ENCODER_PCA;
	else if (processor_name == "one_hot")
		return FTR_PROCESS_ONE_HOT;
	else if (processor_name == "get_prob")
		return FTR_PROCESS_GET_PROB;
	else if (processor_name == "predictor_imputer")
		return FTR_PROCESS_PREDICTOR_IMPUTER;
	else
		MTHROW_AND_ERR("feature_processor_name_to_type got unknown processor_name [%s]\n", processor_name.c_str());
}

// Initialization
//.......................................................................................
FeatureProcessor* FeatureProcessor::make_processor(string processor_name) {

	return make_processor(feature_processor_name_to_type(processor_name));
}

//.......................................................................................
FeatureProcessor * FeatureProcessor::make_processor(string processor_name, string init_string) {

	FeatureProcessorTypes type = feature_processor_name_to_type(processor_name);
	return make_processor(type, init_string);
}

//.......................................................................................
void *FeatureProcessor::new_polymorphic(string dname)
{
	CONDITIONAL_NEW_CLASS(dname, MultiFeatureProcessor);
	CONDITIONAL_NEW_CLASS(dname, FeatureBasicOutlierCleaner);
	CONDITIONAL_NEW_CLASS(dname, FeatureNormalizer);
	CONDITIONAL_NEW_CLASS(dname, FeatureImputer);
	CONDITIONAL_NEW_CLASS(dname, FeatureIterativeImputer);
	CONDITIONAL_NEW_CLASS(dname, DoCalcFeatProcessor);
	CONDITIONAL_NEW_CLASS(dname, UnivariateFeatureSelector);
	CONDITIONAL_NEW_CLASS(dname, MRMRFeatureSelector);
	CONDITIONAL_NEW_CLASS(dname, LassoSelector);
	CONDITIONAL_NEW_CLASS(dname, DgnrtFeatureRemvoer);
	CONDITIONAL_NEW_CLASS(dname, FeaturePCA);
	CONDITIONAL_NEW_CLASS(dname, TagFeatureSelector);
	CONDITIONAL_NEW_CLASS(dname, ImportanceFeatureSelector);
	CONDITIONAL_NEW_CLASS(dname, IterativeFeatureSelector);
	CONDITIONAL_NEW_CLASS(dname, OneHotFeatProcessor);
	CONDITIONAL_NEW_CLASS(dname, GetProbFeatProcessor);
	CONDITIONAL_NEW_CLASS(dname, PredictorImputer);
	MTHROW_AND_ERR("Warning in FeatureProcessor::new_polymorphic - Unsupported class %s\n", dname.c_str());
	return NULL;
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
	else if (processor_type == FTR_PROCESS_ITERATIVE_IMPUTER)
		return new FeatureIterativeImputer;
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
	else if (processor_type == FTR_PROCESS_ENCODER_PCA)
		return new FeaturePCA;
	else if (processor_type == FTR_PROCESSOR_TAGS_SELECTOR)
		return new TagFeatureSelector;
	else if (processor_type == FTR_PROCESSOR_IMPORTANCE_SELECTOR)
		return new ImportanceFeatureSelector;
	else if (processor_type == FTR_PROCESSOR_ITERATIVE_SELECTOR)
		return new IterativeFeatureSelector;
	else if (processor_type == FTR_PROCESS_ONE_HOT)
		return new OneHotFeatProcessor;
	else if (processor_type == FTR_PROCESS_GET_PROB)
		return new GetProbFeatProcessor;
	else if (processor_type == FTR_PROCESS_PREDICTOR_IMPUTER)
		return new PredictorImputer;
	else
		MTHROW_AND_ERR("make_processor got unknown processor type [%d]\n", processor_type);

}

//.......................................................................................
FeatureProcessor * FeatureProcessor::make_processor(FeatureProcessorTypes processor_type, string init_string) {

	FeatureProcessor *newProcessor = make_processor(processor_type);
	if (newProcessor->init_from_string(init_string) < 0)
		MTHROW_AND_ERR("Cannot init FeatureProcessor of type %d with init string \'%s\'\n", processor_type, init_string.c_str());
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
	return _apply(features, temp);
}

//.......................................................................................
int FeatureProcessor::apply(MedFeatures& features, unordered_set<int>& ids) {
	return _apply(features, ids);
}

//.......................................................................................
int FeatureProcessor::apply(MedFeatures& features, unordered_set<string>& req_features) {

	// All Ids - mark as an empty set
	unordered_set<int> temp;
	return _conditional_apply(features, temp, req_features);
}

//.......................................................................................
int FeatureProcessor::apply(MedFeatures& features, unordered_set<int>& ids, unordered_set<string>& req_features) {
	return _conditional_apply(features, ids, req_features);
}

//.......................................................................................
int FeatureProcessor::_conditional_apply(MedFeatures& features, unordered_set<int>& ids, unordered_set<string>& req_features) {

	if (are_features_affected(req_features))
		return _apply(features, ids);
	return 0;
}

//.......................................................................................
string FeatureProcessor::resolve_feature_name(MedFeatures& features, string substr) {

	// Exact name ?
	if (features.data.find(substr) != features.data.end())
		return substr;
	else {
		vector<string> names;
		for (auto& me : features.data)
			names.push_back(me.first);
		return names[find_in_feature_names(names, substr)];
	}
}

// (De)Serialize
//.......................................................................................
size_t FeatureProcessor::get_processor_size() {
	return sizeof(processor_type) + get_size();
}

//.......................................................................................
size_t FeatureProcessor::processor_serialize(unsigned char *blob) {

	size_t ptr = 0;
	if (processor_type == FTR_PROCESS_LAST)
		MTHROW_AND_ERR("programmer error: trying to serialize a feature_processor with an undefined processor_type, must define it in init_defaults and call from ctor\n");
	memcpy(blob + ptr, &processor_type, sizeof(FeatureProcessorTypes)); ptr += sizeof(FeatureProcessorTypes);
	ptr += serialize(blob + ptr);

	return ptr;
}


//.......................................................................................
void FeatureProcessor::dprint(const string &pref, int fp_flag)
{
	if (fp_flag > 0) {
		MLOG("%s :: FP type %d(%s) : feature_name %s \n", pref.c_str(), processor_type, my_class_name().c_str(), feature_name.c_str());
	}
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
		string tp_name = "";
		if (!processors.empty())
			tp_name = processors.back()->my_class_name();
		MLOG("MultiFeautreProcessor - using duplicate to create %zu processors of type %d(%s)\n",
			features_to_process.size(), members_type, tp_name.c_str());
	}

	int RC = 0;

	// Allow nested parallelism if one processor
	if (processors.size() == 1) {
		use_parallel_learn = false;
		use_parallel_apply = false;
	}

	/*if (!use_parallel_learn && !processors.empty())
		MLOG("no threads for processor %s\n", processors.front()->my_class_name().c_str());*/

#pragma omp parallel for schedule(dynamic) if (use_parallel_learn && processors.size()>1)
	for (int j = 0; j < processors.size(); j++) {
		int rc = processors[j]->Learn(features, ids);
#pragma omp critical
		if (rc < 0) RC = -1;
	}

	return RC;
}

//.......................................................................................
int MultiFeatureProcessor::_apply(MedFeatures& features, unordered_set<int>& ids) {
	int RC = 0;
	if (processors.size() == 1) {
		use_parallel_learn = false;
		use_parallel_apply = false;
	}
#pragma omp parallel for schedule(dynamic) if (use_parallel_apply && processors.size() > 1)
	for (int j = 0; j < processors.size(); j++) {
		int rc = processors[j]->_apply(features, ids);
#pragma omp critical
		if (rc < 0) RC = -1;
	}

	return RC;
}

//.......................................................................................
int MultiFeatureProcessor::_conditional_apply(MedFeatures& features, unordered_set<int>& ids, unordered_set<string>& req_features) {

	int RC = 0;
#pragma omp parallel for schedule(dynamic) if (use_parallel_apply && processors.size() > 1)
	for (int j = 0; j < processors.size(); j++) {
		int rc = processors[j]->_conditional_apply(features, ids, req_features);
#pragma omp critical
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

// Clear
//.......................................................................................
void MultiFeatureProcessor::clear()
{
	for (auto pfp : processors) {
		if (pfp != NULL) {
			delete pfp;
			pfp = NULL;
		}
	}
	processors.clear();
}

/// check if a set of features is affected by the current processor
//.......................................................................................
bool MultiFeatureProcessor::are_features_affected(unordered_set<string>& out_req_features) {

	// Empty set == all features
	if (out_req_features.empty())
		return true;

	for (auto& processor : processors) {
		if (processor->are_features_affected(out_req_features))
			return true;
	}

	return false;
}

/// update sets of required as input according to set required as output to processor
//.......................................................................................
void MultiFeatureProcessor::update_req_features_vec(unordered_set<string>& out_req_features, unordered_set<string>& in_req_features) {

	in_req_features.clear();
	for (auto& processor : processors) {
		unordered_set<string> _in_req_features;
		processor->update_req_features_vec(out_req_features, _in_req_features);
		for (string ftr : _in_req_features)
			in_req_features.insert(ftr);
	}
}

// Init 
//.......................................................................................
int MultiFeatureProcessor::init(map<string, string>& mapper) {

	for (auto entry : mapper) {
		string field = entry.first;
		//! [MultiFeatureProcessor::init]
		if (field == "tag") tag = entry.second;
		if (field == "use_parallel_learn") use_parallel_learn = med_stoi(entry.second) > 0;
		if (field == "use_parallel_apply") use_parallel_apply = med_stoi(entry.second) > 0;
		//! [MultiFeatureProcessor::init]
	}

	return 0;
}

//.......................................................................................
void MultiFeatureProcessor::dprint(const string &pref, int fp_flag)
{
	if (fp_flag > 0) {
		MLOG("%s :: FP MULTI type %d : name %s \n", pref.c_str(), processor_type, feature_name.c_str());
		int ind = 0;
		for (auto& proc : processors) {
			proc->dprint("\t" + pref + "-in-MULTI[" + to_string(ind) + "]", fp_flag);
			++ind;
		}
	}
}

//=======================================================================================
// FeatureBasicOutlierCleaner
//=======================================================================================
// Init from map
//.......................................................................................
int FeatureBasicOutlierCleaner::init(map<string, string>& mapper)
{
	init_defaults();

	for (auto entry : mapper) {
		string field = entry.first;
		//! [FeatureBasicOutlierCleaner::init]
		if (field == "name") feature_name = entry.second;
		//! [FeatureBasicOutlierCleaner::init]
	}

	return MedValueCleaner::init(mapper);
}
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
	int rc = get_iterative_min_max(values);
	if (num_samples_after_cleaning == 0)
		MWARN("FeatureBasicOutlierCleaner::iterativeLearn feature [%s] tried learning cleaning params from an empty vector\n", resolved_feature_name.c_str());
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
int FeatureBasicOutlierCleaner::_apply(MedFeatures& features, unordered_set<int>& ids) {

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

	int n;
	medial::stats::get_mean_and_std(values, missing_value, n, mean, sd);

	// Handle constant vector
	if (sd == 0 && values.size()) {
		MWARN("Got constant (%f) vector in feature %s....\n", feature_name.c_str());
		sd = 1.0;
	}
	else  if (sd == 1)
		MLOG("got sd=1.0 in feature %s....\n", feature_name.c_str());

	if (sd == 0)
		MTHROW_AND_ERR("FeatureNormalizer learn sd: %f mean: %f size: %d", sd, mean, (int)values.size());

	//MLOG("FeatureNormalizer::Learn() done for feature %s , mean %f sd %f size %d\n", feature_name.c_str(), mean, sd, (int)values.size());

	return 0;
}

// Apply
//.......................................................................................
int FeatureNormalizer::_apply(MedFeatures& features, unordered_set<int>& ids) {

	// Resolve
	resolved_feature_name = resolve_feature_name(features, feature_name);

	// Attribute
#pragma omp critical
	{
		features.attributes[resolved_feature_name].normalized = true;
		if (fillMissing)
			features.attributes[resolved_feature_name].imputed = true;
		features.attributes[resolved_feature_name].denorm_mean = mean;
		features.attributes[resolved_feature_name].denorm_sdv = sd;
	}

	// treat resolution
	float multiplier = 1;
	if (resolution > 0)
		multiplier = (float)pow(10, resolution);

	// Clean
	bool empty = ids.empty();
	vector<float>& data = features.data[resolved_feature_name];
	for (unsigned int i = 0; i < features.samples.size(); i++) {
		if ((empty || ids.find(features.samples[i].id) != ids.end())) {
			if (data[i] != missing_value) {
				data[i] -= mean;
				if (normalizeSd)
					data[i] /= sd;
				if (resolution > 0) {
					data[i] = roundf(data[i] * multiplier) / multiplier;
					if (resolution_only)
						data[i] = data[i] * sd + mean;
				}
			}
			else if (fillMissing)
				data[i] = 0;
		}
		if (!isfinite(data[i]))
			MTHROW_AND_ERR("FeatureNormalizer sd: %f mean: %f", sd, mean);
	}

	//MLOG("FeatureNormalizer::Apply() done for feature %s , mean %f sd %f size %d flags: normalized %d imputed %d\n", 
	//	feature_name.c_str(), mean, sd, (int)data.size(), (int)features.attributes[resolved_feature_name].normalized, (int)features.attributes[resolved_feature_name].imputed);

	return 0;
}

// Init
//.......................................................................................
int FeatureNormalizer::init(map<string, string>& mapper) {

	init_defaults();

	for (auto entry : mapper) {
		string field = entry.first;
		//! [FeatureNormalizer::init]
		if (field == "missing_value") missing_value = stof(entry.second);
		else if (field == "normalizeSd") normalizeSd = (med_stoi(entry.second) != 0);
		else if (field == "resolution_only") resolution_only = (med_stoi(entry.second) != 0);
		else if (field == "fillMissing") fillMissing = (med_stoi(entry.second) != 0);
		else if (field == "max_samples") max_samples = med_stoi(entry.second);
		else if (field == "resolution") resolution = med_stoi(entry.second);
		else if (field != "names" && field != "fp_type" && field != "tag")
			MLOG("Unknonw parameter \'%s\' for FeatureNormalizer\n", field.c_str());
		//! [FeatureNormalizer::init]
	}

	return 0;
}

//=======================================================================================
// FeatureImputer
//=======================================================================================
void FeatureImputer::print()
{
	if (moment_type == IMPUTE_MMNT_SAMPLE) {
		MLOG("Imputer: Feat: %s nHistograms: %d :: ", feature_name.c_str(), histograms.size());
		for (unsigned int i = 0; i < histograms.size(); i++) {
			for (auto& pair : histograms[i])
				MLOG("%d %f L %f", i, pair.first, pair.second);
		}
		MLOG("\n");
	}
	else {
		MLOG("Imputer: Feat: %s nMoments: %d :: ", feature_name.c_str(), moments.size());
		for (auto moment : moments)
			MLOG("%f ", moment);
		MLOG("\n");
	}
}

// Convert partial feature names to full names (including FTR_...)
//.......................................................................................
void FeatureImputer::check_stratas_name(MedFeatures& features, map <string, string> &strata_name_conversion)
{
	for (int i = 0; i < imputerStrata.nStratas(); i++) {
		if (strata_name_conversion.find(imputerStrata.stratas[i].name) != strata_name_conversion.end())
			// already mapped
			continue;
		strata_name_conversion[imputerStrata.stratas[i].name] = resolve_feature_name(features, imputerStrata.stratas[i].name);
	}
}

// Learn
//.......................................................................................
int FeatureImputer::Learn(MedFeatures& features, unordered_set<int>& ids) {
	// Resolve
	resolved_feature_name = resolve_feature_name(features, feature_name);
	default_moment = missing_value; //initialize
	map <string, string> strata_name_conversion;
	check_stratas_name(features, strata_name_conversion);

	// Get all values
	vector<float> values;
	get_all_values(features, resolved_feature_name, ids, values, max_samples);
	// Get all strata values
	vector<vector<float> > strataValues(imputerStrata.nStratas());
	for (int i = 0; i < imputerStrata.nStratas(); i++) {
		string resolved_strata_name = resolve_feature_name(features, imputerStrata.stratas[i].name);
		get_all_values(features, resolved_strata_name, ids, strataValues[i], max_samples);
	}

	// Collect
	imputerStrata.getFactors();

	vector<vector<float> > stratifiedValues(imputerStrata.nValues());
	vector<float> all_existing_values;
	for (int j = 0; j < values.size(); j++) {
		if (values[j] != missing_value) {
			all_existing_values.push_back(values[j]);
			int index = 0;
			for (int i = 0; i < imputerStrata.nStratas(); i++)
				index += imputerStrata.factors[i] * imputerStrata.stratas[i].getIndex(strataValues[i][j], missing_value);
			stratifiedValues[index].push_back(values[j]);
		}
	}
	//for (int j = 0; j < stratifiedValues.size(); j++)
		//MLOG("collected %d %d\n", j, stratifiedValues[j].size());

	// Get moments
	if (moment_type == IMPUTE_MMNT_SAMPLE)
		histograms.resize(stratifiedValues.size());
	else
		moments.resize(stratifiedValues.size());

	strata_sizes.resize(stratifiedValues.size());
	int too_small_stratas = 0;
	for (unsigned int i = 0; i < stratifiedValues.size(); i++) {

		strata_sizes[i] = (int)stratifiedValues[i].size();
		if (strata_sizes[i] < min_samples) { // Not enough values to make valid imputation
			too_small_stratas++;
			if (moment_type == IMPUTE_MMNT_SAMPLE)
				histograms[i].push_back({ missing_value,(float)1.0 });
			else
				moments[i] = missing_value;
		}
		else if (moment_type == IMPUTE_MMNT_MEAN)
			moments[i] = medial::stats::mean_without_cleaning(stratifiedValues[i]);
		else if (moment_type == IMPUTE_MMNT_MEDIAN) {
			if (stratifiedValues[i].size() > 0)
				moments[i] = medial::stats::median_without_cleaning(stratifiedValues[i]);
			else
				moments[i] = missing_value;
		}
		else if (moment_type == IMPUTE_MMNT_COMMON)
			moments[i] = medial::stats::most_common_without_cleaning(stratifiedValues[i]);
		else if (moment_type == IMPUTE_MMNT_SAMPLE) {
			medial::stats::get_histogram_without_cleaning(stratifiedValues[i], histograms[i]);
		}
		else MTHROW_AND_ERR("Unknown moment type %d for imputing %s\n", moment_type, feature_name.c_str());
	}
	if (all_existing_values.size() < min_samples) {
		MLOG("WARNING: FeatureImputer::Learn found only %d < %d samples over all for [%s], will not learn to impute it\n",
			all_existing_values.size(), min_samples, feature_name.c_str());
		if (moment_type == IMPUTE_MMNT_SAMPLE)
			default_histogram.push_back({ missing_value,(float)1.0 });
		else
			default_moment = missing_value;
	}
	else {
		if (too_small_stratas > 0) {
			if (!leave_missing_for_small_stratas)
			{
				MLOG("WARNING: FeatureImputer::Learn found less than %d samples for %d/%d stratas for [%s], will learn to impute them using all values\n",
					min_samples, too_small_stratas, stratifiedValues.size(), feature_name.c_str());
				if (moment_type == IMPUTE_MMNT_MEAN)
					default_moment = medial::stats::mean_without_cleaning(all_existing_values);
				else if (moment_type == IMPUTE_MMNT_MEDIAN)
					default_moment = medial::stats::median_without_cleaning(all_existing_values);
				else if (moment_type == IMPUTE_MMNT_COMMON)
					default_moment = medial::stats::most_common_without_cleaning(all_existing_values);
				else if (moment_type == IMPUTE_MMNT_SAMPLE)
					medial::stats::get_histogram_without_cleaning(all_existing_values, default_histogram);
			}
			else {
				// leave_missing_for_small_stratas = true
				MLOG("WARNING: FeatureImputer::Learn found less than %d samples for %d/%d stratas for [%s], will NOT impute them using all values\n",
					min_samples, too_small_stratas, stratifiedValues.size(), feature_name.c_str());
				if (moment_type == IMPUTE_MMNT_SAMPLE)
					default_histogram.push_back({ missing_value,(float)1.0 });
				else
					default_moment = missing_value;
			}
		}
	}
	//for (int j = 0; j < moments.size(); j++)
		//MLOG("moment %d = [%f]\n", j, moments[j]);

//#pragma omp critical
//	print();
	if (verbose_learn)
		print();
	return 0;
}

// Apply
//.......................................................................................
int FeatureImputer::_apply(MedFeatures& features, unordered_set<int>& ids) {

	// Resolve
	resolved_feature_name = resolve_feature_name(features, feature_name);

	map <string, string> strata_name_conversion;
	check_stratas_name(features, strata_name_conversion);
	// Attribute
#pragma omp critical
	features.attributes[resolved_feature_name].imputed = true;

	// Impute
	imputerStrata.getFactors();
	vector<float>& data = features.data[resolved_feature_name];
	vector<vector<float> *> strataData(imputerStrata.nStratas());
	for (int j = 0; j < imputerStrata.nStratas(); j++) {
		string resolved_strata_name = resolve_feature_name(features, imputerStrata.stratas[j].name);
		strataData[j] = &(features.data[resolved_strata_name]);
	}

	int missing_cnt = 0;
	for (unsigned int i = 0; i < features.samples.size(); i++) {
		if (data[i] == missing_value) {
			int index = 0;
			for (int j = 0; j < imputerStrata.nStratas(); j++)
				index += imputerStrata.factors[j] * imputerStrata.stratas[j].getIndex((*strataData[j])[i], missing_value);
			if (moment_type == IMPUTE_MMNT_SAMPLE) {
				if (strata_sizes[index] < min_samples)
					data[i] = medial::stats::sample_from_histogram(default_histogram);
				else
					data[i] = medial::stats::sample_from_histogram(histograms[index]);
			}
			else {
				if (strata_sizes[index] < min_samples)
					data[i] = default_moment;
				else
					data[i] = moments[index];
			}
			if (!isfinite(data[i]))
				MTHROW_AND_ERR("[%s] imputed illegal value for row %d moment_type %d index %d strata_sizes[index] %d %f\n",
					resolved_feature_name.c_str(), i, moment_type, index, strata_sizes[index], default_moment);
			++missing_cnt;
		}
	}

	if (verbose && missing_cnt > 0) {
		MLOG_D("FeatureImputer::%s:: with %d imputations out of %zu(%2.2f%%)\n",
			resolved_feature_name.c_str(), missing_cnt, data.size(), 100.0 * missing_cnt / double(data.size()));
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
		//! [FeatureImputer::init]
		if (field == "name") feature_name = entry.second;
		else if (field == "min_samples") min_samples = med_stoi(entry.second);
		else if (field == "moment_type") moment_type = getMomentType(entry.second);
		else if (field == "max_samples") max_samples = med_stoi(entry.second);
		else if (field == "strata") {
			boost::split(strata, entry.second, boost::is_any_of(":"));
			for (string& stratum : strata) addStrata(stratum);
		}
		else if (field == "verbose")
			verbose = stoi(entry.second) > 0;
		else if (field == "verbose_learn")
			verbose_learn = stoi(entry.second) > 0;
		else if (field == "leave_missing_for_small_stratas") leave_missing_for_small_stratas = (bool)med_stoi(entry.second);
		else if (field != "names" && field != "fp_type" && field != "tag")
			MLOG("Unknown parameter \'%s\' for FeatureImputer\n", field.c_str());
		//! [FeatureImputer::init]
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
	else if (entry == "3" || entry == "sample")
		return IMPUTE_MMNT_SAMPLE;
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

/// update sets of required as input according to set required as output to processor
//.......................................................................................
void FeatureImputer::update_req_features_vec(unordered_set<string>& out_req_features, unordered_set<string>& in_req_features) {

	in_req_features = out_req_features;

	// Check if imputer is actually applied
	if (out_req_features.find(feature_name) != out_req_features.end()) {
		// Add signals required for imputation
		for (int i = 0; i < imputerStrata.nStratas(); i++)
			in_req_features.insert(imputerStrata.stratas[i].name);
	}
}

void FeatureImputer::dprint(const string &pref, int fp_flag) {
	if (fp_flag > 0) {
		MLOG("%s :: FP type %d(%s) : feature_name %s :: default_moment %f \n", pref.c_str(),
			processor_type, my_class_name().c_str(), feature_name.c_str(), default_moment);
	}
}

//=======================================================================================
// OneHotFeatProcessor
//=======================================================================================
int OneHotFeatProcessor::init(map<string, string>& mapper) {

	for (auto entry : mapper) {
		string field = entry.first;
		//! [OneHotFeatProcessor::init]
		if (field == "name") feature_name = entry.second;
		else if (field == "prefix") index_feature_prefix = entry.second;
		else if (field == "remove_origin") rem_origin = (med_stoi(entry.second) != 0);
		else if (field == "add_other") add_other = (med_stoi(entry.second) != 0);
		else if (field == "allow_other") allow_other = (med_stoi(entry.second) != 0);
		else if (field == "remove_last") remove_last = (med_stoi(entry.second) != 0);
		else if (field == "max_values") max_values = med_stoi(entry.second);
		else if (field != "names" && field != "fp_type" && field != "tag")
			MLOG("Unknown parameter \'%s\' for OneHotFeatProcessor\n", field.c_str());
		//! [OneHotFeatProcessor::init]
	}

	// Set output names
	if (index_feature_prefix == "")
		index_feature_prefix = feature_name;

	return 0;
}

int OneHotFeatProcessor::Learn(MedFeatures& features, unordered_set<int>& ids) {

	// Resolve
	resolved_feature_name = resolve_feature_name(features, feature_name);
	string out_prefix = resolved_feature_name;
	boost::replace_first(out_prefix, feature_name, index_feature_prefix);

	// Get all values
	vector<float> values;
	get_all_values(features, resolved_feature_name, ids, values, 0);

	// Build value2feature
	unordered_set<float> all_values(values.begin(), values.end());
	if (all_values.size() > max_values)
		MTHROW_AND_ERR("Found %zd different values for %s. More than allowed %d\n", all_values.size(), feature_name.c_str(), max_values);

	;
	for (float value : all_values) {
#pragma omp critical
		// value2feature[value] = "FTR_" + int_to_string_digits(++MedFeatures::global_serial_id_cnt, 6) + "." + index_feature_prefix + "_";
		value2feature[value] = out_prefix + "_";

		if (features.attributes[resolved_feature_name].value2Name.empty())
			value2feature[value] += to_string(value);
		else if (value == features.medf_missing_value)
			value2feature[value] += "MISSING_VALUE";
		else {
			if (features.attributes[resolved_feature_name].value2Name.find(value) == features.attributes[resolved_feature_name].value2Name.end())
				MTHROW_AND_ERR("Cannot find value %f for in feature %s value2Name\n", value, resolved_feature_name.c_str());
			value2feature[value] += features.attributes[resolved_feature_name].value2Name[value];
		}
	}

	other_feature_name = "FTR_" + int_to_string_digits(++MedFeatures::global_serial_id_cnt, 6) + "." + index_feature_prefix + "_other";

	// Remove last one
	if (remove_last && !value2feature.empty())
		removed_feature_name = value2feature.rbegin()->second;

	return 0;
}

int OneHotFeatProcessor::_apply(MedFeatures& features, unordered_set<int>& ids) {


	// Prepare new Features
	int samples_size = (int)features.samples.size();
	for (auto& rec : value2feature) {
		string feature_name = rec.second;
		if (feature_name != removed_feature_name)
#pragma omp critical
		{
			features.data[feature_name].clear();
			features.data[feature_name].resize(samples_size, 0.0);
			// Attributes
			features.attributes[feature_name].normalized = false;
			features.attributes[feature_name].imputed = true;
		}
	}


	if (add_other) {
#pragma omp critical
		{
			features.data[other_feature_name].clear();
			features.data[other_feature_name].resize(samples_size, 0.0);
			// Attributes
			features.attributes[other_feature_name].normalized = false;
			features.attributes[other_feature_name].imputed = true;
		}
	}

	// Fill it up
	for (int i = 0; i < samples_size; i++) {
		if (ids.empty() || ids.find(features.samples[i].id) != ids.end()) {
			float value = features.data[resolved_feature_name][i];
			if (value2feature.find(value) != value2feature.end()) {
				if (value2feature[value] != removed_feature_name)
					features.data[value2feature[value]][i] = 1.0;
			}
			else {
				if (add_other)
					features.data[other_feature_name][i] = 1.0;
				else if (!allow_other)
					MTHROW_AND_ERR("Unknown value %f for feature %s\n", value, feature_name.c_str());
			}
		}
	}

	// Remove original, if required
#pragma omp critical
	if (rem_origin) {
		features.data.erase(resolved_feature_name);
		features.attributes.erase(resolved_feature_name);
	}

	return 0;
}

/// check if a set of features is affected by the current processor
//.......................................................................................
bool OneHotFeatProcessor::are_features_affected(unordered_set<string>& out_req_features) {

	// If empty = all features are required
	if (out_req_features.empty())
		return true;

	// Otherwise - check in generated features
	for (auto& rec : value2feature) {
		string feature_name = rec.second;
		if (out_req_features.find(feature_name) != out_req_features.end())
			return true;
	}

	if (add_other &&out_req_features.find(other_feature_name) != out_req_features.end())
		return true;

	return false;
}

/// update sets of required as input according to set required as output to processor
//.......................................................................................
void OneHotFeatProcessor::update_req_features_vec(unordered_set<string>& out_req_features, unordered_set<string>& in_req_features) {

	// If empty, keep as is
	if (out_req_features.empty())
		in_req_features.clear();
	else {
		in_req_features = out_req_features;
		// If active, than add original 
		if (are_features_affected(out_req_features))
			in_req_features.insert(resolved_feature_name);
	}
}




//=======================================================================================
// GetProbFeatProcessor
//=======================================================================================
//.......................................................................................
int GetProbFeatProcessor::Learn(MedFeatures& features, unordered_set<int>& ids) {

	// Sanity
	if (!target_labels.empty() && all_labels)
		MTHROW_AND_ERR("GetProbFeatProcessor Error: both all_labels and target_labels given\n");

	// Resolve
	resolved_feature_name = resolve_feature_name(features, feature_name);

	// Fill target labels
	if (all_labels) {
		unordered_set<float> all_labels_set;
		for (auto& sample : features.samples)
			all_labels_set.insert(sample.outcome);

		int idx = 0;
		for (float label : all_labels_set)
			target_labels[label] = idx++;


	}

	// Get all values
	vector<float> values;
	get_all_values(features, resolved_feature_name, ids, values, (int)features.samples.size());

	// Learn Probs
	int nlabels = target_labels.empty() ? 1 : (int)target_labels.size();
	map<float, int> nums;
	vector<map<float, int>> pos_nums(nlabels);
	int overall_num = 0;
	vector<int> overall_pos_num(nlabels);

	if (target_labels.empty()) { // Binary outcome

		for (unsigned int i = 0; i < values.size(); i++) {
			if (values[i] != missing_value) {
				nums[values[i]] ++;
				overall_num++;

				if (features.samples[i].outcome) {
					pos_nums[0][values[i]] ++;
					overall_pos_num[0]++;
				}
			}
		}
	}
	else { // Multi-categorical
		for (unsigned int i = 0; i < values.size(); i++) {
			if (values[i] != missing_value) {
				nums[values[i]] ++;
				overall_num++;

				float outcome = features.samples[i].outcome;
				if (target_labels.find(outcome) != target_labels.end()) {
					pos_nums[target_labels[outcome]][values[i]] ++;
					overall_pos_num[target_labels[outcome]]++;
				}
			}
		}

		for (auto& rec : target_labels)
			feature_names[rec.first] = resolved_feature_name + "_" + to_string(rec.first);
	}

	if (overall_num == 0)
		MTHROW_AND_ERR("Cannot learn Get-Prob feature processor on an empty vector for %s\n", feature_name.c_str());

	overall_prob.resize(nlabels);
	probs.resize(nlabels);
	for (int i = 0; i < nlabels; i++) {
		overall_prob[i] = (overall_pos_num[i] + 0.0) / overall_num;
		for (auto& rec : nums)
			if (rec.second >= min_obs)
				probs[i][rec.first] = (pos_nums[i][rec.first] + overall_count * overall_prob[i]) / (nums[rec.first] + overall_count);
			else
				probs[i][rec.first] = overall_prob[i];
	}

	return 0;
}

// Apply
//.......................................................................................
int GetProbFeatProcessor::_apply(MedFeatures& features, unordered_set<int>& ids) {

	//cerr << "Apply\n";
	// Resolve
	resolved_feature_name = resolve_feature_name(features, feature_name);

	// Transform
	bool empty = ids.empty();
	vector<float>& data = features.data[resolved_feature_name];

	if (target_labels.empty()) { // Single outcome. inplace	
		for (unsigned int i = 0; i < features.samples.size(); i++) {
			if ((empty || ids.find(features.samples[i].id) != ids.end())) {
				if (data[i] == missing_value || probs[0].find(data[i]) == probs[0].end())
					data[i] = overall_prob[0];
				else
					data[i] = probs[0][data[i]];
			}
		}
	}
	else { // Multiple outcomes. new features

		// Prepare new Features
		int samples_size = (int)features.samples.size();
		for (auto& rec : feature_names) {
			string feature_name = rec.second;
#pragma omp critical
			{
				features.data[feature_name].clear();
				features.data[feature_name].resize(samples_size, 0.0);
				// Attributes
				features.attributes[feature_name].normalized = false;
				features.attributes[feature_name].imputed = true;
			}
		}

		// Fill
		for (unsigned int i = 0; i < features.samples.size(); i++) {
			if ((empty || ids.find(features.samples[i].id) != ids.end())) {
				if (data[i] == missing_value || probs[0].find(data[i]) == probs[0].end()) {
					for (auto& rec : feature_names)
						features.data[rec.second][i] = overall_prob[target_labels[rec.first]];
				}
				else {
					for (auto& rec : feature_names)
						features.data[rec.second][i] = probs[target_labels[rec.first]][data[i]];
				}
			}
		}

		// Remove original, if required
#pragma omp critical
		if (remove_origin) {
			features.data.erase(resolved_feature_name);
			features.attributes.erase(resolved_feature_name);
		}
	}


	return 0;
}

// Init
//.......................................................................................
int GetProbFeatProcessor::init(map<string, string>& mapper) {

	init_defaults();

	for (auto entry : mapper) {
		string field = entry.first;
		//! [GetProbFeatProcessor::init]
		if (field == "name") feature_name = entry.second;
		else if (field == "missing_value") missing_value = stof(entry.second);
		else if (field == "overall_count") overall_count = med_stoi(entry.second);
		else if (field == "min_obs") min_obs = med_stoi(entry.second);
		else if (field == "remove_origin") remove_origin = (med_stoi(entry.second) != 0);
		else if (field == "target_labels") {
			vector<string> labels;
			boost::split(labels, entry.second, boost::is_any_of(","));
			for (int i = 0; i < (int)labels.size(); i++)
				target_labels[stof(labels[i])] = i;
		}
		else if (field == "all_labels") all_labels = (med_stoi(entry.second) != 0);
		else if (field != "names" && field != "fp_type" && field != "tag")
			MLOG("Unknonw parameter \'%s\' for GetProbFeatProcessor\n", field.c_str());
		//! [GetProbFeatProcessor::init]
	}

	return 0;
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
			jump = size / max_sample;

		vector<float>& dataVec = features.data[signalName];
		for (int i = 0; i < size; i += jump)
			values.push_back(dataVec[i]);

	}
	else {
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

		for (int i = 0; i < size; i += jump)
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

	float f = (float)nBins / (float)reqNbins;
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
