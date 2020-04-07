#include "PredictorImputer.h"
#include <algorithm>
#include <random>

using namespace std;

#define LOCAL_SECTION LOG_MED_MODEL
#define LOCAL_LEVEL	LOG_DEF_LEVEL

void PredictorImputer::init_defaults() {
	processor_type = FTR_PROCESS_PREDICTOR_IMPUTER;
	missing_value = MED_MAT_MISSING_VALUE;
	feature_name = "";
	resolved_feature_name = "";
	tag_search = "";
	generator_args = "";
	sampling_args = "";
	verbose_learn = true;
	verbose_apply = false;
	gen_type = GeneratorType::GIBBS;

	random_device rd;
	gen = mt19937(rd());
}

void PredictorImputer::post_deserialization() {
	init_sampler(false);
}

void PredictorImputer::init_sampler(bool with_sampler) {
	switch (gen_type)
	{
	case GIBBS:
		if (with_sampler) {
			_gibbs.init_from_string(generator_args);
			_sampler = unique_ptr<SamplesGenerator<float>>(new GibbsSamplesGenerator<float>(_gibbs, true));
		}
		_gibbs_sample_params.init_from_string(sampling_args);
		_gibbs_sample_params.samples_count = 1;
		sampler_sampling_args = &_gibbs_sample_params;
		break;
	case GAN:
		if (with_sampler) {
			_sampler = unique_ptr<SamplesGenerator<float>>(new MaskedGAN<float>);
			static_cast<MaskedGAN<float> *>(_sampler.get())->read_from_text_file(generator_args);
			static_cast<MaskedGAN<float> *>(_sampler.get())->mg_params.init_from_string(sampling_args);
		}
		break;
	case MISSING:
		if (with_sampler)
			_sampler = unique_ptr<SamplesGenerator<float>>(new MissingsSamplesGenerator<float>(missing_value));
		break;
	case RANDOM_DIST:
		if (with_sampler)
			_sampler = unique_ptr<SamplesGenerator<float>>(new RandomSamplesGenerator<float>(0, 5));
		sampler_sampling_args = &n_masks;
		break;
	case UNIVARIATE_DIST:
		if (with_sampler) {
			_sampler = unique_ptr<SamplesGenerator<float>>(new UnivariateSamplesGenerator<float>);
			if (!generator_args.empty())
				_sampler->init_from_string(generator_args);
		}
		break;
	default:
		MTHROW_AND_ERR("Error in ShapleyExplainer::init_sampler() - Unsupported Type %d\n", gen_type);
	}
}

int PredictorImputer::init(map<string, string>& mapper) {

	for (auto &it : mapper)
	{
		//! [PredictorImputer::init]
		if (it.first == "missing_value")
			missing_value = med_stof(it.second);
		else if (it.first == "generator_args")
			generator_args = it.second;
		else if (it.first == "verbose_learn")
			verbose_learn = med_stoi(it.second) > 0;
		else if (it.first == "sampling_args")
			sampling_args = it.second;
		else if (it.first == "verbose_apply")
			verbose_apply = med_stoi(it.second) > 0;
		else if (it.first == "tag_search" || it.first == "tag")
			tag_search = it.second;
		else if (it.first == "gen_type")
			gen_type = GeneratorType_fromStr(it.second);
		else if (it.first == "fp_type" || it.first == "use_parallel_learn" || it.first == "use_parallel_apply") {}
		else
			MTHROW_AND_ERR("Error in PredictorImputer::init - unsupported argument %s\n", it.first.c_str());
		//! [PredictorImputer::init]
	}

	init_sampler();

	return 0;
}

int PredictorImputer::Learn(MedFeatures& features, unordered_set<int>& ids) {
	//Work on tag_search features:
	vector<string> all_names;
	features.get_feature_names(all_names);
	for (const auto &it : features.data)
	{
		const unordered_set<string> &f_tags = features.tags.at(it.first);
		if (it.first == tag_search || f_tags.find(tag_search) != f_tags.end() || tag_search.empty())
			impute_features.push_back(it.first);
	}
	//add option to learn only on impute_features:
	_sampler->learn(features.data, impute_features, true);

	return 0;
}

void PredictorImputer::load_GIBBS(const GibbsSampler<float> &gibbs, const GibbsSamplingParams &sampling_args) {
	_gibbs = gibbs;
	_gibbs_sample_params = sampling_args;

	sampler_sampling_args = &_gibbs_sample_params;
	_sampler = unique_ptr<SamplesGenerator<float>>(new GibbsSamplesGenerator<float>(_gibbs, true));

	gen_type = GeneratorType::GIBBS;
}

void PredictorImputer::load_GAN(const string &gan_path) {
	_sampler = unique_ptr<SamplesGenerator<float>>(new MaskedGAN<float>);
	static_cast<MaskedGAN<float> *>(_sampler.get())->read_from_text_file(gan_path);
	static_cast<MaskedGAN<float> *>(_sampler.get())->mg_params.init_from_string(sampling_args);

	gen_type = GeneratorType::GAN;
}

void PredictorImputer::load_MISSING() {
	_sampler = unique_ptr<SamplesGenerator<float>>(new MissingsSamplesGenerator<float>(missing_value));
	gen_type = GeneratorType::MISSING;
}

void PredictorImputer::load_sampler(unique_ptr<SamplesGenerator<float>> &&generator) {
	_sampler = move(generator);
}

int PredictorImputer::_apply(MedFeatures& features, unordered_set<int>& ids) {
	if (features.data.empty())
		return 0;
	int nsamples = (int)features.data.begin()->second.size();
	_sampler->prepare(sampler_sampling_args);
	vector<string> all_names;
	features.get_feature_names(all_names);
	unordered_set<string> imput_set(impute_features.begin(), impute_features.end());
	vector<vector<float> *> fast_access(all_names.size());
	for (size_t i = 0; i < all_names.size(); ++i)
		fast_access[i] = &features.data.at(all_names[i]);

	//impute samples:
	if (verbose_apply)
		MLOG("Start PredictorImputer::apply for %d samples\n", nsamples);
	if (_sampler->use_vector_api) {
		vector<bool> mask(features.data.size(), true);
		MedMat<float> res; //the result matrix
		MedMat<float> mat_inp;
		features.get_as_matrix(mat_inp);
		vector<vector<bool>> masks(nsamples);
		for (size_t i = 0; i < nsamples; ++i)
		{
			vector<bool> &b_vec = masks[i]; //imput falses
			b_vec.resize(features.data.size());
			for (size_t k = 0; k < features.data.size(); ++k)
				b_vec[k] = imput_set.find(all_names[k]) == imput_set.end() || mat_inp(i, k) != missing_value;
		}

		_sampler->get_samples(res, 1, sampler_sampling_args, masks, mat_inp, gen);
		//copy to features!
		for (size_t i = 0; i < nsamples; ++i) {
			const vector<bool> &b_vec = masks[i]; //imput falses
			for (size_t k = 0; k < features.data.size(); ++k)
				if (!b_vec[k])
					fast_access[k]->at(i) = res(i, k);
		}
	}
	else {
		random_device rd;
		int N_TH = omp_get_max_threads();
		vector<mt19937> gens(N_TH);
		for (size_t i = 0; i < N_TH; ++i)
			gens[i] = mt19937(rd());
		MedProgress progress("PredictorImputer::apply", nsamples, 30, 10);
#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < nsamples; ++i)
		{
			int th_n = omp_get_thread_num();
			mt19937 &thread_gen = gens[th_n];
			vector<float> x(features.data.size());
			vector<bool> mask(features.data.size(), true);
			map<string, vector<float>> gen_matrix;
			for (size_t k = 0; k < x.size(); ++k)
			{
				x[k] = fast_access[k]->at(i);
				mask[k] = (x[k] != missing_value) ||
					imput_set.find(all_names[k]) == imput_set.end();
			}

			_sampler->get_samples(gen_matrix, sampler_sampling_args, mask, x, thread_gen);
			//copy to features!
//#pragma omp critical
			for (size_t k = 0; k < x.size(); ++k)
				if (!mask[k])
					fast_access[k]->at(i) = gen_matrix.at(all_names[k])[0];
			progress.update();
		}

	}
	return 0;
}