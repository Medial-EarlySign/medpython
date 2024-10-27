#ifndef __MED__MPBOOTSTRAP__H__
#define __MED__MPBOOTSTRAP__H__

#include "MedPyCommon.h"
#include "MPSamples.h"

class MedBootstrap;

class MPStringBtResultMap {
	bool o_owned = true;
	std::map<std::string, std::map<std::string, float> >* o;
public:
	MPStringBtResultMap();
	MEDPY_IGNORE(MPStringBtResultMap(const MPStringBtResultMap& other));
	MEDPY_IGNORE(MPStringBtResultMap(std::map<std::string, std::map<std::string, float> >* ptr));
	~MPStringBtResultMap();
	int __len__();
	MPStringFloatMapAdaptor __getitem__(const std::string &key);
	void __setitem__(const std::string &key, MPStringFloatMapAdaptor& val);
	std::vector<std::string> keys();
	MEDPY_IGNORE(MPStringBtResultMap& operator=(const MPStringBtResultMap& other));

	void to_file(const std::string &file_path);
};

class MPBootstrap {
public:
	MEDPY_IGNORE(MedBootstrap* o);
	MPBootstrap();
	~MPBootstrap();
	//vector<int> pids;
	//vector<int> split;

	int MEDPY_GET_sample_per_pid();
	void MEDPY_SET_sample_per_pid(int _sample_per_pid);

	int MEDPY_GET_nbootstrap();
	void MEDPY_SET_nbootstrap(int _nbootstrap);

	int MEDPY_GET_ROC_score_min_samples();
	void MEDPY_SET_ROC_score_min_samples(int _ROC_score_min_samples);

	int MEDPY_GET_ROC_score_bins();
	void MEDPY_SET_ROC_score_bins(int _ROC_score_bins);

	float MEDPY_GET_ROC_max_diff_working_point();
	void MEDPY_SET_ROC_max_diff_working_point(float _ROC_max_diff_working_point);

	float MEDPY_GET_ROC_score_resolution();
	void MEDPY_SET_ROC_score_resolution(float _ROC_score_resolutiont);

	bool MEDPY_GET_ROC_use_score_working_points();
	void MEDPY_SET_ROC_use_score_working_points(bool _ROC_use_score_working_points);

	std::vector<float> MEDPY_GET_ROC_working_point_FPR();
	void MEDPY_SET_ROC_working_point_FPR(const std::vector<float> &_ROC_working_point_FPR);

	std::vector<float> MEDPY_GET_ROC_working_point_PR();
	void MEDPY_SET_ROC_working_point_PR(const std::vector<float> &_ROC_working_point_PR);

	std::vector<float> MEDPY_GET_ROC_working_point_SENS();
	void MEDPY_SET_ROC_working_point_SENS(const std::vector<float> &_ROC_working_point_SENS);

	std::vector<float> MEDPY_GET_ROC_working_point_Score();
	void MEDPY_SET_ROC_working_point_Score(const std::vector<float> &_ROC_working_point_Score);

	std::vector<int> MEDPY_GET_ROC_working_point_TOPN();
	void MEDPY_SET_ROC_working_point_TOPN(const std::vector<int> &MEDPY_GET_ROC_working_point_TOPN);

	void parse_cohort_line(const string &line);

	void get_cohort_from_arg(const string &single_cohort);

	void parse_cohort_file(const string &cohorts_path);

	void init_from_str(const string &str);

	MPStringBtResultMap bootstrap_cohort(MPSamples *samples, const string &rep_path, const string &json_model, const string &cohorts_path);

	MPStringFloatMapAdaptor _bootstrap(const std::vector<float> &preds,const std::vector<float> &labels);

	MPStringFloatMapAdaptor _bootstrap_pid(const std::vector<float> &pids, const std::vector<float> &preds, const std::vector<float> &labels);
};

#endif // !__MED__MPSPLIT__H__
