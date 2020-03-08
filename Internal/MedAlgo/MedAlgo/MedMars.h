#pragma once
#include <MedAlgo/MedAlgo/MedAlgo.h>
//======================================================================================
// Mars: C++ version of Mars
//======================================================================================
#define MED_MARS_DEF_MAXTERMS			50
#define MED_MARS_DEF_DEGREE				2
#define MED_MARS_DEF_PENALTY			3
#define MED_MARS_DEF_THRESH				0.001
#define MED_MARS_DEF_MINSPAN			0
#define MED_MARS_DEF_PRUNE				1
#define MED_MARS_DEF_FASTK				20
#define MED_MARS_DEF_FASTBETA			0
#define MED_MARS_DEF_NEWVAR_PENALTY		0
#define MED_MARS_DEF_USEBETACACHE		1
#define MED_MARS_DEF_TRACE				3

struct MedMarsParams {

	// Required
	int MaxTerms;	///< Maximal number of Terms in final model
	int MaxDegree;	///< Model max linear degree : 1 - linear, 2 - square, 3 - cubic, etc...
	double Penalty;
	double Thresh;
	int MinSpan;
	bool Prune;
	int FastK;
	double FastBeta;
	double NewVarPenalty;
	bool UseBetaCache;
	double Trace;	///< debug prints during algorithm run (recommended): 0: no prints , 3: print all (recommended)
};

class MedMars : public MedPredictor {
public:
	// Model 
	int nMaxTerms;
	int nTerms;
	int nPreds;
	bool *BestSet;			///< size: nMaxTerms
	vector<int> Dirs;		///< size: nMaxTerms*nPreds
	vector<double> Cuts;	///< size: nMaxTerms*nPreds
	vector<double> Betas;	///< size: nMaxTerms

							// Model Inner quality measures
	double BestGcv;
	vector<double> bx;
	vector<double> Residuals;

	// Parameters
	MedMarsParams params;

	// Function
	MedMars();
	MedMars(void *params);
	MedMars(MedMarsParams& params);
	int init(void *params);
	/// The parsed fields from init command.
	/// @snippet MedMars.cpp MedMars::init
	virtual int set_params(map<string, string>& mapper);
	void init_defaults();

	int Learn(float *x, float *y, const float *w, int nsamples, int nftrs);
	int Predict(float *x, float *&preds, int nsamples, int nftrs) const;

	//int denormalize_model(float *f_avg, float *f_std, float lavel_avg, float label_std) {return 0;};

	// (De)Desrialize - virtual class methods that do the actuale (De)Serializing. Should be created for each predictor
	ADD_CLASS_NAME(MedMars)
		size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

	// Print
	void print(FILE *fp, const string& prefix, int level = 0) const;

	// Predictions per sample
	int n_preds_per_sample() const;
};

// Initialization of parameters
void init_default_mars_params(MedMarsParams& _params);

MEDSERIALIZE_SUPPORT(MedMars)