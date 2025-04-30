//
// This file holds calculator implementations 
// Currently for the class RepCalcSimpleSignals

#include <MedProcessTools/MedProcessTools/RepProcess.h>
#include <MedUtils/MedUtils/MedUtils.h>
#include <MedUtils/MedUtils/MedMedical.h>
#include <cmath>

#define LOCAL_SECTION LOG_REPCLEANER
#define LOCAL_LEVEL	LOG_DEF_LEVEL
void *SimpleCalculator::new_polymorphic(string derived_class_name) {
	CONDITIONAL_NEW_CLASS(derived_class_name, RatioCalculator);
	CONDITIONAL_NEW_CLASS(derived_class_name, eGFRCalculator);
	CONDITIONAL_NEW_CLASS(derived_class_name, KfreCalculator);
	CONDITIONAL_NEW_CLASS(derived_class_name, logCalculator);
	CONDITIONAL_NEW_CLASS(derived_class_name, SumCalculator);
	CONDITIONAL_NEW_CLASS(derived_class_name, RangeCalculator);
	CONDITIONAL_NEW_CLASS(derived_class_name, MultiplyCalculator);
	CONDITIONAL_NEW_CLASS(derived_class_name, SetCalculator);
	CONDITIONAL_NEW_CLASS(derived_class_name, ExistsCalculator);
	CONDITIONAL_NEW_CLASS(derived_class_name, EmptyCalculator);
	CONDITIONAL_NEW_CLASS(derived_class_name, ConstantValueCalculator);

	MTHROW_AND_ERR("Warning in SimpleCalculator::new_polymorphic - Unsupported class %s\n", derived_class_name.c_str());
}

//....................................Empty Calculator..................................
void EmptyCalculator::list_output_signals(const vector<string> &input_signals, vector<pair<string, string>> &_virtual_signals, const string &output_type) {

	_virtual_signals.push_back(pair<string, string>("EMPTY_DEFAULT", output_type));
}
void EmptyCalculator::validate_arguments(const vector<string> &input_signals, const vector<string> &output_signals) const {
	if (output_signals.size() != 1)
		MTHROW_AND_ERR("Error EmptyCalculator::validate_arguments - Requires 1 output signals \n");
}

bool EmptyCalculator::do_calc(const vector<float> &vals, float &res) const {
	res = missing_value;
	return false;
}
//....................................Ratio Calculator..................................
int RatioCalculator::init(map<string, string>& mapper) {

	for (auto it = mapper.begin(); it != mapper.end(); ++it)
	{
		//! [RatioCalculator::init]
		if (it->first == "factor")
			factor = stof(it->second);
		else if (it->first == "power_base")
			power_base = stof(it->second);
		else if (it->first == "power_mone")
			power_mone = stof(it->second);
		else if (it->first == "keep_only_in_range")
			keep_only_in_range = stoi(it->second) > 0;
		else
			MTHROW_AND_ERR("Error in RatioCalculator::init - Unsupported argument \"%s\"\n",
				it->first.c_str());
		//! [RatioCalculator::init]
	}

	return 0;
}

void RatioCalculator::validate_arguments(const vector<string> &input_signals, const vector<string> &output_signals) const {
	if (input_signals.size() == 2 && output_signals.size() == 1)
		return;
	MTHROW_AND_ERR("Error RatioCalculator::validate_arguments - Requires 2 input signals and 1 output signal\n");
}

void RatioCalculator::list_output_signals(const vector<string> &input_signals, vector<pair<string, string>> &_virtual_signals, const string &output_type) {
	string o_name = input_signals.front();
	if (power_mone != 1)
		o_name += "^" + medial::print::print_obj(power_mone, "%2.3f");
	for (size_t i = 1; i < input_signals.size(); ++i) {
		o_name += "_over_" + input_signals[i];
		if (power_base != 1)
			o_name += "^" + medial::print::print_obj(power_base, "%2.3f");
	}

	_virtual_signals.push_back(pair<string, string>(o_name, output_type));
}

bool RatioCalculator::do_calc(const vector<float> &vals, float &res) const {
	//assumes vals is 2 values vector = {V1, V2} for calculating := V1 / V2 * factor
	res = missing_value;
	if (vals[1] != 0)
		res = pow(vals[0], power_mone) / pow(vals[1], power_base) * factor;
	return true;
}
//.......................................KFRE Calculator................................
int KfreCalculator::init(map<string, string>& mapper) {

	for (auto it = mapper.begin(); it != mapper.end(); ++it)
	{
		if (it->first == "n_variables")
			n_variables = stoi(it->second);
		else if (it->first == "prediction_years")
			prediction_years = stoi(it->second);
		else if (it->first == "kfre_version")
			kfre_version = stoi(it->second);
		else if (it->first == "region") {
			region = it->second;
			region_id = region2id[region];
		}

		else if (it->first == "discard_range_check")
			discard_range_check = stoi(it->second) > 0;
		else if (it->first == "keep_only_in_range") {
			keep_only_in_range = stoi(it->second) > 0;
			MLOG_V("Set keep_only_in_range=%d\n", keep_only_in_range);
		}
		else
			MTHROW_AND_ERR("Error in KfreCalculator::init - Unsupported argument \"%s\"\n", it->first.c_str());
	}

	return 0;
}

void KfreCalculator::validate_arguments(const vector<string> &input_signals, const vector<string> &output_signals) const {

	/*
	Verify that the number of the input and output parameters is correct
	and that names of the input parameters are correct.
	*/

	// Validate the number of output variables
	if (!(output_signals.size() == 1))
		MTHROW_AND_ERR("Error KfreCalculator::validate_arguments - Requires exactly 1 output signal\n");

	// Validate the number of input variables
	if (!(input_signals.size() == n_variables))
		MTHROW_AND_ERR("Error KfreCalculator::validate_arguments - Requires %d input signals -  got \"%d\"\n", n_variables, (int)input_signals.size());

	// Verify that variable names are as expected
	if (n_variables >= 3)
	{
		// NOTE: "break" statements are missing on purpose
		if (input_signals[0] != "BDATE")
			MTHROW_AND_ERR("Error KfreCalculator::validate_arguments - First signal should be BDATE - got \"%s\"\n", input_signals[0].c_str());
		if (input_signals[1] != "GENDER")
			MTHROW_AND_ERR("Error KfreCalculator::validate_arguments - Second signal should be GENDER- got \"%s\"\n", input_signals[1].c_str());
		if (input_signals[2] != "eGFR_CKD_EPI")
			MTHROW_AND_ERR("Error KfreCalculator::validate_arguments - Third signal should be eGFR_CKD_EPI- got \"%s\"\n", input_signals[2].c_str());
	}
	if (n_variables >= 4)
	{
		if (input_signals[3] != "UrineAlbumin_over_Creatinine")
			MTHROW_AND_ERR("Error KfreCalculator::validate_arguments - Fourth signal should be UrineAlbumin_over_Creatinine- got \"%s\"\n", input_signals[3].c_str());
	}

	if (n_variables >= 8)
	{
		if (input_signals[4] != "Ca")
			MTHROW_AND_ERR("Error KfreCalculator::validate_arguments - Fifth signal should be Ca - got \"%s\"\n", input_signals[4].c_str());
		if (input_signals[5] != "Phosphore")
			MTHROW_AND_ERR("Error KfreCalculator::validate_arguments - Sixth signal should be Phosphore - got \"%s\"\n", input_signals[5].c_str());
		if (input_signals[6] != "Albumin")
			MTHROW_AND_ERR("Error KfreCalculator::validate_arguments - Seventh signal should be Albumin - got \"%s\"\n", input_signals[6].c_str());
		if (input_signals[7] != "Bicarbonate")
			MTHROW_AND_ERR("Error KfreCalculator::validate_arguments - Eighth signal should be Bicarbonate - got \"%s\"\n", input_signals[7].c_str());
	}
}

void KfreCalculator::list_output_signals(const vector<string> &input_signals, vector<pair<string, string>> &_virtual_signals, const string &output_type) {
	//ignores output_type
	if (work_channel == 0)
		_virtual_signals.push_back(pair<string, string>("calc_KFRE", "T(i),V(f)"));
	else
		MTHROW_AND_ERR("ERROR in KfreCalculator::list_output_signals - Unsupported work_channel=%d\n", work_channel);
}

bool KfreCalculator::do_calc(const vector<float> &vals, float &res) const {

	int bdate, year, time;
	float age = 0;
	int   gender = 0;
	float eGFR = 0;
	float UACR = 0;
	float Calcium = 0;
	float Phosphorus = 0;
	float Albumin = 0;
	float Bicarbonate = 0;

	res = missing_value;

	// Verify that the model id is valid
	if (n_variables != 3 && n_variables != 4 && n_variables != 8)
		MTHROW_AND_ERR("Error KfreCalculator::do_calc - n_variables should be either 3,4 or 8\n");

	// Verify that the number of arguments is correct
	// (one more than n_variables, due to the "time"
	if (vals.size() != n_variables + 1)
		MTHROW_AND_ERR("Error KfreCalculator::do_calc -wrong number of arguments. expected %d, got %zu\n",
			n_variables + 1, vals.size());

	// Prepare and validate input arguments

	bool range_check_failed = false;

	// BDATE
	if (vals[0] <= 0)
		return !keep_only_in_range;

	bdate = (int)vals[0]; //BDATE
	year = int(bdate / 10000);

	if (n_variables >= 3)
	{
		// age
		time = (int)vals.back(); //TIME
		age = med_time_converter.get_age(time, MedTime::Date, year);

		if (age < 18 || age > 90)
			range_check_failed = true;

		// gender [unless otherwise stated gender is 1 for males and 2 for females]
		if (vals[1] == 1)
			gender = 1.;
		else
			gender = 0.;

		eGFR = vals[2];
		if (eGFR < 10 || eGFR > 60)
			range_check_failed = true;
	}

	if (n_variables >= 4)
	{
		UACR = vals[3];

		// Convert UACR from mg/mmol to mg/g by multiplying by 8.84
		UACR *= 8.84;

		if (UACR < 10 || UACR > 3000)
			range_check_failed = true;
	}

	if (n_variables >= 8)
	{
		Calcium = vals[4];
		if (Calcium < 7.5 || Calcium > 10.5)
			range_check_failed = true;

		Phosphorus = vals[5];
		if (Phosphorus < 3 || Phosphorus > 6.5)
			range_check_failed = true;

		Albumin = vals[6];
		if (Albumin < 1 || Albumin > 4)
			range_check_failed = true;

		Bicarbonate = vals[7];
		if (Bicarbonate < 15 || Bicarbonate > 28)
			range_check_failed = true;
	}

	// when we perform range check AND it fails
	// there are TWO options: do NOT return a value (keep_only_in_range=TRUE) or return a MISSING value (keep_only_in_range=FALSE[default]) 
	if (!discard_range_check)
		if (range_check_failed)
		{
			//MLOG_V("Observed keep_only_in_range=%d\n", keep_only_in_range);
			res = -1.;
			return !keep_only_in_range;
		}

	bool valid = true;

	if (kfre_version == 0) {

		// Legacy implementation of KFRE.v1 based on 2011 article
		// The code is kept for testing purposes, but can be later removed
		// along with functions get_KFRE_Model_2(), get_KFRE_Model_3, get_KFRE_Model_6

		switch (n_variables)
		{
		case 3:
			res = get_KFRE_Model_2(age, gender, eGFR);
			break;
		case 4:
			res = get_KFRE_Model_3(age, gender, eGFR, UACR);
			break;
		case 8:
			valid = get_KFRE_Model_6(
				res,
				age,
				gender,
				eGFR,
				UACR,
				Calcium,
				Phosphorus,
				Albumin,
				Bicarbonate);
			if (!valid)
			{
				// we get here when intermediate computation overflows
				res = missing_value;
				return false;
			}
			break;
		}
	}
	else if (kfre_version == 1) {

		// Compute KFRE.v1 based on 2011 article

		double baseline;
		vector<double> Coeff;
		vector<double> Xbar;

		FetchCoefficients_v1(
			n_variables,
			baseline,
			Coeff,
			Xbar
		);

		switch (n_variables)
		{
		case 3:
			valid = get_KFRE3(
				res,
				baseline,
				Coeff,
				Xbar,
				age,
				gender,
				eGFR
			);
			if (!valid)
			{
				// we get here when intermediate computation overflows
				res = missing_value;
				return false;
			}
			break;
		case 4:
			valid = get_KFRE4(
				res,
				baseline,
				Coeff,
				Xbar,
				age,
				gender,
				eGFR,
				UACR
			);
			if (!valid)
			{
				// we get here when intermediate computation overflows
				res = missing_value;
				return false;
			}
			break;
		case 8:
			valid = get_KFRE8(
				res,
				baseline,
				Coeff,
				Xbar,
				age,
				gender,
				eGFR,
				UACR,
				Calcium,
				Phosphorus,
				Albumin,
				Bicarbonate
			);
			if (!valid)
			{
				// we get here when intermediate computation overflows
				res = missing_value;
				return false;
			}
			break;
		default:
			MTHROW_AND_ERR("Error KfreCalculator::do_calc - KFRE.v1 only supports 3, 4 and 8 variables\n");
		}
	}
	else if (kfre_version == 2) {
		// Compute KFRE.v2, based on article published in 2016

		//======================================================
		// Prepare coefficients
		//======================================================

		double baseline;
		vector<double> Coeff;
		vector<double> Xbar;

		FetchCoefficients(
			n_variables,
			prediction_years,
			region_id,
			baseline,
			Coeff,
			Xbar
		);

		//======================================================
		// Calculate KFRE score
		//======================================================

		switch (n_variables)
		{
		case 4:
			valid = get_KFRE4(
				res,
				baseline,
				Coeff,
				Xbar,
				age,
				gender,
				eGFR,
				UACR
			);
			if (!valid)
			{
				// we get here when intermediate computation overflows
				res = missing_value;
				return false;
			}
			break;
		case 8:
			valid = get_KFRE8(
				res,
				baseline,
				Coeff,
				Xbar,
				age,
				gender,
				eGFR,
				UACR,
				Calcium,
				Phosphorus,
				Albumin,
				Bicarbonate
			);
			if (!valid)
			{
				// we get here when intermediate computation overflows
				res = missing_value;
				return false;
			}
			break;
		default:
			MTHROW_AND_ERR("Error KfreCalculator::do_calc - KFRE.v2 only implemented for 4 and 8 variables\n");
		}
	}

	return true;
}



//.......................................eGFR Calculator................................
int eGFRCalculator::init(map<string, string>& mapper) {
	for (auto it = mapper.begin(); it != mapper.end(); ++it)
	{
		//! [eGFRCalculator::init]
		if (it->first == "ethnicity")
			ethnicity = stoi(it->second);
		else if (it->first == "mdrd")
			mdrd = stoi(it->second) > 0;
		else if (it->first == "keep_only_in_range")
			keep_only_in_range = stoi(it->second) > 0;
		else
			MTHROW_AND_ERR("Error in eGFRCalculator::init - Unsupported argument \"%s\"\n",
				it->first.c_str());
		//! [eGFRCalculator::init]
	}
	if (mdrd)
		calculator_name = "eGFR_MDRD";
	return 0;
}

void eGFRCalculator::validate_arguments(const vector<string> &input_signals, const vector<string> &output_signals) const {
	//Signals order: "Creatinine", "GENDER", "BDATE"
	if (!(input_signals.size() == 3 && output_signals.size() == 1))
		MTHROW_AND_ERR("Error eGFRCalculator::validate_arguments - Requires 3 input signals and 1 output signal\n");
	if (input_signals[1] != "GENDER")
		MTHROW_AND_ERR("Error eGFRCalculator::validate_arguments - Second signal should be GENDER- got \"%s\"\n", input_signals[1].c_str());
	if (input_signals[2] != "BDATE")
		MTHROW_AND_ERR("Error eGFRCalculator::validate_arguments - Third signal should be BDATE - got \"%s\"\n", input_signals[2].c_str());
}

void eGFRCalculator::list_output_signals(const vector<string> &input_signals, vector<pair<string, string>> &_virtual_signals, const string &output_type) {
	//ignores output_type
	if (work_channel == 0)
		_virtual_signals.push_back(pair<string, string>("calc_eGFR", "T(i),V(f)"));
	else
		MTHROW_AND_ERR("ERROR in eGFRCalculator::list_output_signals - Unsupported work_channel=%d\n", work_channel);
}

bool eGFRCalculator::do_calc(const vector<float> &vals, float &res) const {
	res = missing_value;
	if (vals[0] <= 0)
		return !keep_only_in_range;
	if (vals.size() != 4)
		MTHROW_AND_ERR("Error eGFRCalculator::do_calc -wrong number of arguments. expected 4, got %zu\n",
			vals.size());
	//input: creatinine, gender, byear, time
	int bdate = (int)vals[2]; //BDATE
	int year = int(bdate / 10000);
	int time = (int)vals[3]; //TIME

	float age = med_time_converter.get_age(time, MedTime::Date, year);
	//age, creatinine, gender, ethnicity
	if (!mdrd)
		res = round(get_eGFR_CKD_EPI(age, vals[0], vals[1], ethnicity));
	else
		res = round(get_eGFR_MDRD(age, vals[0], vals[1], ethnicity));
	return true;
}
//.................................LOG CALCULATOR.......................................
int logCalculator::init(map<string, string>& mapper) {
	for (auto it = mapper.begin(); it != mapper.end(); ++it)
	{
		//! [logCalculator::init]
		if (it->first == "keep_only_in_range")
			keep_only_in_range = stoi(it->second) > 0;
		else
			MTHROW_AND_ERR("Error in logCalculator::init - Unsupported argument \"%s\"\n",
				it->first.c_str());
		//! [logCalculator::init]
	}
	return 0;
}
void logCalculator::validate_arguments(const vector<string> &input_signals, const vector<string> &output_signals) const {
	if (!(input_signals.size() == 1 && output_signals.size() == 1))
		MTHROW_AND_ERR("Error logCalculator::validate_arguments - Requires 1 input signals and 1 output signal\n");
}
void logCalculator::list_output_signals(const vector<string> &input_signals, vector<pair<string, string>> &_virtual_signals, const string &output_type) {
	_virtual_signals.push_back(pair<string, string>("log_" + input_signals[0], output_type));
}
bool logCalculator::do_calc(const vector<float> &vals, float &res) const {
	res = missing_value;
	if (vals[0] > 0)
		res = log(vals[0]);
	else
		return !keep_only_in_range;
	return true; //always return
}
//...................................SUM CALCULATOR.....................................
int SumCalculator::init(map<string, string>& mapper) {
	for (auto it = mapper.begin(); it != mapper.end(); ++it)
	{
		//! [SumCalculator::init]
		if (it->first == "factors") {
			vector<string> tokens;
			boost::split(tokens, it->second, boost::is_any_of(",;"));
			factors.resize(tokens.size());
			for (size_t i = 0; i < tokens.size(); ++i)
				factors[i] = stof(tokens[i]);
		}
		else if (it->first == "b0")
			b0 = stof(it->second);
		else if (it->first == "keep_only_in_range")
			keep_only_in_range = stoi(it->second) > 0;
		else
			MTHROW_AND_ERR("Error in SumCalculator::init - Unsupported argument \"%s\"\n",
				it->first.c_str());
		//! [SumCalculator::init]
	}
	return 0;
}

void SumCalculator::validate_arguments(const vector<string> &input_signals, const vector<string> &output_signals) const {
	if (!(input_signals.size() >= 1 && output_signals.size() == 1))
		MTHROW_AND_ERR("Error SumCalculator::validate_arguments - Reqiures at least 1 input signals and 1 output signal\n");
}

void SumCalculator::list_output_signals(const vector<string> &input_signals, vector<pair<string, string>> &_virtual_signals, const string &output_type) {
	if (factors.size() < input_signals.size())
		factors.resize(input_signals.size(), 1);

	char buff[1000];
	if (factors[0] != 1)
		snprintf(buff, sizeof(buff), "%2.3f_X_%s", factors[0], input_signals[0].c_str());
	else
		snprintf(buff, sizeof(buff), "%s", input_signals[0].c_str());
	string o_name = string(buff);

	for (size_t i = 1; i < input_signals.size(); ++i) {
		if (factors[i] != 1)
			snprintf(buff, sizeof(buff), "%s_plus_%2.3f_X_%s",
				o_name.c_str(), factors[i], input_signals[i].c_str());
		else
			snprintf(buff, sizeof(buff), "%s_plus_%s", o_name.c_str(), input_signals[i].c_str());
		o_name = string(buff);
	}
	_virtual_signals.push_back(pair<string, string>(o_name, output_type));
}

bool SumCalculator::do_calc(const vector<float> &vals, float &res) const {
	//no missing values
	res = b0;
	for (size_t i = 0; i < vals.size(); ++i)
		res += factors[i] * vals[i];
	return !keep_only_in_range || res > 0;
}
//...................................RANGE CALCULATOR...................................
int RangeCalculator::init(map<string, string>& mapper) {
	for (auto it = mapper.begin(); it != mapper.end(); ++it)
	{
		//! [RangeCalculator::init]
		if (it->first == "min_range")
			min_range = stof(it->second);
		else if (it->first == "max_range")
			max_range = stof(it->second);
		else if (it->first == "in_range_val")
			in_range_val = stof(it->second);
		else if (it->first == "out_range_val")
			out_range_val = stof(it->second);
		else
			MTHROW_AND_ERR("Error in SumCalculator::init - Unsupported argument \"%s\"\n",
				it->first.c_str());
		//! [RangeCalculator::init]
	}
	if (min_range == MED_MAT_MISSING_VALUE || max_range == MED_MAT_MISSING_VALUE)
		MTHROW_AND_ERR("ERROR: RangeCalculator::init - must provide min_range,max-range arguments\n");
	return 0;
}

void RangeCalculator::validate_arguments(const vector<string> &input_signals, const vector<string> &output_signals) const {
	if (!(input_signals.size() == 1 && output_signals.size() == 1))
		MTHROW_AND_ERR("Error RangeCalculator::validate_arguments - Requires 1 input signals and 1 output signal\n");
}

void RangeCalculator::list_output_signals(const vector<string> &input_signals, vector<pair<string, string>> &_virtual_signals, const string &output_type) {
	char buff[500];
	snprintf(buff, sizeof(buff), "%s_in_%2.3f_%2.3f", input_signals[0].c_str(), min_range, max_range);
	string o_name = string(buff);

	_virtual_signals.push_back(pair<string, string>(o_name, output_type));
}

bool RangeCalculator::do_calc(const vector<float> &vals, float &res) const {
	res = missing_value;
	float v = vals[0];
	if (v >= min_range && v <= max_range)
		res = in_range_val;
	else
		res = out_range_val;
	return true;
}
//.................................MULTIPLY CALCULATOR...................................
int MultiplyCalculator::init(map<string, string>& mapper) {
	for (auto it = mapper.begin(); it != mapper.end(); ++it)
	{
		//! [MultiplyCalculator::init]
		if (it->first == "b0")
			b0 = stof(it->second);
		else if (it->first == "powers") {
			vector<string> tokens;
			boost::split(tokens, it->second, boost::is_any_of(",;"));
			powers.resize(tokens.size());
			for (size_t i = 0; i < tokens.size(); ++i)
				powers[i] = stof(tokens[i]);
		}
		else
			MTHROW_AND_ERR("Error in MultiplyCalculator::init - Unsupported argument \"%s\"\n",
				it->first.c_str());
		//! [MultiplyCalculator::init]
	}

	return 0;
}

void MultiplyCalculator::validate_arguments(const vector<string> &input_signals, const vector<string> &output_signals) const {
	if (!(input_signals.size() >= 1 && output_signals.size() == 1))
		MTHROW_AND_ERR("Error MultiplyCalculator::validate_arguments - Requires at least 1 input signals and 1 output signal\n");
}

void MultiplyCalculator::list_output_signals(const vector<string> &input_signals, vector<pair<string, string>> &_virtual_signals, const string &output_type) {
	if (powers.size() < input_signals.size())
		powers.resize(input_signals.size(), 1);

	char buff[1000];
	if (powers[0] != 1)
		snprintf(buff, sizeof(buff), "%s_^_%2.3f", input_signals[0].c_str(), powers[0]);
	else
		snprintf(buff, sizeof(buff), "%s", input_signals[0].c_str());
	string o_name = string(buff);

	for (size_t i = 1; i < input_signals.size(); ++i) {
		if (powers[i] != 1)
			snprintf(buff, sizeof(buff), "%s_*_%s_^_%2.3f",
				o_name.c_str(), input_signals[i].c_str(), powers[i]);
		else
			snprintf(buff, sizeof(buff), "%s_*_%s", o_name.c_str(), input_signals[i].c_str());
		o_name = string(buff);
	}

	_virtual_signals.push_back(pair<string, string>(o_name, output_type));

}

bool MultiplyCalculator::do_calc(const vector<float> &vals, float &res) const {
	double res_d = b0;
	for (size_t i = 0; i < vals.size(); ++i)
		res_d *= pow(vals[i], powers[i]);
	res = res_d;
	return true;
}
//.............................SET CALCULATOR.........................................
int SetCalculator::init(map<string, string>& mapper) {
	for (auto it = mapper.begin(); it != mapper.end(); ++it)
	{
		//! [SetCalculator::init]
		if (it->first == "sets")
			boost::split(sets, it->second, boost::is_any_of(","));
		else if (it->first == "sets_file")
			medial::io::read_codes_file(it->second, sets);
		else if (it->first == "in_range_val")
			in_range_val = stof(it->second);
		else if (it->first == "out_range_val")
			out_range_val = stof(it->second);
		else if (it->first == "keep_only_in_range")
			keep_only_in_range = stoi(it->second) > 0;
		else if (it->first == "regex_on_sets")
			regex_on_sets = (bool)stoi(it->second) > 0;

		else
			MTHROW_AND_ERR("Error in SumCalculator::init - Unsupported argument \"%s\"\n",
				it->first.c_str());
		//! [SetCalculator::init]
	}
	if (sets.empty())
		MTHROW_AND_ERR("ERROR: SetCalculator::init - must provide min_range,max-range arguments\n");

	return 0;
}

void SetCalculator::validate_arguments(const vector<string> &input_signals, const vector<string> &output_signals) const {
	if (!(input_signals.size() == 1 && output_signals.size() == 1))
		MTHROW_AND_ERR("Error SetCalculator::validate_arguments - Requires 1 input signals and 1 output signal\n");
}

void SetCalculator::list_output_signals(const vector<string> &input_signals, vector<pair<string, string>> &_virtual_signals, const string &output_type) {
	char buff[500];

	snprintf(buff, sizeof(buff), "%s_in_set", input_signals[0].c_str());
	string o_name = string(buff);

	_virtual_signals.push_back(pair<string, string>(o_name, output_type));
	input_signal = input_signals.front();
}

void SetCalculator::init_tables(MedDictionarySections& dict, MedSignals& sigs, const vector<string> &input_signals) {
	int section_id = dict.section_id(input_signals.front());
	if (regex_on_sets)
	{
		unordered_set<string> aggregated_values;
		for (auto& s : sets)
		{
			vector<string> curr_set;
			dict.dicts[section_id].get_regex_names(".*" + s + ".*", curr_set);
			aggregated_values.insert(curr_set.begin(), curr_set.end());
		}
		sets.clear();
		sets.insert(sets.begin(), aggregated_values.begin(), aggregated_values.end());
	}
	dict.prep_sets_lookup_table(section_id, sets, Flags);
}

bool SetCalculator::do_calc(const vector<float> &vals, float &res) const {
	res = missing_value;
	float v = vals[0];
	if (Flags[v])
		res = in_range_val;
	else {
		if (keep_only_in_range)
			return false;
		else
			res = out_range_val;
	}
	return true;
}

void SetCalculator::get_required_signal_categories(unordered_map<string, vector<string>> &signal_categories_in_use) const {
	signal_categories_in_use[input_signal] = sets;
}
//.............................EXISTS CALCULATOR.........................................
int ExistsCalculator::init(map<string, string>& mapper) {
	for (auto it = mapper.begin(); it != mapper.end(); ++it)
	{
		//! [ExistsCalculator::init]
		if (it->first == "in_range_val")
			in_range_val = stof(it->second);
		else if (it->first == "out_range_val")
			out_range_val = stof(it->second);
		else
			MTHROW_AND_ERR("Error in ExistsCalculator::init - Unsupported argument \"%s\"\n",
				it->first.c_str());
		//! [ExistsCalculator::init]
	}
	return 0;
}
void ExistsCalculator::validate_arguments(const vector<string> &input_signals, const vector<string> &output_signals) const {
	if (!(input_signals.size() == 1 && output_signals.size() == 1))
		MTHROW_AND_ERR("Error ExistsCalculator::validate_arguments - Requires 1 input signals and 1 output signal\n");
}
void ExistsCalculator::list_output_signals(const vector<string> &input_signals, vector<pair<string, string>> &_virtual_signals, const string &output_type) {
	_virtual_signals.push_back(pair<string, string>("exists_" + input_signals[0], output_type));
}
bool ExistsCalculator::do_calc(const vector<float> &vals, float &res) const {
	res = out_range_val;
	if (vals.back() > 0)
		res = in_range_val;
	else
		return !keep_only_in_range;
	return true; //always return
}
//.......................................................................................

int ConstantValueCalculator::init(map<string, string>& mapper) {
	for (auto it = mapper.begin(); it != mapper.end(); ++it)
	{
		//! [ConstantValueCalculator::init]
		if (it->first == "is_numeric")
			is_numeric = stoi(it->second) > 0;
		else if (it->first == "value")
			value = it->second;
		else if (it->first == "additional_dict_vals")
			boost::split(additional_dict_vals, it->second, boost::is_any_of(","));
		else
			MTHROW_AND_ERR("Error in ConstantValueCalculator::init - Unsupported argument \"%s\"\n",
				it->first.c_str());
		//! [ConstantValueCalculator::init]
	}
	if (value.empty())
		MTHROW_AND_ERR("Error in ConstantValueCalculator::init - value must be given\n");
	if (is_numeric)
		numeric_val = stof(value);
	return 0;
}
void ConstantValueCalculator::validate_arguments(const vector<string> &input_signals, const vector<string> &output_signals) const {
	if (output_signals.size() != 1)
		MTHROW_AND_ERR("Error ConstantValueCalculator::validate_arguments - Requires 1 output signals \n");
}
void ConstantValueCalculator::list_output_signals(const vector<string> &input_signals, vector<pair<string, string>> &_virtual_signals, const string &output_type) {
	_virtual_signals.push_back(pair<string, string>("DEFAULT_CONSTANT_VALUE", output_type));
}
bool ConstantValueCalculator::do_calc(const vector<float> &vals, float &res) const {
	res = numeric_val;
	//We will need to register this value "1" as categorical if non numerci
	return true; //always return
}

void ConstantValueCalculator::fit_for_repository(MedPidRepository& rep, vector<pair<string, string>> &_virtual_signals) {
	if (is_numeric)
		return;
	//Start new virtual signal with categories:
	for (size_t i = 0; i < _virtual_signals.size(); ++i)
	{
		string vsig_name = _virtual_signals[i].first;
		if (rep.sigs.sid(vsig_name) < 0) {
			MLOG_D("ConstantValueCalculator:: Adding signal %s", vsig_name.c_str());
			rep.sigs.insert_virtual_signal(vsig_name, _virtual_signals[i].second);
			//Store as "categorical"
			int vsig_id = rep.sigs.Name2Sid[vsig_name];
			MLOG_D(" - %d\n", vsig_id);
			rep.sigs.Sid2Info[vsig_id].is_categorical_per_val_channel[0] = 1;
			if (rep.dict.SectionName2Id.find(output_signal_names[0]) == rep.dict.SectionName2Id.end()) {
				MLOG_D("ConstantValueCalculator:: Adding section %s\n", output_signal_names[0].c_str());
				rep.dict.add_section(output_signal_names[0]);
				//sections_names

			}
			int add_section = rep.dict.section_id(vsig_name);
			rep.dict.connect_to_section(output_signal_names[0], add_section);
		}


		int vsig_id = rep.sigs.Name2Sid[vsig_name];
		int add_section = rep.dict.section_id(vsig_name);
		rep.dict.dicts[add_section].Name2Id[vsig_name] = vsig_id;
		rep.dict.dicts[0].Name2Id[vsig_name] = vsig_id;
		rep.dict.dicts[add_section].Id2Name[vsig_id] = vsig_name;
		rep.dict.dicts[add_section].Id2Names[vsig_id] = { vsig_name };
		rep.sigs.Sid2Info[vsig_id].time_unit = rep.sigs.my_repo->time_unit;
		//rep.dict.SectionName2Id[vsig_name] = 0;
		MLOG_D("updated dict %d : %d\n", add_section, rep.dict.dicts[add_section].id(vsig_name));
	}


	int section_id = rep.dict.section_id(output_signal_names[0]);
	//Add categorical value "1" as described:
	MLOG_D("Adding value %s for signal %s in section %d\n", value.c_str(), output_signal_names[0].c_str(), section_id);
	rep.dict.dicts[section_id].push_new_def(value, (int)1);

	//add new value:
	int max_id = rep.dict.dicts[section_id].Id2Name.rbegin()->first;
	for (size_t i = 0; i < additional_dict_vals.size(); ++i)
	{
		++max_id;
		rep.dict.dicts[section_id].push_new_def(additional_dict_vals[i], max_id);
	}
}

//.......................................................................................

//.......................................................................................
//mode 1 (https://en.wikipedia.org/wiki/Model_for_End-Stage_Liver_Disease): 3.78 ln[serum bilirubin (mg/dL)] + 11.2 ln[INR] + 9.57 ln[serum creatinine (mg/dL)] + 6.43
//mode 2 (http://www.thecalculator.co/health/MELD-Na-Score-Calculator-846.html): MELD_mode1 ? Na ? [0.025  MELD_mode1  (140 ? Na)] + 140
// mode2 := meld1 - na - 0.025F * meld1 * (140.0F - na) + 140.0F;
// please use logCalculator with SumCalculator, for mode2 need also MultiplyCalculator

//.......................................................................................
//BMI - defined as weight / height^2. expects weight and height
//please use - RatioCalculator

//.......................................................................................
//APRI - defined as ast / platelets * 2.5, expects ast, platelets.see http ://www.hepatitisc.uw.edu/page/clinical-calculators/apri
//please use - RatioCalculator

//.......................................................................................
//SIDa - defined as [Na+] - [Cl-] (mode 1) or [Na+] + [K+] - [Cl-] see Amland paper on sepsis
//expects Na, K, and Cl
//please use - SumCalculator

//.......................................................................................
//please use - RatioCalculator

//.......................................................................................
//return nervous SOFA
//x = Glasgow Coma Scale score
//please use - 5 RangeCalculator for each range and SumCalculator to unit all
// 4: x [3 , 5], 3: x [6, 9], 2: [10, 12], 1: [13, 14], 0: [15, 15]

//.......................................................................................
//return liver SOFA
//x = serum bilirubin level mg/dl
//please use - 5 RangeCalculator for each range and SumCalculator to unit all

//.......................................................................................
//return coagulation SOFA
//x = plateletCount 10^9/L
//please use - 5 RangeCalculator for each range and SumCalculator to unit all

//combine two sources of medication: first (med2kg) is normalized by kg, the other (med) is not, and weight (kg) in kg is provided
//of course, the dosages should be otherwise in the same units
//please use - RatioCalculator and SumCalculator in second step

//.......................................................................................
//calculate SIRS score. 2 out of the following 4 - Temp > 38 or < 36; Heart Rate > 90; Repoiratory Rate > 20 per min or PaCO2 < 32 mm Hg; 
//WBC > 12 or WBC < 4 or > 10% bands
//expects Temperature, Heart_Rate, Resp_Rate, Art_PaCO2, WBC, Neutrophils_Bands. 
//mode 1: output is 1 if certainly >=2, 0 if certainly < 2, badValue() if undecided
//mode 2: raw score if all constituents present, badValue() otherwise.
//mode 3: optimistic
//mode 4: pessimistic
//Please use RangeCalculator and SumCalculator or write new SIRSCalculator that allow missing values
float calc_hosp_SIRS(float temp, float hr, float resp, float paco2, float wbc, float bands, int mode) {
	int maxSIRS = 0, minSIRS = 0;
	// Check Temperature
	if (temp == MED_MAT_MISSING_VALUE) maxSIRS++; else if (temp < 36 || temp > 38) { maxSIRS++; minSIRS++; }

	// Check Heart Rate
	if (hr == MED_MAT_MISSING_VALUE) maxSIRS++; else if (hr > 90) { maxSIRS++; minSIRS++; }

	// Check Respiratory Rate and PaCO2
	bool c1 = (resp != MED_MAT_MISSING_VALUE && resp > 20);
	bool c2 = (paco2 != MED_MAT_MISSING_VALUE && paco2 < 32);

	if (resp == MED_MAT_MISSING_VALUE || paco2 == MED_MAT_MISSING_VALUE || c1 || c2)
		maxSIRS++;

	if (c1 || c2)
		minSIRS++;

	// Check WBC and bands
	c1 = (wbc != MED_MAT_MISSING_VALUE && (wbc > 12 || wbc < 4));
	c2 = (bands != MED_MAT_MISSING_VALUE && bands > 10);

	if (wbc == MED_MAT_MISSING_VALUE || bands == MED_MAT_MISSING_VALUE || c1 || c2)
		maxSIRS++;

	if (c1 || c2)
		minSIRS++;

	if (mode == 1) {
		if (minSIRS >= 2)
			return 1.0F;
		else if (maxSIRS < 2)
			return 0.0F;
		else
			return MED_MAT_MISSING_VALUE;
	}
	else if (mode == 2) {
		if (maxSIRS == minSIRS)
			return (float)maxSIRS;
		else
			return MED_MAT_MISSING_VALUE;
	}
	else if (mode == 3) {
		return (float)minSIRS;
	}
	else {//no check if actually 4...
		return (float)maxSIRS;
	}
}

//.......................................................................................

//pressure adjusted hr: HR * CVP / MAP
//expects hr, cvp, map
//please use MultiplyCalculator and RatioCalculator

//.......................................................................................
//MODS score. expects paO2, FiO2, platelets, serum bilirubin, hr, cvp, map, gcs, serum creatinine
//http://reference.medscape.com/calculator/mods-score-multiple-organ-dysfunction
//see also Marshall et al. "Multiple organ dysfunction score : a reliable descriptor of a complex clinical outcome"
// write new calculator or use RangeCalculator with SumCalculators
float calc_hosp_MODS(float paO2, float fiO2, float plt, float bili, float hr, float cvp, float map, float gcs, float cre) {
	if (map == 0.0F || fiO2 == 0.0F)
		return MED_MAT_MISSING_VALUE;

	float x;
	int s = 0;

	x = 100.0F * paO2 / fiO2; if (x < 76) s += 4; else if (x < 151) s += 3; else if (x < 226) s += 2; else if (x < 301) s += 1;
	x = plt; if (x < 21) s += 4; else if (x < 51) s += 3; else if (x < 81) s += 2; else if (x < 121) s += 1;
	x = bili; if (x > 14) s += 4; else if (x > 7) s += 3; else if (x > 3.5) s += 2; else if (x > 1.2) s += 1;
	x = hr * cvp / map; if (x > 30) s += 4; else if (x > 20) s += 3; else if (x > 15) s += 2; else if (x > 10) s += 1;
	x = gcs; if (x < 7) s += 4; else if (x < 10) s += 3; else if (x < 13) s += 2; else if (x < 15) s += 1;
	x = cre; if (x > 5.7) s += 4; else if (x > 4.0) s += 3; else if (x > 2.3) s += 2; else if (x > 1.1) s += 1;

	return (float)s;
}

//.......................................................................................
//shockIndex - defined as heart rate divided by systolic blood pressure (mode 1) or mean pressure (mode=2). 
//expects heart rate, systolic bp, diastolic bp
//Use RatioCalculator in mode1 (hr / sysBp), mode2(h1 / (2.0F/3.0F*diaBp + 1.0/3.0*sysBp)): SumCalculator and RatioCalculator

//.......................................................................................
//pulse pressure: systolic bp - diastolic bp
//expects systolic bp, diastolic bp
//use SumCalculator

//.......................................................................................
//feature for various types of eGFR. mode 1: MDRD, mode 2: CKD EPI, mode 3: KeGFR (not implemented)
//useEthnicity - a boolean flag
//expects: creatinine, age, gender (1-male/2-female/missing), isAfricanAmerican (0/1/missing), 
//important: KeGFR looks 72 hours back. 
//implement into eGFR
float calc_hosp_eGFR(float cr, float age, float gender, float isAfricanAmerican, bool useEthnicity, int mode) {
	bool isAfAm = (isAfricanAmerican == 1.0);
	bool isFemale = (gender == 2.0);

	//double MDRD = log(175.0) - 1.154*log((double)cr) - 0.203*log(age) + (isFemale ? log(0.742) : 0.0) + ((isAfAm && useEthnicity) ? log(1.212) : 0.0);
	//kdigo says 186, not 175
	double MDRD = log(186.0) - 1.154*log((double)cr) - 0.203*log(age) + (isFemale ? -0.298406 : 0.0) + ((isAfAm && useEthnicity) ? 0.1922719 : 0.0);
	MDRD = exp(MDRD);

	if (mode == 1) {
		return (float)MDRD;
	}
	else if (mode == 2) {
		double a = (isFemale ? -0.329 : -0.411);
		double k = (isFemale ? 0.7 : 0.9);

		double ckdepi = log(141.0) + a * log(min((double)cr / k, 1.0)) - 1.209 * log(max((double)cr / k, 1.0)) + age * log(0.993) + (isFemale ? log(1.018) : 0.0) +
			((isAfAm && useEthnicity) ? log(1.159) : 0.0);

		ckdepi = exp(ckdepi);
		return (float)ckdepi;
	}

	return -1;
}

//.......................................................................................
//x - mmHg PaO2 / FiO2
//y - is mechanically ventilated, 1 if true, 0 otherwise			
//optimistic - if true, treats missing data as optimistically as possible, otherwise treats it pessimistically.
// use RangeCalculator, write new Calculator
float calc_hosp_SOFA_respiratory(float x, float y, bool optimistic) {
	int s = 0;
	if (x < 100.0 && y > 0.0) s = 4; else if (x < 200.0 && y > 0.0) s = 3; else if (x < 300.0) s = 2; else if (x < 400.0) s = 1;
	return (float)s;
}

//.......................................................................................
//return renal SOFA
//x = serum creatinine level mg/dl
//y = urine output, ml/min
//use range_calculator, write new Calculator
float calc_hosp_SOFA_renal(float x, float y, bool optimistic) {
	int s = 0;
	if (x > 5.0 || y < 200.0) s = 4; else if (x >= 3.5 || y < 500.0) s = 3; else if (x >= 2.0) s = 2; else if (x >= 1.2) s = 1;
	return (float)s;
}

//.......................................................................................
//return cardiovascular SOFA
//x - mean arterial pressure (mmHg)
//y - dopamine, vasopressor g/kg/min
//z - epinephrine, vasopressor g/kg/min
//u - norepinephrine, vasopressor g/kg/min
//v - dobutamine, vasopressor g/kg/min
//write new calculator
float calc_hosp_SOFA_cardio(float x, float y, float z, float u, float v, bool optimistic) {
	int s = 0;
	if (y > 15.0 || z > 0.1 || u > 0.1) s = 4;
	else if (y > 5.0 || z > 0.0 || u > 0.0) s = 3;
	else if (y > 0.0 || v > 0.0) s = 2;
	else if (x < 70.0) s = 1;

	return (float)s;
}

//.......................................................................................
//SOFA score: mode 1 is sum of constituents, mode 2 is max.
//it is assumed that all constituents were calculated similarly - all optimistic or all pessimistic
//as a result, there should be no missing values
float calc_hosp_SOFA(float snerve, float sliver, float scoag, float sresp, float srenal, float scardio, int mode) {
	return (mode == 1 ? (snerve + sliver + scoag + sresp + srenal + scardio) : max(max(max(max(max(snerve, sliver), scoag), sresp), srenal), scardio));
}


//.......................................................................................
//return qSOFA
//gcs - glasgow coma score
//sBp - systolic bp
//resp - respiration
//mode 1: if at least 2 of 3 conditions are fulfilled, return 1, if not, return 0, otherwise (undecided) return badValue
//mode 2 - return the raw number of conditions met, bad value if not clear
//write new calculator
float calc_hosp_qSOFA(float gcs, float sBp, float resp, int mode) {
	int sMin = 0, sMax = 0;
	if (gcs == MED_MAT_MISSING_VALUE)
		sMax++;
	else if (gcs <= 13) {
		sMin++; sMax++;
	}

	if (sBp == MED_MAT_MISSING_VALUE)
		sMax++;
	else if (sBp <= 100) {
		sMin++; sMax++;
	}

	if (resp == MED_MAT_MISSING_VALUE)
		sMax++;
	else if (resp >= 22) {
		sMin++; sMax++;
	}

	if (mode == 1) {
		if (sMax < 2)
			return (float)0;
		else if (sMin >= 2)
			return (float)1;
		else
			return MED_MAT_MISSING_VALUE;
	}
	else { //no check if actually 2...
		if (sMax == sMin)
			return (float)sMin;
		else
			return MED_MAT_MISSING_VALUE;
	}
}

//.......................................................................................
//.......................................................................................