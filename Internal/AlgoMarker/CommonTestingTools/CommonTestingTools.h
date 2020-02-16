#ifndef __COMMONTESTINGTOOLS_H
#define __COMMONTESTINGTOOLS_H

#include <AlgoMarker/AlgoMarker/AlgoMarker.h>
#include <AlgoMarker/DynAMWrapper/DynAMWrapper.h>
#include <AlgoMarker/CommonTestingTools/CommonTestingTools.h>
#include <Logger/Logger/Logger.h>
#include <MedProcessTools/MedProcessTools/MedModel.h>
#include <MedProcessTools/MedProcessTools/MedSamples.h>
#include <MedIO/MedIO/MedIO.h>
#include <json/json.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string.hpp>
#include <algorithm>
#include <string>
#include <vector>
#include "DataLoader.h"


#ifdef __linux__ 
#include <wordexp.h>
#elif _WIN32
#include "windows.h" 
#endif

namespace CommonTestingTools {

	class DataLoader;

	const map<int, int> code_to_status_tbl = {
		{ 300, 2 },
	{ 301, 2 },
	{ 310, 2 },
	{ 311, 2 },
	{ 320, 1 },
	{ 321, 2 },
	{ 390, 0 },
	{ 391, 1 },
	{ 392, 2 }
	};

	const map<string, string> units_tbl = {
		{ "BMI" , "kg/m^2" },
	{ "Glucose" , "mg/dL" },
	{ "HbA1C" , "%" },
	{ "HDL" , "mg/dL" },
	{ "Triglycerides" , "mg/dL" },
	{ "ALT" , "U/L" },
	{ "RBC" , "10^6/uL" },
	{ "Na" , "mmol/L" },
	{ "Weight" , "Kg" },
	{ "WBC" , "10^3/uL" },
	{ "Basophils#" , "#" },
	{ "Basophils%" , "%" },
	{ "Eosinophils#" , "#" },
	{ "Eosinophils%" , "%" },
	{ "Hematocrit" , "%" },
	{ "Hemoglobin" , "g/dL" },
	{ "Lymphocytes#" , "#" },
	{ "Lymphocytes%" , "%" },
	{ "MCH" , "pg/cell" },
	{ "MCHC-M" , "g/dL" },
	{ "MCV" , "fL" },
	{ "Monocytes#" , "#" },
	{ "Monocytes%" , "%" },
	{ "MPV" , "mic*3" },
	{ "Neutrophils#" , "#" },
	{ "Neutrophils%" , "%" },
	{ "Platelets" , "10^3/uL" },
	{ "RDW" , "%" },
	{ "MSG" , "#" }
	};


	using namespace std;

	string precision_float_to_string(float val);
	json read_json_array_next_chunk(ifstream& infile, bool& in_array);

	//Expand string with embedded Environment variables in it
	string expandEnvVars(const string &str);

	// convert a C++ vector of strings to a char**
	class charpp_adaptor : public vector<string> {
	protected:
		char** charpp_arr;
		char* charpp_buf;
	public:
		void init() {
			charpp_arr = (char**)malloc(1);
			charpp_buf = (char*)malloc(1);
		}
		~charpp_adaptor() {
			free(charpp_arr);
			free(charpp_buf);
		};
		charpp_adaptor() : vector<string>() { init(); };
		charpp_adaptor(int capacity) : vector<string>(capacity) { init(); };
		charpp_adaptor(const charpp_adaptor& other) : vector<string>(other) { init(); };

		char** get_charpp();
	};

	// get a malloc'ed read-write copy of a vector's .data() pointer
	template<typename T>
	class get_volatile_data_adaptor {
	protected:
		T * data_p;
		size_t data_p_size;
		int n_elem;
	public:
		~get_volatile_data_adaptor() {
			free(data_p);
		};

		T* from_vec(const vector<T>& orig) {
			n_elem = (int)orig.size();
			data_p_size = n_elem * sizeof(T);
			if (data_p_size == 0)
				return nullptr;
			data_p = (T*)realloc(data_p, data_p_size);
			memcpy(data_p, orig.data(), data_p_size);
			return data_p;
		};

		get_volatile_data_adaptor(const vector<T>& orig) : n_elem(0), data_p_size(0) { data_p = (T*)malloc(1); from_vec(orig); }

		T* get_volatile_data() {
			if (data_p_size == 0)
				return nullptr;
			return data_p;
		}
	};

	json json_AddData(const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, float* Values, int n_time_channels, int n_val_channels);
	json json_AddDataStr(const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, char** Values, int n_time_channels, int n_val_channels);

	int get_preds_from_algomarker(AlgoMarker *am, vector<MedSample> &res, bool print_msgs, DataLoader& d, bool force_add_data, ofstream& msgs_stream, vector<string> ignore_sig, bool extended_score=false);
	int get_preds_from_algomarker_single(AlgoMarker *am, vector<MedSample> &res, bool print_msgs, DataLoader& d, bool force_add_data, ofstream& msgs_stream, vector<string> ignore_sig, ofstream& json_reqfile_stream, bool extended_score=false);
	void save_sample_vec(vector<MedSample> sample_vec, const string& fname);
}

#endif // __COMMONTESTINGTOOLS_H