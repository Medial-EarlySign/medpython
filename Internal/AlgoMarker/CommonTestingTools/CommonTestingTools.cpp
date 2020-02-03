#include <AlgoMarker/CommonTestingTools/CommonTestingTools.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>

#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL	LOG_DEF_LEVEL


string CommonTestingTools::precision_float_to_string(float val) {
	stringstream ss;
	ss << std::setprecision(10) << val;
	return ss.str();
}

//Expand string with embedded Environment variables in it
string CommonTestingTools::expandEnvVars(const string &str) {
	string ret = "";
#ifdef __linux__ 
	wordexp_t p;
	char** w;
	wordexp(str.c_str(), &p, 0);
	w = p.we_wordv;
	for (size_t i = 0; i < p.we_wordc; i++) ret += w[i];
	wordfree(&p);
#elif _WIN32
	DWORD max_str_len = 4 * 1024;
	auto buf = new char[max_str_len];
	DWORD req_len = ExpandEnvironmentStrings(str.c_str(), buf, max_str_len);
	if (req_len > max_str_len) {
		delete buf;
		buf = new char[req_len];
		req_len = ExpandEnvironmentStrings(str.c_str(), buf, req_len);
	}
	if (req_len > 0)
		ret = buf;
	delete buf;
#endif
	return ret;
}


char** CommonTestingTools::charpp_adaptor::get_charpp() {
	if (this->size() == 0)
		return nullptr;
	size_t charpp_arr_sz = this->size() * sizeof(char*);
	size_t charpp_buf_sz = 0;
	for (auto& str : *this) {
		charpp_buf_sz += str.size() + 1;
	}
	charpp_buf_sz *= sizeof(char);

	charpp_arr = (char**)realloc(charpp_arr, charpp_arr_sz);
	charpp_buf = (char*)realloc(charpp_buf, charpp_buf_sz);

	char** charpp_arr_i = charpp_arr;
	char* charpp_buf_i = charpp_buf;

	for (auto& str : *this)
	{
		*charpp_arr_i = charpp_buf_i;
		charpp_arr_i++;
		for (int i = 0; i < str.size(); ++i) {
			*charpp_buf_i = str[i];
			charpp_buf_i++;
		}
		*charpp_buf_i = '\0';
		charpp_buf_i++;
	}
	return charpp_arr;
}

json CommonTestingTools::read_json_array_next_chunk(ifstream& infile, bool& in_array) {
	char prev_c = '\0';
	char c = '\0';
	bool in_string = false;
	string ret_str = "";
	int block_depth = 0;
	while (infile.get(c)) {
		switch (c) {
		case '"':
			if (!in_string)
				in_string = true;
			else if (prev_c != '\\')
				in_string = false;
			break;
		case '{':
			if (!in_array)
				throw runtime_error("File should be a JSON array containing objects");
			if (!in_string)
				block_depth++;
			break;
		case '}':
			if (!in_array)
				throw runtime_error("Did not expect a '}'");
			if (!in_string)
				block_depth--;
			if (block_depth < 0)
				throw runtime_error("Did not expect a '}'");
			break;
		}
		if (c == '[' && !in_array) {
			in_array = true;
			continue;
		}
		if ((c == ']' || c == ',') && in_array && block_depth == 0)
			break;
		ret_str += c;
		prev_c = c;
	}
	json ret;
	if (std::all_of(ret_str.begin(), ret_str.end(), [](char c) { return c == ' ' || c == '\n' || c == '\t' || c == '\r'; }))
		return ret;
	try {
		ret = json::parse(ret_str);
	}
	catch (...) {
		MERR("Error parsing chunk: \n'%s'\n", ret_str.c_str());
	}
	return ret;
}

json CommonTestingTools::json_AddData(const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, float* Values, int n_time_channels, int n_val_channels) {
	json json_sig = json({ { "code", signalName },{ "data", json::array() } });
	if (units_tbl.count(signalName) != 0)
		json_sig["unit"] = units_tbl.at(signalName);
	int nelem = 0;
	if (TimeStamps_len != 0)
		nelem = TimeStamps_len / n_time_channels;
	else nelem = Values_len / n_val_channels;
	for (int i = 0; i < nelem; i++) {
		json json_sig_data_item = json({ { "timestamp" , json::array() },{ "value" , json::array() } });

		for (int j = 0; j < n_time_channels; j++) {
			json_sig_data_item["timestamp"].push_back(*TimeStamps);
			TimeStamps++;
		}

		for (int j = 0; j < n_val_channels; j++) {
			json_sig_data_item["value"].push_back(*Values);
			Values++;
		}

		json_sig["data"].push_back(json_sig_data_item);
	}
	return json_sig;
}

json CommonTestingTools::json_AddDataStr(const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, char** Values, int n_time_channels, int n_val_channels) {
	json json_sig = json({ { "code", signalName },{ "data", json::array() } });
	if (units_tbl.count(signalName) != 0)
		json_sig["unit"] = units_tbl.at(signalName);
	int nelem = 0;
	if (TimeStamps_len != 0)
		nelem = TimeStamps_len / n_time_channels;
	else nelem = Values_len / n_val_channels;
	for (int i = 0; i < nelem; i++) {
		json json_sig_data_item = json({ { "timestamp" , json::array() },{ "value" , json::array() } });

		for (int j = 0; j < n_time_channels; j++) {
			json_sig_data_item["timestamp"].push_back(*TimeStamps);
			TimeStamps++;
		}

		for (int j = 0; j < n_val_channels; j++) {
			json_sig_data_item["value"].push_back(*Values);
			Values++;
		}

		json_sig["data"].push_back(json_sig_data_item);
	}
	return json_sig;
}


