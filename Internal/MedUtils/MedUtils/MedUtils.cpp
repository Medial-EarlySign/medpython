#define _CRT_SECURE_NO_WARNINGS
#include "MedUtils.h"

#define LOCAL_SECTION LOG_MED_UTILS
#define LOCAL_LEVEL	LOG_DEF_LEVEL


// Pearson Correlation: A pure C function
double  get_corr_pearson(float *v1, float *v2, int len)
{
	double sx, sy, sxy, sxx, syy, n;

	sx = sy = sxy = sxx = syy = 0;

	for (int i = 0; i < len; i++) {
		sx += v1[i];
		sy += v2[i];
		sxx += v1[i] * v1[i];
		syy += v2[i] * v2[i];
		sxy += v1[i] * v2[i];
	}

	n = (double)len;

	sx /= n;
	sy /= n;
	sxx /= n;
	syy /= n;
	sxy /= n;

	double c1 = sxy - sx*sy;
	double c2 = sxx - sx*sx;
	double c3 = syy - sy*sy;

	double epsilon = 1e-8;
	if (c2 < epsilon || c3 < epsilon)
		return 0;

	return (float)(c1 / (sqrt(c2)*sqrt(c3)));
}

// Mutual information of binned-vectors
int get_mutual_information(vector<int>& x, vector<int>& y, double& mi, int &n) {

	if (x.size() != y.size()) {
		MERR("Size mismatch. Quitting\n");
		return -1;
	}

	map<int, int> x_counts, y_counts;
	map<pair<int, int>, int> xy_counts;
	n = 0;

	for (unsigned int i = 0; i < x.size(); i++) {
		if (x[i] != -1 && y[i] != -1) {
			x_counts[x[i]]++;
			y_counts[y[i]]++;
			xy_counts[pair<int, int>(x[i], y[i])]++;
			n++;
		}
	}

	if (n < 2) {
		MLOG_V("Not enough common non-missing entries for mutual information.\n");
		mi = -1;
	}
	else {
		mi = get_mutual_information(x_counts, n, y_counts, n, xy_counts, n);
	}

	return 0;
}

// Counts from binned vectors
void get_counts(vector<int>& x, map<int, int>& counts, int& n) {

	n = 0;
	counts.clear();
	for (unsigned int i = 0; i < x.size(); i++) {
		if (x[i] != -1) {
			counts[x[i]] ++;
			n++;
		}
	}
}

int get_co_counts(vector<int>& x, vector<int>& y, map<pair<int, int>, int>& counts, int& n) {

	if (x.size() != y.size()) {
		MERR("Size mismatch (%d vs %d). Quitting\n", x.size(), y.size());
		return -1;
	}

	n = 0;
	counts.clear();
	for (unsigned int i = 0; i < x.size(); i++) {
		if (x[i] != -1 && y[i] != -1) {
			counts[pair<int, int>(x[i], y[i])] ++;
			n++;
		}
	}

	return 0;
}

// Mutual information from counts
double get_mutual_information(map<int, int>& x_count, int nx, map<int, int>& y_count, int ny, map<pair<int, int>, int>& xy_count, int n) {

	double mi = 0;
	for (auto it = xy_count.begin(); it != xy_count.end(); it++) {
		double p = (it->second + 0.0) / n;
		double px = (x_count[it->first.first] + 0.0) / nx;
		double py = (y_count[it->first.second] + 0.0) / ny;

		mi += p * log(p / px / py) / log(2.0);
	}

	return mi;
}

// Get moments of a vector
int get_moments(vector<float>& values, vector<float>& wgts, float missing_value, float& mean, float&sd, bool do_missing) {

	return get_moments(&(values[0]), &(wgts[0]), (int)values.size(), missing_value, mean, sd, do_missing);
}

int get_moments(float *values, float* wgts, int size, float missing_value, float& mean, float&sd, bool do_missing) {
	double  n = 0;
	double s = 0;

	for (int i = 0; i < size; i++) {
		if (!do_missing || values[i] != missing_value) {
			n += wgts[i];
			s += wgts[i] * values[i];
		}
	}
	if (n == 0) {
		mean = 0;
		sd = 1.0;
		return 0;
	}

	mean = (float)(s / n);

	s = 0.0;
	for (int i = 0; i < size; i++) {
		if (!do_missing || values[i] != missing_value)
			s += wgts[i] * (values[i] - mean)*(values[i] - mean);
	}

	if (n > 1)
		sd = (float)sqrt((s / n));
	else
		sd = (float) 1.0;

	if (sd == 0) {
		MWARN("get_moments for all-zeros vector, fixing SD from 0.0 to 1.0\n");
		sd = 1.0;
	}
	return (int)n;
}

template<class T> string medial::print::print_obj(T obj, const string &format) {
	//return to_string((round(num * 1000) / 1000));
	char res[50];
	snprintf(res, sizeof(res), format.c_str(), obj);
	return string(res);
}
template string medial::print::print_obj<float>(float obj, const string &format);
template string medial::print::print_obj<double>(double obj, const string &format);
template string medial::print::print_obj<const char *>(const char *obj, const string &format);
template string medial::print::print_obj<int>(int obj, const string &format);
template<class T> void medial::process::prctils(const vector<T> &x, const vector<double> &prc,
	vector<T> &res, const vector<float> *weights) {
	if (x.size() == 0)
		MTHROW_AND_ERR("x is empty\n");
	bool has_weights = weights != NULL && !weights->empty();
	if (has_weights && x.size() != weights->size())
		MTHROW_AND_ERR("x and weights are not same size\n");

	if (!has_weights) {
		vector<T> cp(x);
		T *data = cp.data();
		sort(cp.begin(), cp.end());

		res.resize((int)prc.size());
		for (size_t i = 0; i < res.size(); ++i)
		{
			double pos = prc[i] * (x.size() - 1);
			int pos_a = (int)pos;
			T r = data[pos_a];
			res[i] = r;
			if (pos_a + 1 < x.size())
				res[i] = T(res[i] * (1 - (pos - pos_a)) + (pos - pos_a) * data[pos_a + 1]);
		}
	}
	else {
		vector<pair<T, float>> cp(x.size());
		for (size_t i = 0; i < x.size(); ++i)
		{
			cp[i].first = x[i];
			cp[i].second = (*weights)[i];
		}
		pair<T, float> *data = cp.data();
		sort(cp.begin(), cp.end());
		vector<float> w(weights->size());
		w[0] = cp.front().second;
		for (size_t i = 1; i < cp.size(); ++i)
			w[i] = w[i - 1] + cp[i].second;
		float total_w = w.back();

		res.resize((int)prc.size());
		for (size_t i = 0; i < res.size(); ++i)
		{
			float pos = float(prc[i] * total_w);
			int pos_a = medial::process::binary_search_position(w.data(), w.data() + (int)w.size() - 1, pos);
			pos_a = min(pos_a, (int)cp.size() - 1);
			//int pos_a = (int)pos;
			T r = data[pos_a].first;
			res[i] = r;
			//if (pos_a + 1 < x.size())
			//	res[i] = T(res[i] * (1 - (pos - pos_a)) + (pos - pos_a) * data[pos_a + 1]);
		}
	}
}
template void medial::process::prctils<float>(const vector<float> &x, const vector<double> &prc, vector<float> &res, const vector<float> *weights);
template void medial::process::prctils<double>(const vector<double> &x, const vector<double> &prc, vector<double> &res, const vector<float> *weights);

template<class T> void medial::print::print_vec(const vector<T> &vec, const string &title, const string &format) {
	if (vec.empty()) {
		MLOG("%s: EMPTY\n", title.c_str());
		return;
	}
	string bf = print_obj(vec[0], format);
	for (size_t i = 1; i < vec.size(); ++i)
		bf += ", " + print_obj(vec[i], format);

	MLOG("%s: [%s]\n", title.c_str(), bf.c_str());
}
template void medial::print::print_vec<double>(const vector<double> &vec, const string &title, const string &format);
template void medial::print::print_vec<float>(const vector<float> &vec, const string &title, const string &format);

template<class T> void medial::print::print_hist_vec(const vector<T> &vec, const string &title,
	const string &format, const vector<double> *prctile_samples) {
	vector<double> default_prctiles = { 0, 0.1, 0.2, 0.3,0.4,0.5,0.6,0.7,0.8,0.9 ,1 };
	if (prctile_samples == NULL)
		prctile_samples = &default_prctiles;

	if (vec.size() <= 10)
		print_vec(vec, title, format);
	vector<T> prcs;
	map<T, int> uniq_vals;
	for (size_t i = 0; i < vec.size(); ++i)
		++uniq_vals[vec[i]];
	char res[500];

	if (uniq_vals.size() > 10) {
		medial::process::prctils(vec, *prctile_samples, prcs);
		snprintf(res, sizeof(res), ("%2.1f%%:" + format).c_str(), 100.0*(*prctile_samples)[0], prcs[0]);
		string bf = string(res);

		for (size_t i = 1; i < prcs.size(); ++i) {
			snprintf(res, sizeof(res), (", %2.1f%%:" + format).c_str(), 100.0*(*prctile_samples)[i], prcs[i]);
			bf += string(res);
		}
		MLOG("%s: HISTOGRAM[%s]\n", title.c_str(), bf.c_str());
	}
	else {
		auto ii = uniq_vals.begin();
		snprintf(res, sizeof(res), (format + ":%2.2f%%").c_str(), ii->first,
			100 * ii->second / double(vec.size()));
		string bf = string(res);
		++ii;
		for (; ii != uniq_vals.end(); ++ii) {
			snprintf(res, sizeof(res), (", " + format + ":%2.2f%%").c_str(), ii->first,
				100 * ii->second / double(vec.size()));
			bf += string(res);
		}
		MLOG("%s: VALUES[%s]\n", title.c_str(), bf.c_str());
	}
}
template void medial::print::print_hist_vec<double>(const vector<double> &vec, const string &title, const string &format, const vector<double> *prctile_samples);
template void medial::print::print_hist_vec<float>(const vector<float> &vec, const string &title, const string &format, const vector<double> *prctile_samples);

template<typename T> int medial::process::binary_search_index(const T *begin, const T *end, T val) {
	int maxSize = (int)(end - begin) + 1;
	int mid = int((maxSize - 1) / 2);
	if (maxSize <= 2) {
		if (*begin == val) {
			return 0;
		}
		else if (*end == val) {
			return 1;
		}
		else {
			return -1;
		}
	}

	if (begin[mid] == val) {
		//return first Index if there are more then one
		if (mid > 0 && begin[mid - 1] == val)
		{
			if (*begin == val) {
				return 0;
			}
			else {
				return binary_search_index(begin, begin + mid + 1, val);
			}
		}
		else {
			return mid;
		}
	}
	else {
		if (val > begin[mid]) {
			int p = binary_search_index(begin + mid, end, val);
			if (p > 0) {
				return mid + p;
			}
			else {
				return -1;
			}
		}
		else {
			return binary_search_index(begin, begin + mid, val);
		}
	}
}
template<typename T> int medial::process::binary_search_position(const T *begin, const T *end, T val, bool reversed) {
	int maxSize = (int)(end - begin) + 1;
	int mid = int((maxSize - 1) / 2);
	if (maxSize <= 2) {
		if (!reversed) {
			if (val < *begin) {
				return 0;
			}
			else if (val < *end) {
				return 1;
			}
			else {
				return maxSize;
			}
		}
		else {
			if (val >= *begin) {
				return 0;
			}
			else if (val >= *end) {
				return 1;
			}
			else {
				return maxSize;
			}
		}
	}

	if (!reversed) {
		if (val <= begin[mid]) {
			return binary_search_position(begin, begin + mid, val, reversed);
		}
		else {
			return mid + binary_search_position(begin + mid, end, val, reversed);
		}
	}
	else {
		if (val >= begin[mid]) {
			return binary_search_position(begin, begin + mid, val, reversed);
		}
		else {
			return mid + binary_search_position(begin + mid, end, val, reversed);
		}
	}
}
template int medial::process::binary_search_position<int>(const int *begin, const int *end, int val, bool reversed);
template int medial::process::binary_search_position<double>(const double *begin, const double *end, double val, bool reversed);
template int medial::process::binary_search_position<float>(const float *begin, const float *end, float val, bool reversed);

template int medial::process::binary_search_index(const float *begin, const float *end, float val);
template int medial::process::binary_search_index(const int *begin, const int *end, int val);
template int medial::process::binary_search_index(const double *begin, const double *end, double val);
template int medial::process::binary_search_index(const string *begin, const string *end, string val);

string medial::print::print_any(po::variable_value &a) {
	if (a.value().type() == typeid(string)) {
		return a.as<string>();
	}
	else if (a.value().type() == typeid(int)) {
		return to_string(a.as<int>());
	}
	else if (a.value().type() == typeid(float)) {
		return to_string(a.as<float>());
	}
	else if (a.value().type() == typeid(bool)) {
		return to_string(a.as<bool>());
	}

	return "";
}

void medial::io::ProgramArgs_base::init(po::options_description &prg_options, const string &app_l) {
	po::options_description general_options("Program General Options",
		(unsigned int)po::options_description::m_default_line_length * 2);
	general_options.add_options()
		("help,h", "help & exit")
		("base_config", po::value<string>(&base_config), "config file with all arguments - in CMD we override those settings")
		("debug", po::bool_switch(&debug), "set debuging verbose");
	desc.add(general_options);
	desc.add(prg_options);
	if (!app_l.empty())
		app_logo = app_l;
	debug = false;
	init_called = true;
}

int medial::io::ProgramArgs_base::parse_parameters(int argc, char *argv[]) {
	if (!init_called)
		MTHROW_AND_ERR("ProgramArgs_base::init function wasn't called\n");
	po::variables_map vm;
	po::options_description desc_file(desc);
	po::variables_map vm_config;

	auto parsed_args = po::parse_command_line(argc, argv, desc,
		po::command_line_style::style_t::default_style);
	po::store(parsed_args, vm);
	if (vm.count("help") || vm.count("h")) {
		MLOG("%s\n", app_logo.c_str());
		cout << desc << endl;
		return -1;
	}

	if (vm.count("base_config") > 0) {
		std::ifstream ifs(vm["base_config"].as<string>(), std::ifstream::in);
		if (!ifs.good())
			MTHROW_AND_ERR("IO Error: can't read \"%s\"\n", vm["base_config"].as<string>().c_str());
		auto parsed = po::parse_config_file(ifs, desc_file, true);
		po::store(parsed, vm_config);
		ifs.close();
	}
	//iterate on all values and override defaults in desc:
	bool has_valus = false;
	for (auto it = vm_config.begin(); it != vm_config.end(); ++it)
	{
		if (it->second.defaulted()) {
			continue;
		}
		has_valus = true;
		if (vm.find(it->first) == vm.end() || vm[it->first].defaulted()) {
			//should not happended

			if (vm.find(it->first) == vm.end()) {
				vm.insert(pair<string, po::variable_value>(it->first, it->second));
			}
			vm.at(it->first) = it->second;

		}
	}

	po::notify(vm);

	post_process();

	if (debug) {
		MLOG("Debug Running With:\n");
		string full_params = string(argv[0]);
		char buffer[1000];
		for (auto it = vm.begin(); it != vm.end(); ++it) {
			MLOG("%s = %s\n", it->first.c_str(), medial::print::print_any(it->second).c_str());
			string val = medial::print::print_any(it->second);
			if (val.empty() || val.find(";") != string::npos)
				val = "\"" + val + "\"";
			snprintf(buffer, sizeof(buffer), " --%s %s", it->first.c_str(), val.c_str());
			full_params += string(buffer);
		}
		MLOG("######################################\n\n%s\n", app_logo.c_str());
		MLOG("######################################\n");
		MLOG("Full Running Command:\n%s\n", full_params.c_str());
		MLOG("######################################\n");
	}

	return 0;
}

template<class T> string medial::io::get_list(const unordered_map<string, T> &ls) {
	string res = "";
	for (auto it = ls.begin(); it != ls.end(); ++it)
		if (it == ls.begin())
			res += it->first;
		else
			res += "," + it->first;
	return res;
}
template string medial::io::get_list<int>(const unordered_map<string, int> &ls);

template<class T> string medial::io::get_list_op(const unordered_map<T, string> &ls) {
	string res = "";
	for (auto it = ls.begin(); it != ls.end(); ++it)
		if (it == ls.begin())
			res += it->second;
		else
			res += "," + it->second;
	return res;
}
template string medial::io::get_list_op<int>(const unordered_map<int, string> &ls);

float med_stof(const string& _Str) {
	try {
		return stof(_Str);
	}
	catch (exception e) {
		MTHROW_AND_ERR("invalid stof argument [%s]\n", _Str.c_str());
	}
}

int med_stoi(const string& _Str) {
	try {
		return stoi(_Str);
	}
	catch (exception e) {
		MTHROW_AND_ERR("invalid stoi argument [%s]\n", _Str.c_str());
	}
}
