#define _CRT_SECURE_NO_WARNINGS
#include "MedUtils.h"
#include <unordered_set>
#include <sstream>
#include "MedGitVersion.h"

#define LOCAL_SECTION LOG_MED_UTILS
#define LOCAL_LEVEL	LOG_DEF_LEVEL

string run_current_path = "";

template<typename T> void medial::sort_ops::get_sort_indexes(const vector<T> &x, const bool descending_order, bool const abs_val, vector<int> &indexes)
{
	vector<pair<T, int>> v(x.size());

	if (abs_val) {
		for (int i = 0; i < x.size(); i++) {
			v[i].first = abs(x[i]);
			v[i].second = i;
		}
	}
	else {
		for (int i = 0; i < x.size(); i++) {
			v[i].first = x[i];
			v[i].second = i;
		}
	}

	if (descending_order) {
		sort(v.begin(), v.end(), [](const pair<T, int> &v1, const pair<T, int> &v2) {return (v1.first > v2.first); });
	}
	else {
		sort(v.begin(), v.end(), [](const pair<T, int> &v1, const pair<T, int> &v2) {return (v1.first < v2.first); });
	}

	indexes.resize(x.size());
	for (int i = 0; i < x.size(); i++)
		indexes[i] = v[i].second;

}

template void medial::sort_ops::get_sort_indexes<double>(const vector <double> &x, const bool descending_order, bool const abs_val, vector<int> &indexes);
template void medial::sort_ops::get_sort_indexes<float>(const vector<float> &x, const bool descending_order, bool const abs_val, vector<int> &indexes);
template void medial::sort_ops::get_sort_indexes<int>(const vector<int> &x, const bool descending_order, bool const abs_val, vector<int> &indexes);
template void medial::sort_ops::get_sort_indexes<short>(const vector<short> &x, const bool descending_order, bool const abs_val, vector<int> &indexes);

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

template<class T> void medial::print::print_vec(const vector<T> &vec, const string &title, const string &format, const string &delimeter) {
	if (vec.empty()) {
		MLOG("%s: EMPTY\n", title.c_str());
		return;
	}
	string bf = print_obj(vec[0], format);
	for (size_t i = 1; i < vec.size(); ++i)
		bf += delimeter + print_obj(vec[i], format);

	MLOG("%s: [%s]\n", title.c_str(), bf.c_str());
}
template void medial::print::print_vec<double>(const vector<double> &vec, const string &title, const string &format, const string &delimeter);
template void medial::print::print_vec<float>(const vector<float> &vec, const string &title, const string &format, const string &delimeter);
template void medial::print::print_vec<const char *>(const vector<const char *> &vec, const string &title, const string &format, const string &delimeter);
void medial::print::print_vec(const vector<string> &vec, const string &title, const string &delimeter) {
	if (vec.empty()) {
		MLOG("%s: EMPTY\n", title.c_str());
		return;
	}
	string bf = vec[0];
	for (size_t i = 1; i < vec.size(); ++i)
		bf += delimeter + vec[i];

	MLOG("%s: [%s]\n", title.c_str(), bf.c_str());
}

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

	if (uniq_vals.size() > prctile_samples->size()) {
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
template void medial::print::print_hist_vec<int>(const vector<int> &vec, const string &title, const string &format, const vector<double> *prctile_samples);

template<class T> int medial::process::binary_search_position(const vector<T> &v, T search, int start, int end) {
	int maxSize = (end - start) + 1;
	while (maxSize > 2) {
		int mid = start + int((maxSize - 1) / 2);
		if (search <= v[mid])
			end = mid;
		else
			start = mid;
		maxSize = (end - start) + 1;
	}

	//finish, when has 1-2 elements from start to end
	if (search <= v[start])
		return start;
	else if (search <= v[end])
		return end;
	else
		return end + 1;
}
template<class T> int medial::process::binary_search_position(const vector<T> &v, T search) {
	return binary_search_position(v, search, 0, (int)v.size() - 1);
}

template<typename T> int medial::process::binary_search_index(const T *begin, const T *end, T val) {
	int maxSize = (int)(end - begin) + 1;
	int pos = 0;
	while (maxSize > 2) {
		int mid = int((maxSize - 1) / 2);
		if (val < begin[mid])
			end = begin + mid;
		else if (val > begin[mid]) {
			begin = begin + mid;
			pos += mid;
		}
		else { //equals - find first positions
			if (begin[mid - 1] == val)
				end = begin + mid;
			else
				return pos + mid;
		}
		maxSize = (end - begin) + 1;
	}
	//maxSize <= 2
	if (*begin == val)
		return pos;
	else if (*end == val)
		return pos + 1;
	else
		return -1; //not found
}
template<typename T> int medial::process::binary_search_position(const T *begin, const T *end, T val, bool reversed) {
	int maxSize = (int)(end - begin) + 1;
	int pos = 0;
	while (maxSize > 2) {
		int mid = int((maxSize - 1) / 2);
		if (!reversed) {
			if (val <= begin[mid])
				end = begin + mid;
			else {
				begin = begin + mid;
				pos += mid;
			}
		}
		else {
			if (val >= begin[mid])
				end = begin + mid;
			else {
				begin = begin + mid;
				pos += mid;
			}
		}
		maxSize = (end - begin) + 1;
	}

	//maxSize <= 2) 
	if (!reversed) {
		if (val <= *begin)
			return pos;
		else if (val <= *end)
			return pos + 1;
		else
			return pos + maxSize;
	}
	else {
		if (val >= *begin)
			return pos;
		else if (val >= *end)
			return pos + 1;
		else
			return pos + maxSize;
	}



}
template<typename T> int medial::process::binary_search_position_last(const T *begin, const T *end, T val, bool reversed) {
	int maxSize = (int)(end - begin) + 1;
	int pos = 0;
	while (maxSize > 2) {
		int mid = int((maxSize - 1) / 2);
		if (!reversed) {
			if (val < begin[mid])
				end = begin + mid;
			else {
				begin = begin + mid;
				pos += mid;
			}
		}
		else {
			if (val > begin[mid])
				end = begin + mid;
			else {
				begin = begin + mid;
				pos += mid;
			}
		}
		maxSize = (end - begin) + 1;
	}

	//maxSize <= 2
	if (!reversed) {
		if (val < *begin)
			return pos;
		else if (val < *end)
			return pos + 1;
		else
			return pos + maxSize;
	}
	else {
		if (val > *begin)
			return pos;
		else if (val > *end)
			return pos + 1;
		else
			return pos + maxSize;
	}

}
template int medial::process::binary_search_position<int>(const int *begin, const int *end, int val, bool reversed);
template int medial::process::binary_search_position<double>(const double *begin, const double *end, double val, bool reversed);
template int medial::process::binary_search_position<float>(const float *begin, const float *end, float val, bool reversed);
template int medial::process::binary_search_position_last<int>(const int *begin, const int *end, int val, bool reversed);
template int medial::process::binary_search_position_last<double>(const double *begin, const double *end, double val, bool reversed);
template int medial::process::binary_search_position_last<float>(const float *begin, const float *end, float val, bool reversed);
template int medial::process::binary_search_position<int>(const vector<int> &v, int search);
template int medial::process::binary_search_position<double>(const vector<double> &v, double search);
template int medial::process::binary_search_position<float>(const vector<float> &v, float search);

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
	else if (a.value().type() == typeid(double)) {
		return to_string(a.as<double>());
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
		("help_module", po::value<string>(), "help on specific module")
		("base_config", po::value<string>(&base_config), "config file with all arguments - in CMD we override those settings")
		("debug", po::bool_switch(&debug), "set debuging verbose")
		("version", "prints version information of the program")
		;
	desc.add(general_options);
	desc.add(prg_options);
	if (!app_l.empty())
		app_logo = app_l;
	debug = false;
	init_called = true;
}

//finds module section help in full help message
string medial::io::ProgramArgs_base::get_section(const string &full_help, const string &search) {
	stringstream res;
	vector<string> lines;
	boost::split(lines, full_help, boost::is_any_of("\n"));
	bool in_section = false;
	for (size_t i = 0; i < lines.size(); ++i) {
		string ln = boost::trim_copy(lines[i]);
		if (!ln.empty() && ln.at(ln.length() - 1) == ':' && ln.substr(0, 2) != "--") {
			if (lines[i].find(search) != string::npos)
				in_section = true;
			else
				in_section = false;
		}
		if (in_section)
			res << lines[i] << "\n";
	}
	return res.str();
}

void medial::io::ProgramArgs_base::list_sections(const string &full_help, vector<string> &all_sec) {
	vector<string> lines;
	boost::split(lines, full_help, boost::is_any_of("\n"));
	for (size_t i = 0; i < lines.size(); ++i) {
		boost::trim(lines[i]);
		if (!lines[i].empty() && lines[i].at(lines[i].length() - 1) == ':' && lines[i].substr(0, 2) != "--")
			all_sec.push_back(lines[i].substr(0, lines[i].length() - 1));
	}
}

int medial::io::ProgramArgs_base::parse_parameters(int argc, char *argv[]) {
	if (!init_called)
		MTHROW_AND_ERR("ProgramArgs_base::init function wasn't called\n");

	po::options_description desc_file(desc);
	po::variables_map vm_config;

	auto parsed_args = po::parse_command_line(argc, argv, desc,
		po::command_line_style::style_t::default_style);

	po::store(parsed_args, vm);
	if (vm.count("help_module")) {
		string help_search = vm["help_module"].as<string>();
		stringstream help_stream;
		help_stream << desc;
		string full_help = help_stream.str();
		string module_help = get_section(full_help, help_search);
		if (module_help.empty()) {
			vector<string> all_sections;
			list_sections(full_help, all_sections);
			string section_msg = medial::io::get_list(all_sections, "\n");
			cout << "No help on search for module \"" << help_search << "\", Available Sections(" << all_sections.size() << "):\n" << section_msg << endl;
		}
		else
			cout << module_help << endl;

		return -1;
	}

	if (vm.count("help") || vm.count("h")) {
		MLOG("%s\n", app_logo.c_str());
		cout << desc << endl;
		return -1;
	}

	if (vm.count("version")) {
		MLOG("%s\n", app_logo.c_str());
		cout << "Version Info:\n" << medial::get_git_version() << endl;
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
	for (auto it = vm_config.begin(); it != vm_config.end(); ++it)
	{
		if (it->second.defaulted()) {
			continue;
		}
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
		string full_log_format = "$time\t$level\t$section\t%s";
		global_logger.init_format(LOG_APP, full_log_format);
		global_logger.init_format(LOG_DEF, full_log_format);
		global_logger.init_format(LOG_MED_MODEL, full_log_format);
		global_logger.init_format(LOG_MEDALGO, full_log_format);
		MLOG("Version Info:\n%s\n", medial::get_git_version().c_str());
		MLOG("Debug Running With:\n");
		string full_params = string(argv[0]);
		char buffer[1000];
		for (auto it = vm.begin(); it != vm.end(); ++it) {
			MLOG("%s = %s\n", it->first.c_str(), medial::print::print_any(it->second).c_str());
			string val = medial::print::print_any(it->second);
			//gets deafult value when defaulted
			//desc.find(it->first, true).semantic()
			if (it->second.value().type() == typeid(string))
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

template<class T> string medial::io::get_list(const unordered_map<string, T> &ls, const string &delimeter) {
	string res = "";
	for (auto it = ls.begin(); it != ls.end(); ++it)
		if (it == ls.begin())
			res += it->first;
		else
			res += delimeter + it->first;
	return res;
}
template string medial::io::get_list<int>(const unordered_map<string, int> &ls, const string &delimeter);
template<class ContainerType> string medial::io::get_list(const ContainerType &ls, const string &delimeter) {
	string res = "";
	for (auto it = ls.begin(); it != ls.end(); ++it)
		if (it == ls.begin())
			res += *it;
		else
			res += delimeter + *it;
	return res;
}
template string medial::io::get_list(const set<string> &ls, const string &delimeter);
template string medial::io::get_list(const unordered_set<string> &ls, const string &delimeter);
template string medial::io::get_list(const vector<string> &ls, const string &delimeter);

template<class T> string medial::io::get_list_op(const unordered_map<T, string> &ls, const string &delimeter) {
	string res = "";
	for (auto it = ls.begin(); it != ls.end(); ++it)
		if (it == ls.begin())
			res += it->second;
		else
			res += delimeter + it->second;
	return res;
}
template string medial::io::get_list_op<int>(const unordered_map<int, string> &ls, const string &delimeter);

void medial::print::log_with_file(ofstream &fw, const char *format_str, ...) {
	char buff[10000];
	va_list argptr;
	va_start(argptr, format_str);
	vsnprintf(buff, sizeof(buff), format_str, argptr);
	va_end(argptr);
	string final_str = string(buff);

	if (fw.is_open() && fw.good()) {
		fw << final_str;
		fw.flush();
	}
	MLOG("%s", final_str.c_str());
}

void medial::io::read_codes_file(const string &file_path, vector<string> &tokens) {
	tokens.clear();
	ifstream file;
	file.open(file_path);
	if (!file.is_open())
		MTHROW_AND_ERR("Unable to open test indexes file:\n%s\n", file_path.c_str());
	string line;
	//getline(file, line); //ignore first line
	while (getline(file, line)) {
		boost::trim(line);
		if (line.empty())
			continue;
		if (line.at(0) == '#')
			continue;
		if (line.find('\t') != string::npos)
			line = line.substr(0, line.find('\t'));
		tokens.push_back(line);
	}
	file.close();
}

string medial::get_git_version() {
	return GIT_HEAD_VERSION;
}