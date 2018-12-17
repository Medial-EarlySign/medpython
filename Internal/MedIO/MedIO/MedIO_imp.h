// implementation for template functions in MedIO.h

#include <fstream>
#include <cassert>

//#define LOCAL_SECTION LOG_MEDIO
//#define LOCAL_LEVEL	LOG_DEF_LEVEL
//extern MedLogger global_logger;

//-----------------------------------------------------------------------------
template <class T> int write_vector(const string &fname, vector<T> &data)
{
	ofstream of;

	//	MLOG_D("MedIO: write_vector: fname %s data size %d\n",fname.c_str(),(int)data.size());
	of.open(fname, ios::out | ios::binary);
	if (!of)
		return -1;

	if (data.size() > 0) {
		char *d = (char *)(&data[0]);
		//		MLOG_D("MedIO: writing %d bytes to file %s\n",sizeof(T)*data.size(),fname.c_str());
		of.write(d, sizeof(T)*data.size());
	}

	of.close();
	return 0;
}

//-----------------------------------------------------------------------------
template <class T> int read_vector(const string &fname, unsigned long long start_pos, vector<T> &data)
{
	ifstream inf;
	unsigned long long size;

	data.clear();
	size = get_file_size(fname);

	if (start_pos >= size)
		return 0;

	inf.open(fname, ios::in | ios::binary);
	if (!inf)
		return -1;

	inf.seekg(start_pos, ios_base::beg);

	size = size - start_pos;
	if (size > 0) {
		data.resize(size / (sizeof(T)));

		char *d = (char *)&data[0];
		inf.read(d, size);
	}
	return 0;
}

//-----------------------------------------------------------------------------
template <class T> int read_vector(const string &fname, vector<T> &data)
{
	return(read_vector(fname, 0, data));
}

//-----------------------------------------------------------------------------
//a class for reading separated files (csv, tsv, etc.) with mixed data types
//the tuple of specified types should be provided as a template type. Integral types are checked for correctness/missingness
//and replaced with special values accordingly: incorrect as the type's min(), missing as the type's max(). all fields are trimmed for white space.
//int the untyped case, all the fields are read as strings
//in both cases, separators inside quotes are ignored. quotes inside quotes are expected to be doubled.

//Example use:

//with type specification:
/*
using Type = tuple<int, int, int, string, string, int, int, string, string, int, string, string, string, string>;
string fname = "//nas1/Data/eicu/pastHistory.csv";
int monitorStep = 100000;
string fname = //nas1/Data/eicu/pastHistory.csv;
SepFileReader<Type> sr(fname, ',', 0, true, monitorStep);
sr.read();
const vector<Type>& data = sr.getData();
*/

//without type specification:
/*
SepFileReader<string> sr1(fname, ',', 0, true, monitorStep);
sr1.read();
vector<vector<string> > data1 = sr1.getData();
*/
//---------- some helper classes first ----------

template<typename Tuple>
class SepFileReader;

//helper class for _read
//read a single position in a single record string into an element of a tuple, given a pair of positions 
//where a field starts and ends
template <typename Tuple, size_t N>
struct handleSinglePos {
	void operator()(const string& str, const pair<int, int>& pos1, Tuple& e) {

		pair<int, int> pos = pos1;

		//trim whitespace
		const char *buf = str.c_str();
		while (pos.first <= pos.second && isspace(buf[pos.first]))
			pos.first++;

		while (pos.first <= pos.second && isspace(buf[pos.second]))
			pos.second--;

		//handle missing value
		if (pos.first > pos.second)
			SepFileReader<Tuple>::setMissingVal(get<N>(e));
		else {
			stringstream s(str.substr(pos.first, pos.second - pos.first + 1));
			s >> get<N>(e);

			if (s.fail())
				SepFileReader<Tuple>::setBadVal(get<N>(e));
		}
	};
};

//helper class for _read
//read a single record string into a tuple, given pairs of positions 
//where fields start and end
template <typename Tuple, size_t N>
struct handleSepRecord {
	void operator()(const string& str, const vector<pair<int, int> >& pos, Tuple& e) {
		handleSinglePos<Tuple, N>()(str, pos[N], e);
		handleSepRecord<Tuple, N - 1>()(str, pos, e);
	};
};

//specializetion of handleSepRecord for N=0
template <typename Tuple>
struct handleSepRecord<Tuple, 0> {
	void operator()(const string& str, const vector<pair<int, int> >& pos, Tuple& e) {
		handleSinglePos<Tuple, 0>()(str, pos[0], e);
	};
};

//------------- the main class -------------
template <typename Tuple>
class SepFileReader {
public:
	//constructor:
	//sep - separator character
	//maxRecords - max records to read. default 0 (all).
	//skipHeader - whether or not to skip a single header line
	//monitorStep - each <monitorStep> records, print out a '.' to monitor progress
	SepFileReader(const string& fname1, char sep1 = ',', size_t maxRecords1 = 0, bool skipHeader1 = true, int monitorStep1 = 0) :
		fname(fname1), sep(sep1), maxRecords(maxRecords1), skipHeader(skipHeader1), monitorStep(monitorStep1) {}

	//read the file
	//returns - 0 on success, -1 on failure.
	int read() { return _read<tuple_size<Tuple>::value>(); }
	vector<Tuple>& getData() { return v; }
	//set variable to its type's missing value
	//for integral types it's the max(), otherwise nothing is done
	template <typename T>
	static void setMissingVal(T& val) {
		if (is_arithmetic<T>::value)
			val = numeric_limits<T>::max();
	}

	//for arithmetic types, set val to the min of the type, signifying bad val
	//for other types, do nothing
	template <typename T>
	static void setBadVal(T& val) {
		if (is_arithmetic<T>::value)
			val = numeric_limits<T>::min();
	}
private:
	vector<Tuple> v;
	string fname;
	char sep;
	size_t maxRecords;
	bool skipHeader;
	int monitorStep;

	//does main work of read()
	template<size_t N>
	int _read() {
		v.clear();
		const size_t DEFAULT_NUM_RECORDS = 1000000;

		//reserve proper space for v
		if (maxRecords == 0) {
			v.reserve(DEFAULT_NUM_RECORDS);
			maxRecords = numeric_limits<size_t>::max();
		}
		else
			v.reserve(maxRecords);

		//open file
		ifstream inf(fname);
		if (!inf.is_open()) {
			cerr << "Cannot open " << fname << " for reading" << endl;
			return -1;
		}

		std::string str;
		Tuple elem;
		vector<pair<int, int> > positions; //hold the positions of fields (begin and end)

										   //handle header line correctly
		if (skipHeader && !inf.eof())
			getline(inf, str);

		//go over lines
		while (!inf.eof() && v.size() < maxRecords) {

			//monitor progress
			if (monitorStep > 0 && v.size() % monitorStep == 0)
				cout << ".";

			//process a single line
			getline(inf, str);
			if (str.empty())
				break;

			const char *buf = str.c_str();
			positions.clear();

			int curStart = 0;
			bool inQuotes = false;

			for (int i = 0; i < str.length(); ++i) {
				if (buf[i] == '"') {
					inQuotes = !inQuotes;
					continue;
				}

				if (!inQuotes && buf[i] == sep) {
					positions.push_back(pair<int, int>(curStart, i - 1));
					curStart = i + 1;
				}
			}

			positions.push_back(pair<int, int>(curStart, (int)(str.length()) - 1));
			assert(positions.size() == N);

			handleSepRecord<Tuple, N - 1>()(str, positions, elem);
			v.push_back(elem);
		}

		if (monitorStep > 0)
			cout << endl;

		inf.close();
		v.shrink_to_fit();

		return 0;
	}
};

//a specialization of SepFileReader that reads separated files 
//treating all fields read as strings. The data is held as 
//a vector<vector<string> > rather than a vector of a specified tuple
//this setting is preferable if the composition and number of the fields
//is not important. no missing/correctness checks are performed
//
//note: this specialization is quicker, and more memory-expensive.
template<>
class SepFileReader<string> {
public:
	//constructor
	SepFileReader(const string& fname1, char sep1 = ',', size_t maxRecords1 = 0, bool skipHeader1 = true, int monitorStep1 = 0) :
		fname(fname1), sep(sep1), maxRecords(maxRecords1), skipHeader(skipHeader1), monitorStep(monitorStep1) {}

	//count the number of fields in the file, using the first line
	//if file does not open return -1. if zero lines, return -2.
	int getNumberFieldsDirectly() const {
		ifstream inf(fname);
		if (!inf.is_open())
			return -1;

		if (inf.eof())
			return -2;

		string str;
		getline(inf, str);
		inf.close();

		if (str.empty())
			return 0;

		const char *buf = str.c_str();
		int count = 0;
		bool inQuotes = false;

		for (int i = 0; i < str.length(); ++i) {
			if (buf[i] == '"') {
				inQuotes = !inQuotes;
				continue;
			}

			if (!inQuotes && buf[i] == sep)
				count++;
		}

		return count + 1;
	}

	//return const reference to the data read
	vector<vector<string> > & getData() { return v; }

	//read the data
	int read() {
		int numFields = getNumberFieldsDirectly();

		if (numFields < 1)
			return -1;

		v.clear();
		const size_t DEFAULT_NUM_RECORDS = 1000000;

		if (maxRecords == 0) {
			v.reserve(DEFAULT_NUM_RECORDS);
			maxRecords = numeric_limits<size_t>::max();
		}
		else
			v.reserve(maxRecords);

		ifstream inf(fname);
		if (!inf.is_open()) {
			cerr << "Cannot open " << fname << " for reading" << endl;
			return -1;
		}

		string str;
		vector<string> record; //a single record
		record.reserve(numFields);

		//handle header line correctly
		if (skipHeader && !inf.eof())
			getline(inf, str);

		//go over records
		while (!inf.eof() && v.size() < maxRecords) {
			if (monitorStep > 0 && v.size() % monitorStep == 0)
				cout << ".";
			//process a single line
			getline(inf, str);
			if (str.empty())
				break;

			const char *buf = str.c_str();
			record.clear();

			int curStart = 0;
			bool inQuotes = false;

			for (int i = 0; i < str.length(); ++i) {
				if (buf[i] == '"') {
					inQuotes = !inQuotes;
					continue;
				}

				if (!inQuotes && buf[i] == sep) {
					record.push_back(str.substr(curStart, i - curStart));
					curStart = i + 1;
				}
			}

			record.push_back(str.substr(curStart, str.length() - curStart));

			//trim white space
			for (auto& s : record) {
				int i = 0, j = (int)(s.length()) - 1;
				const char *buf = s.c_str();
				while (i <= j && isspace(buf[i])) i++;
				while (i <= j && isspace(buf[j])) j--;
				s = s.substr(i, j - i + 1);
				s.shrink_to_fit();
			}

			assert(record.size() == numFields);
			v.push_back(record);
		}

		if (monitorStep > 0)
			cout << endl;

		inf.close();
		v.shrink_to_fit();

		return 0;
	}

private:
	vector<vector<string> > v;
	string fname;
	char sep;
	size_t maxRecords;
	bool skipHeader;
	int monitorStep;
};

