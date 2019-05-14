#include "MedSparseMat.h"
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <Logger/Logger/Logger.h>
#include <omp.h>
#include <mutex>

#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL	LOG_DEF_LEVEL

using namespace std;

mutex m_sparse_mat;


//======================================================================================================
void MedSparseMat::convert_map_to_line(map<int, float> &_map, vector<pair<int, float>> &line)
{
	line.clear();
	for (auto & e : _map) line.push_back(pair<int, float>(e.first, e.second));
}

//======================================================================================================
void MedSparseMat::convert_map_to_line(map<int, int> &_map, vector<pair<int, float>> &line)
{
	line.clear();
	for (auto & e : _map) line.push_back(pair<int, float>(e.first, (float)e.second));
}


//======================================================================================================
void MedSparseMat::insert_line(SparseMatRowMetaData &_meta, int line_num, vector<pair<int, float>> &line)
{
	{
		lock_guard<mutex> guard(m_sparse_mat);
		if (line_num < 0) line_num = (int)lines.size();
		if (lines.size() < (line_num+1)) {
			lines.resize(line_num+1); meta.resize(line_num+1);
		}
	}
	lines[line_num] = line;
	meta[line_num] = _meta;
}


//======================================================================================================
void MedSparseMat::insert_line(SparseMatRowMetaData &_meta, map<int, float> &_line)
{
	vector<pair<int, float>> line;
	convert_map_to_line(_line, line);
	insert_line(_meta, -1, line);
}

//======================================================================================================
void MedSparseMat::insert_line(SparseMatRowMetaData &_meta, map<int, int> &_line)
{
	vector<pair<int, float>> line;
	convert_map_to_line(_line, line);
	insert_line(_meta, -1, line);
}

//======================================================================================================
int MedSparseMat::write_to_files(string mat_file, string meta_file, string dict_file)
{
	// writing mat file: csv triplets: row, col, data
	if (mat_file != "" && lines.size() > 0) {
		ofstream mat_f;

		mat_f.open(mat_file);

		if (!mat_f.is_open()) {
			MERR("Can't open file %s for writing\n", mat_file.c_str());
			return -1;
		}

		mat_f << "row,col,data\n";
		for (int i=0; i<lines.size(); i++) {
			for (auto &p : lines[i]) {
				int ival = (int)p.second;
				if (p.second == (float)ival)
					mat_f << i << "," << p.first << "," << ival << "\n";
				else
					mat_f << i << "," << p.first << "," << p.second << "\n";
			}
		}
		mat_f.close();

		//write_to_bin_file(mat_file+"_bin");
	}

	// writing meta file
	if (meta_file != "" && meta.size() > 0) {
		ofstream meta_f;
		meta_f.open(meta_file);
		if (!meta_f.is_open()) {
			MERR("Can't open file %s for writing\n", meta_file.c_str());
			return -1;
		}
		meta_f << "line_num,pid,time\n";
		for (int i=0; i<meta.size(); i++) {
			meta_f << i << "," << meta[i].pid << "," << meta[i].time << "\n";
		}
		meta_f.close();
	}

	// writing dictionary
	if (dict_file != "" && dict.size() > 0) {
		ofstream dict_f;
		dict_f.open(dict_file);
		if (!dict_f.is_open()) {
			MERR("Can't open file %s for writing\n", dict_file.c_str());
			return -1;
		}
		dict_f << "val\tname\n";
		for (auto &s : dict)
			dict_f << s.first << "\t" << s.first << "\t" << s.second << "\n";
		dict_f.close();
	}

	return 0;
}

int MedSparseMat::read_from_files(string mat_file, string meta_file)
{
	lines.clear();
	meta.clear();

	// read mat file
	ifstream mat_f(mat_file);

	if (!mat_f) {
		MERR("MedSparseMat: read: Can't open file %s\n", mat_file.c_str());
		return -1;
	}

	string curr_line;
	while (getline(mat_f, curr_line)) {
		if ((curr_line.size() > 1) && (curr_line[0] != 'r')) {
			vector<string> fields;
			boost::split(fields, curr_line, boost::is_any_of(","));
			int row = stoi(fields[0]);
			int col = stoi(fields[1]);
			int val = stof(fields[2]);
			if (lines.size() < row+1) lines.resize(row+1);
			lines[row].push_back(pair<int, float>(col, (float)val));
		}
	}

	mat_f.close();

	// read meta file
	ifstream meta_f(meta_file);

	if (!meta_f) {
		MERR("MedSparseMat: read: Can't open file %s\n", meta_file.c_str());
		return -1;
	}

	while (getline(meta_f, curr_line)) {
		if ((curr_line.size() > 1) && (curr_line[0] != 'l')) {
			vector<string> fields;
			boost::split(fields, curr_line, boost::is_any_of(","));
			int pid= stoi(fields[1]);
			int time = stoi(fields[2]);
			SparseMatRowMetaData m;
			m.pid = pid;
			m.time = time;
			meta.push_back(m);
		}
	}

	meta_f.close();
	return 0;
}


// the following writes only the mat to a bin file with 3 arrays of len for row , col and val (each with length at start).
// this allows for an easier and faster transfer to python.
int MedSparseMat::write_to_bin_file(string bin_file)
{

	vector<int> rows, cols;
	vector<float> vals;

	for (int r=0; r<lines.size(); r++)
		for (auto &p : lines[r]) {
			rows.push_back(r);
			cols.push_back(p.first);
			vals.push_back(p.second);
		}


	MLOG("NedSparseMat::write_to_bin_file : rows %d cols %d vals %d\n", (int)rows.size(), (int)cols.size(), (int)vals.size());
	size_t size = MedSerialize::get_size(rows, cols, vals);
	MLOG("NedSparseMat::write_to_bin_file : get_size %d\n", (int)size);
	vector<unsigned char> data(size);
	MedSerialize::serialize(&data[0], rows, cols, vals);
	return MedSerialize::write_binary_data(bin_file, &data[0], size);
}


//---------------------------------------------------------------------------------------------------------------------------------------
void MedSparseMat::get_col_stat(int &nrows, int &ncols, vector<int> &nonz_counts)
{
	ncols = 0;
	nrows = (int)lines.size();
	nonz_counts.resize(10000000, 0); // assuming always less than 10M columns

	for (auto &line : lines) {
		for (auto &e : line) {
			if (e.second != 0)
			nonz_counts[e.first]++;
			if (e.first > ncols) ncols = e.first;
		}
	}
	ncols++;
	nonz_counts.resize(ncols);
}

//---------------------------------------------------------------------------------------------------------------------------------------
int MedSparseMat::write_col_stat_file(string f_stat)
{
	ofstream stat_f;

	stat_f.open(f_stat);

	if (!stat_f.is_open()) {
		MERR("Can't open file %s for writing\n", f_stat.c_str());
		return -1;
	}

	int ncols, nrows;
	vector<int> nz_counts;
	get_col_stat(nrows, ncols, nz_counts);
	stat_f << "col,count,n_lines,prob\n";
	for (int i = 0; i < nz_counts.size(); i++)
		stat_f << i << "," << nz_counts[i] << "," << nrows << "," << (double)nz_counts[i] / (double)nrows << "\n";
	stat_f.close();
	return 0;
}