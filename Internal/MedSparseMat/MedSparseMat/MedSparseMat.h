#pragma once

// 
// A small lib to handle sparse matrices
// Mainly as input and use for embedding algorithms
//

#include <SerializableObject/SerializableObject/SerializableObject.h>

struct SparseMatRowMetaData : public SerializableObject {
	int pid;
	int time;
	ADD_CLASS_NAME(SparseMatRowMetaData)
	ADD_SERIALIZATION_FUNCS(pid, time)
};

// MedSparseMat is not templated ... currently float only
// Assuming the sparsness is on the cols dimension
class MedSparseMat : public SerializableObject {

public:
	int n_rows = 0;
	int n_cols = 0;

	// actual data
	vector<vector<pair<int, float>>> lines;

	// line meta data
	vector<SparseMatRowMetaData> meta;

	// dictionary for categories values
	map<int, string> dict;

	void clear() { lines.clear(); meta.clear(); dict.clear(); n_rows=0; n_cols=0; }

	void init(int _n_rows, int _n_cols) { n_rows = _n_rows; n_cols = _n_cols; lines.resize(n_rows); }

	// line_num < 0 is a signal to add as the next available line
	void insert_line(SparseMatRowMetaData &_meta, int line_num, vector<pair<int, float>> &line);
	void insert_line(SparseMatRowMetaData &_meta, map<int, float> &line);
	void insert_line(SparseMatRowMetaData &_meta, map<int, int> &line);

	static void convert_map_to_line(map<int, float> &_map, vector<pair<int, float>> &line);
	static void convert_map_to_line(map<int, int> &_map, vector<pair<int, float>> &line);


	void get_col_stat(int &nrows, int &ncols, vector<int> &nonz_counts);
	int write_col_stat_file(string f_stat);

	void insert_dict_item(int _val, string &_name) { dict[_val] = _name; }

	int write_to_files(string mat_file, string meta_file, string dict_file);

	int read_from_files(string mat_file, string meta_file);

	// the following writes only the mat to a bin file with len at the start and then 3 arrays of length len for row , col and val.
	// this allows for an easier and faster transfer to python.
	int write_to_bin_file(string bin_file);

	ADD_CLASS_NAME(MedSparseMat)
	ADD_SERIALIZATION_FUNCS(n_rows, n_cols, lines, meta, dict)
};

MEDSERIALIZE_SUPPORT(SparseMatRowMetaData)
MEDSERIALIZE_SUPPORT(MedSparseMat)
