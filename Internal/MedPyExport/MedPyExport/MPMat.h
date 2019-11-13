#ifndef __MED__MPMAT__H__
#define __MED__MPMAT__H__

#include "MedPyCommon.h"
#include "MPSerializableObject.h"

template <class T>
class MedMat;

class MPMat {
	bool o_owned = true;
public:
#ifndef SWIG
	MedMat<float>* o;
	MPMat(MedMat<float>* from_ptr);
	MPMat(const MPMat& other);
#endif

	static const int Normalize_Cols;
	static const int Normalize_Rows;
	static const int MISSING_VALUE;

	MPMat();
	MPMat(int n_rows, int n_cols);
	MPMat(float *IN_ARRAY2, unsigned long long DIM1, unsigned long long DIM2);

	~MPMat();
	
	/*NumPy export/import*/
	
	/* "Unsafe" means the same buffer is being used for both the NumPy array and the MedMat.
	 * This, as you might have guessed, can cause the segfault-making kind of issues when data is manipulated and changed.
	 */
	
	void get_numpy_copy(float** ARGOUTVIEWM_ARRAY2, unsigned long long* DIM1, unsigned long long* DIM2);
	void get_numpy_view_unsafe(float** ARGOUTVIEW_ARRAY2, unsigned long long* DIM1, unsigned long long* DIM2);
	void test(float** ARGOUTVIEWM_ARRAY2, unsigned long long* DIM1, unsigned long long* DIM2, int n_rows, int n_cols);
	void load_numpy(float* IN_ARRAY2, unsigned long long DIM1, unsigned long long DIM2);

	static int MEDPY_GET_Normalize_Cols();
	static int MEDPY_GET_Normalize_Rows();

	int MEDPY_GET_nrows();
	int MEDPY_GET_ncols();
	unsigned long long __len__();

	bool is_valid(bool output = false);

	float __getitem__(const std::vector<int>& index);
	void __setitem__(const std::vector<int>& index, float val);

	void MEDPY_GET_row_ids(MEDPY_NP_OUTPUT(int** out_row_ids_buf, unsigned long long* out_row_ids_buf_len));
	void MEDPY_SET_row_ids(MEDPY_NP_INPUT(int* row_ids_buf, unsigned long long row_ids_buf_len));

	std::vector<std::string> MEDPY_GET_signals();
	
	void MEDPY_GET_avg(MEDPY_NP_OUTPUT(int** avg_buf, unsigned long long* avg_buf_len));
	void MEDPY_GET_std(MEDPY_NP_OUTPUT(int** std_buf, unsigned long long* std_buf_len));

	int MEDPY_GET_normalized_flag();
	void MEDPY_SET_normalized_flag(int newval);

	int MEDPY_GET_transposed_flag();
	void MEDPY_SET_transposed_flag(int newval);

	float MEDPY_GET_missing_value();
	void MEDPY_SET_missing_value(float newval);

	void zero();
	void set_val(float val);
	void clear();

	void resize(int n_rows, int n_cols);

	int read_from_bin_file(const string &fname);
	int write_to_bin_file(const string &fname);
	int write_to_csv_file(const string &fname);
	int read_from_csv_file(const string &fname, int titles_line_flag);

	void transpose();

	void get_sub_mat(vector<int> &rows_to_take, vector<int> &cols_to_take);
	void get_sub_mat_by_flags(vector<int> &rows_to_take_flag, vector<int> &cols_to_take_flag);
	//void reorder_by_row(vector<int> &row_order);
	//void reorder_by_col(vector<int> &col_order);

	void add_rows(MPMat& m_add);
	void add_rows(MEDPY_NP_INPUT(float* m_add, unsigned long long nrows_to_add));
	void add_cols(MPMat& m_add);
	void add_cols(MEDPY_NP_INPUT(float* m_add, unsigned long long ncols_to_add));

	// get a row or a column to a vector
	void get_row(int i_row, MEDPY_NP_OUTPUT(float** rowv, unsigned long long* rowv_n)) const;
	void get_col(int i_col, MEDPY_NP_OUTPUT(float** colv, unsigned long long* colv_n)) const;

	void normalize(int norm_type, MEDPY_NP_INPUT(float* wgts, unsigned long long wgts_n));
	void normalize(int norm_type = Normalize_Cols);

	void normalize(MEDPY_NP_INPUT(float* external_avg, unsigned long long external_avg_n), MEDPY_NP_INPUT(float* external_std, unsigned long long external_std_n), int norm_type = 1);

	void get_cols_avg_std(MEDPY_NP_OUTPUT(float** buf_avg, unsigned long long* buf_avg_n), MEDPY_NP_OUTPUT(float** buf_std, unsigned long long* buf_std_n));
	void set_signals(std::vector<std::string> & sigs);

	void load(float *IN_ARRAY2, unsigned long long DIM1, unsigned long long DIM2);
	void load_transposed(float *IN_ARRAY2, unsigned long long DIM1, unsigned long long DIM2);
	void load(MPMat& x);
	MPSerializableObject asSerializable();

/********
 *
 *  Missing stuff:
 *
	// metadata holders
	vector<RecordData> recordsMetadata;
	template <class S> MedMat(vector<S> &x, int n_cols) { clear(); load(x, n_cols); }
	template <class S> void load(vector<S> &x, int n_cols);
	T *data_ptr() { if (m.size()>0) return &m[0]; else return NULL; }
	T *data_ptr(int r, int c) { if (m.size()>r*ncols + c) return &m[(unsigned long long)r*ncols + c]; else return NULL; }
	void print_row(FILE *fout, const string &prefix, const string &format, int i_row);
 *
 *
**/

};


#endif //!__MED__MPMAT__H__
