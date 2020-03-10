//
// General routines for easier manipulation of a general 2d matrix
// To be used to hold data points, split them (train, test, cv), 
// and create features matrices.
//
// Besides the major MedMat class, contains also several routines
// to handle Matrices and vectors.
//

#ifndef __MED_MAT_H__
#define __MED_MAT_H__

#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif



#include <vector>
#include <math.h>
#include <SerializableObject/SerializableObject/SerializableObject.h>
#include <MedIO/MedIO/MedIO.h>
#include <MedTime/MedTime/MedTime.h>
#include <MedMat/MedMat/MedMatConstants.h>
#include <MedUtils/MedUtils/MedGlobalRNG.h>

using namespace std;

/// TODO: this class should be refactored and merged with MedSample
class RecordData : public SerializableObject {
public:
	RecordData() {};
	RecordData(int id, int date, long outcomeTime, int split, float weight, float label, float pred) {
		this->id = id;
		this->date = date;
		this->outcomeTime = outcomeTime;
		this->split = split;
		this->weight = weight;
		this->label = label;
		this->pred = pred;
	}
	int id;
	int date;
	int outcomeTime;
	int split;
	float label;

	float pred;
	float weight;

	ADD_CLASS_NAME(RecordData)
		ADD_SERIALIZATION_FUNCS(id, date, outcomeTime, split, label, pred, weight)
};


// General Class for a 2d Marix
// The matrix can contain elements of any type <T> (typically <T> is float or double)
// There are several basic methods to construct, load, and append to a matrix (even from a general other type <S>)
// Also reordering, choosing a sub matrix or transposing (while keeping the transpose state) are possible.
// One can serialize/deserialize a matrix and then use other utils to IO it, or use a direct method.
// mat(i,j) can be used to access the (i,j) element in the matrix for read or write.
template <class T>
class MedMat : public SerializableObject {
private:
	vector<T> m;
public:

	const static int Normalize_Cols = 1;
	const static int Normalize_Rows = 2;

	// data holders (major)
	//vector<T> m;
	unsigned long long nrows = 0;
	unsigned long long ncols = 0;
	unsigned long long size() const { return (unsigned long long)nrows*ncols; }
	vector<T> &get_vec() { return m; }
	void copy_vec(vector<T> &c_vec) { c_vec.clear(); c_vec = m; }

	// metadata holders
	vector<int> row_ids;
	vector<RecordData> recordsMetadata;
	vector<string> signals;
	//int time_unit;

	// Normalization/DeNormalization
	vector<T> avg;
	vector<T> std;

	int normalized_flag = 0; // 0 - non normalized 1 - cols normalized, 2 - rows normalized
	int transposed_flag = 0; // 0 - as was when matrix was created/loaded, 1 - transpose of it

	T missing_value;

	// get/set
	inline T operator ()(size_t i, size_t j) const { return m[i*ncols + j]; }

	inline T &operator ()(size_t i, size_t j) { return m[i*ncols + j]; }

	inline T get(size_t i, size_t j) const { return m[i*ncols + j]; }

	inline T& set(size_t i, size_t j) { return m[i*ncols + j]; }  // use var.set(i,j) = .... 

	inline void copy_header(MedMat<T> &other) {
		if (this != &other) {
			nrows = other.nrows;
			ncols = other.ncols;
			signals = other.signals;
			avg = other.avg;
			std = other.std;
			normalized_flag = other.normalized_flag;
			transposed_flag = other.transposed_flag;
			missing_value = other.missing_value;
		}
	}

	inline MedMat<T>& operator=(MedMat<T> &&other) noexcept {
		if (this != &other) {
			m = move(other.m);
			nrows = other.nrows;
			ncols = other.ncols;
			row_ids = move(other.row_ids);
			signals = move(other.signals);
			avg = move(other.avg);
			std = move(other.std);
			normalized_flag = other.normalized_flag;
			transposed_flag = other.transposed_flag;
			missing_value = other.missing_value;
			recordsMetadata = move(other.recordsMetadata);
		}
		return *this;
	}
	inline MedMat<T>& operator=(const MedMat<T> &other) noexcept {
		if (this != &other) {
			m = other.m;
			nrows = other.nrows;
			ncols = other.ncols;
			row_ids = other.row_ids;
			signals = other.signals;
			avg = other.avg;
			std = other.std;
			normalized_flag = other.normalized_flag;
			transposed_flag = other.transposed_flag;
			missing_value = other.missing_value;
			recordsMetadata = other.recordsMetadata;
		}
		return *this;
	}

	// init
	MedMat() { clear(); }
	MedMat(int n_rows, int n_cols) { clear(); nrows = n_rows; ncols = n_cols; m.resize(nrows*ncols); zero(); };
	MedMat(const MedMat<T> &other) {
		m = other.m;
		nrows = other.nrows;
		ncols = other.ncols;
		row_ids = other.row_ids;
		signals = other.signals;
		avg = other.avg;
		std = other.std;
		normalized_flag = other.normalized_flag;
		transposed_flag = other.transposed_flag;
		missing_value = other.missing_value;
		recordsMetadata = other.recordsMetadata;
	}
	template <class S> MedMat(S *x, int n_rows, int n_cols) { clear(); load(x, n_rows, n_cols); }
	template <class S> MedMat(const vector<S> &x, int n_cols) { clear(); load(x, n_cols); }
	template <class S> MedMat(MedMat<S> &x) { clear(); load(x); }

	template <class S> void load(S *x, int n_rows, int n_cols);
	template <class S> void load_transposed(S *x, int n_rows, int n_cols);
	template <class S> void load(const vector<S> &x, int n_cols);
	template <class S> void load(MedMat<S> &x);

	void zero() { fill(m.begin(), m.end(), (T)0); }
	void set_val(T val) { fill(m.begin(), m.end(), val); } // set all matrix to a certain value.

	// basic 
	void clear() { m.clear(); row_ids.clear(); signals.clear(); recordsMetadata.clear(); nrows = 0; ncols = 0; normalized_flag = 0; transposed_flag = 0; missing_value = (T)MED_MAT_MISSING_VALUE; }
	T *data_ptr() { if (m.size() > 0) return &m[0]; else return NULL; }
	const T *data_ptr() const { if (m.size() > 0) return &m[0]; else return NULL; }
	T *data_ptr(size_t r, size_t c) { if (m.size() > r*ncols + c) return &m[r*ncols + c]; else return NULL; }
	int get_nrows() { return (int)nrows; }
	int get_ncols() { return (int)ncols; }
	void resize(int n_rows, int n_cols) { nrows = n_rows; ncols = n_cols; m.resize(nrows*ncols); }

	// i/o from specific format files
	int read_from_bin_file(const string &fname);
	int write_to_bin_file(const string &fname);
	int write_to_csv_file(const string &fname);
	int read_from_csv_file(const string &fname, int titles_line_flag);
	//int read_from_csv_file(const string &fname, int titles_line_flag, vector<string>& fields_out);

	// serialize(), deserialize()
	//size_t get_size();
	//size_t serialize(unsigned char *buf);
	//size_t deserialize(unsigned char *buf);
	ADD_SERIALIZATION_FUNCS(m, nrows, ncols, row_ids, recordsMetadata, signals, avg, std, normalized_flag, transposed_flag, missing_value);

	// simple handling options
	void transpose();
	void get_sub_mat(vector<int> &rows_to_take, vector<int> &cols_to_take); // empty list means - take them  all
	void get_sub_mat_by_flags(vector<int> &rows_to_take_flag, vector<int> &cols_to_take_flag);
	void reorder_by_row(vector<int> &row_order);
	void reorder_by_col(vector<int> &col_order);

	void random_split_mat_by_ids(MedMat<T> &mat_0, MedMat<T> &mat_1, float p0, vector<int> &inds0, vector<int> &inds1);

	template <class S> void add_rows(MedMat<S> &m_add);
	template <class S> void add_rows(S *m_add, int nrows_to_add);
	template <class S> void add_rows(vector<S> &m_add);
	template <class S> void add_cols(MedMat<S> &m_add);
	template <class S> void add_cols(S *m_add, int ncols_to_add); // packed as nrows x ncols_to_add 
	template <class S> void add_cols(vector<S> &m_add);

	// get a row or a column to a vector
	void get_row(int i_row, vector<T> &rowv) const;
	void get_col(int i_col, vector<T> &colv) const;

	// normalization (norm_type = 1 for cols (default), 2 for rows)
	void normalize(int norm_type, float *wgts);
	void normalize(int norm_type, vector<float> &wgts) { return normalize(norm_type, &wgts[0]); }
	void normalize(int norm_type = Normalize_Cols) { return normalize(norm_type, NULL); }

	template <class S> void normalize(const vector<S>& external_avg, const vector<S>& external_std, int norm_type = 1);

	void get_cols_avg_std(vector<T>& _avg, vector<T>& _std);

	void print_row(FILE *fout, const string &prefix, const string &format, int i_row);

	void set_signals(vector<string> & sigs);

	//return true iff the matrix contains only valid floating point vals (not nan/infinite)
	//if the type of the matrix is not floating point, always returns true
	//if output = true, output the first invalid entry encountered to cerr
	bool is_valid(bool output = false) {
		if (std::is_floating_point<T>::value == false)
			return true;

		for (size_t i = 0; i < nrows; i++) {
			for (size_t j = 0; j < ncols; j++) {
				double x = (double)(m[i*ncols + j]);
				if (!isfinite(x)) {
					if (output)
						cerr << "invalid element(" << i << ", " << j << ") = " << x << endl;

					return false;
				}
			}
		}

		return true;
	}

	string my_class_name() const { return "MedMat"; }
};

// a few related util functions
void flags_to_indexes(vector<int> &flags, vector<int> &inds);

#include "MedMat_imp.h"

// a few basic tools for MedMat<float> mats 
int get_rand_medmat(MedMat<float> &A); // fills mat with uniform 0-1 numbers
int multiply_medmat(MedMat<float> &A, MedMat<float> &B, MedMat<float> &C); // A:n x m B:m x k --> gets C=A*B C:n x k
int fast_multiply_medmat_(const MedMat<float> &A, const MedMat<float> &B, MedMat<float> &C); // A:n x m B:m x k --> gets C=A*B C:n x k :: Uses Eigen to get performance
int fast_multiply_medmat(const MedMat<float> &A, const MedMat<float> &B, MedMat<float> &C, float s); // A:n x m B:m x k --> gets C=s*A*B C:n x k :: Uses Eigen to get performance
int fast_sum_medmat_rows(MedMat<float> &A, MedMat<float> &Asum, float factor); // A: n x m , output: Asum: 1 x m , summing all rows, done with matrix mult with factor (more efficient this way)
int fast_sum_medmat_cols(MedMat<float> &A, MedMat<float> &Asum, float factor); // A: n x m , output: Asum: n x 1 , summing all cols, done with matrix mult with factor (more efficient this way)
int fast_multiply_medmat_transpose(const MedMat<float> &A, const MedMat<float> &B, MedMat<float> &C, int transpose_flag); // A:n x m B:m x k --> gets C=A*B C:n x k , but allows transposing each mat
int fast_multiply_medmat_transpose(const MedMat<float> &A, const MedMat<float> &B, MedMat<float> &C, int transpose_flag, float s); // A:n x m B:m x k --> gets C=s*A*B C:n x k , but allows transposing each mat

int fast_multiply_scalar_vector(vector<float> &v, float s); //v = s * v; 
int fast_multiply_scalar_vector(vector<float> &v, float s, vector<float> &w); //w = s * v; 
int fast_element_dot_vector_vector(vector<float> &v, vector<float> &u, vector<float> &w); //w = v * u elementwise
int fast_element_dot_vector_vector(vector<float> &v, vector<float> &u); //v = v * u elementwise
int fast_element_affine_scalar(vector<float> &v, float s, vector<float> &u); // v = v + s*u element wise
int fast_element_affine_scalar(float s1, vector<float> &v, float s2, vector<float> &u); // v = s1*v + s2*u element wise
int fast_element_affine_scalar(vector<float> &v, float s, vector<float> &u, vector<float> &w); // w = v + s*u element wise

// gets inds[j] = i iff flags[i] != 0 


// pearson correlation between two columns (A,B can be the same mat)
//template <class T> double corr_mats_cols(MedMat<T> &A, int Acol, MedMat &B, int Bcol);

// next are useful to split matrices randomly to train/test
//void get_rand_binary_vec(vector<int> &v, double p, int len);
//template <class T> void split_mat_by_rows(MedMat<T> &A, vector<int> &flag, MedMat<T> &B, MedMat<T> &C);
//template <class T> void get_mat_by_rows(MedMat<T> &A, vector<int> &flag, MedMat<T> &B);
//void split_vector_by_flag(vector<int> &A, vector<int>flag, vector<int> &B, vector<int> &C);

// summasions
//template <class T> template <class S> void medmat_sum_rows(MedMat<T> &A, MedMat<S> &B);
//template <class T> template <class S> void medmat_sum_cols(MedMat<T> &A, MedMat<S> &B);
//template <class T> template <class S> void medmat_avg_rows(MedMat<T> &A, MedMat<S> &B);
//template <class T> template <class S> void medmat_avg_cols(MedMat<T> &A, MedMat<S> &B);
//template <class T> template <class S> void medmat_scalar_mult(MedMat<T> &A, S &s);


//========================================================
// Joining the MedSerialize Wagon
//========================================================
MEDSERIALIZE_SUPPORT(RecordData);
MEDSERIALIZE_SUPPORT(MedMat<float>);
MEDSERIALIZE_SUPPORT(MedMat<int>);
MEDSERIALIZE_SUPPORT(MedMat<double>);

#endif
