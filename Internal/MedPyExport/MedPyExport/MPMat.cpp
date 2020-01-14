#include "MPMat.h"
#include "MedMat/MedMat/MedMat.h"

const int MPMat::Normalize_Cols = MedMat<float>::Normalize_Cols;
const int MPMat::Normalize_Rows = MedMat<float>::Normalize_Rows;
const int MPMat::MISSING_VALUE = MED_MAT_MISSING_VALUE;

MPMat::MPMat() { o = new MedMat<float>(); }
MPMat::MPMat(int n_rows, int n_cols) { o = new MedMat<float>(n_rows, n_cols); };
MPMat::MPMat(float *IN_ARRAY2, unsigned long long DIM1, unsigned long long DIM2) { o = new MedMat<float>(IN_ARRAY2, DIM1, DIM2); };
MPMat::MPMat(MedMat<float>* from_ptr) { o_owned = false; o = from_ptr; }
MPMat::MPMat(const MPMat& other)
{
	o_owned = other.o_owned;
	if (!other.o_owned) {
		o = other.o;
	}
	else {
		o = new MedMat<float>(*(other.o));
	}
};

MPMat::~MPMat() { if (o_owned) delete o; }

void MPMat::get_numpy_copy(float** ARGOUTVIEWM_ARRAY2, unsigned long long* DIM1, unsigned long long* DIM2) {
	auto sz = o->size();
	*ARGOUTVIEWM_ARRAY2 = (float*)malloc(sz*sizeof(float));
	if (*ARGOUTVIEWM_ARRAY2 == nullptr)
		throw runtime_error("Out of memory while creating a copy of MedMat");
	memcpy(*ARGOUTVIEWM_ARRAY2, o->get_vec().data(), sz*sizeof(float));
	*DIM1 = o->nrows;
	*DIM2 = o->ncols;
}

void MPMat::get_numpy_view_unsafe(float** ARGOUTVIEW_ARRAY2, unsigned long long* DIM1, unsigned long long* DIM2) {
	*ARGOUTVIEW_ARRAY2 = o->get_vec().data();
	*DIM1 = o->nrows;
	*DIM2 = o->ncols;
}

void MPMat::load_numpy(float* IN_ARRAY2, unsigned long long DIM1, unsigned long long DIM2) {
	o->load(IN_ARRAY2, DIM1, DIM2);
}


void MPMat::test(float** ARGOUTVIEWM_ARRAY2, unsigned long long* DIM1, unsigned long long* DIM2,int n_rows, int n_cols) {
	MedMat<float> ret(n_rows, n_cols);
	auto sz = n_rows*n_cols;
	*ARGOUTVIEWM_ARRAY2 = (float*)malloc(sz*sizeof(float));
	if (*ARGOUTVIEWM_ARRAY2 == nullptr)
		throw runtime_error("Out of memory while creating a copy of MedMat");
	for (int row = 0; row < n_rows;++row)
		for (int col = 0; col < n_cols; ++col)
		{
			(*ARGOUTVIEWM_ARRAY2)[((unsigned long long)row)*n_cols + col] = col*1000+row;
		}
	*DIM1 = n_rows;
	*DIM2 = n_cols;

	o->load(*ARGOUTVIEWM_ARRAY2, n_rows, n_cols);
}

int MPMat::MEDPY_GET_Normalize_Cols() { return MedMat<float>::Normalize_Cols; }
int MPMat::MEDPY_GET_Normalize_Rows() { return MedMat<float>::Normalize_Rows; }

int MPMat::MEDPY_GET_nrows() { return o->nrows; };
int MPMat::MEDPY_GET_ncols() { return o->ncols; };

unsigned long long MPMat::__len__() { return o->size(); };

float MPMat::__getitem__(const vector<int>& index) { return o->get(index[0], index[1]); }
void MPMat::__setitem__(const vector<int>& index, float val) { o->set(index[0], index[1]) = val; }

void MPMat::MEDPY_GET_row_ids(MEDPY_NP_OUTPUT(int** out_row_ids_buf, unsigned long long* out_row_ids_buf_len)) {
	vector_to_buf(o->row_ids, out_row_ids_buf, out_row_ids_buf_len);
}
void MPMat::MEDPY_SET_row_ids(MEDPY_NP_INPUT(int* row_ids_buf, unsigned long long row_ids_buf_len)) {
	buf_to_vector(row_ids_buf, row_ids_buf_len, o->row_ids);
}

void MPMat::MEDPY_GET_avg(MEDPY_NP_OUTPUT(int** avg_buf, unsigned long long* avg_buf_len)) {
	vector_to_buf(o->avg, avg_buf, avg_buf_len);
}

void MPMat::MEDPY_GET_std(MEDPY_NP_OUTPUT(int** std_buf, unsigned long long* std_buf_len)) {
	vector_to_buf(o->std, std_buf, std_buf_len);
}

void MPMat::set_signals(std::vector<std::string> & sigs) { o->set_signals(sigs); };

std::vector<std::string> MPMat::MEDPY_GET_signals() { return o->signals; }

int MPMat::MEDPY_GET_normalized_flag() { return o->normalized_flag; };
void MPMat::MEDPY_SET_normalized_flag(int newval) { o->normalized_flag = newval; };

int MPMat::MEDPY_GET_transposed_flag() { return o->transposed_flag; };
void MPMat::MEDPY_SET_transposed_flag(int newval) { o->transposed_flag = newval; };

float MPMat::MEDPY_GET_missing_value() { return o->missing_value; };
void MPMat::MEDPY_SET_missing_value(float newval) { o->missing_value = newval; };

void MPMat::zero() { o->zero(); };
void MPMat::set_val(float val) { o->set_val(val); };
void MPMat::clear() { o->clear(); };

void MPMat::resize(int n_rows, int n_cols) { o->resize(n_rows, n_cols); };

int MPMat::read_from_bin_file(const string &fname) { return o->read_from_bin_file(fname); };
int MPMat::write_to_bin_file(const string &fname) { return o->write_to_bin_file(fname); };
int MPMat::write_to_csv_file(const string &fname) { return o->write_to_csv_file(fname); };
int MPMat::read_from_csv_file(const string &fname, int titles_line_flag) { return o->read_from_csv_file(fname, titles_line_flag); };

void MPMat::transpose() { o->transpose(); };

void MPMat::get_sub_mat(vector<int> &rows_to_take, vector<int> &cols_to_take) { o->get_sub_mat(rows_to_take,cols_to_take); }
void MPMat::get_sub_mat_by_flags(vector<int> &rows_to_take_flag, vector<int> &cols_to_take_flag) { o->get_sub_mat_by_flags(rows_to_take_flag, cols_to_take_flag); };
//void MPMat::reorder_by_row(vector<int> &row_order) { o->reorder_by_row(row_order); };
//void MPMat::reorder_by_col(vector<int> &col_order) { o->reorder_by_col(col_order); };

void MPMat::add_rows(MPMat& m_add) { o->add_rows(*(m_add.o)); };
void MPMat::add_rows(MEDPY_NP_INPUT(float* m_add, unsigned long long nrows_to_add)) { o->add_rows(m_add, nrows_to_add); };
void MPMat::add_cols(MPMat& m_add) { o->add_cols(*(m_add.o)); };
void MPMat::add_cols(MEDPY_NP_INPUT(float* m_add, unsigned long long ncols_to_add)) { o->add_cols(m_add, ncols_to_add); };

void MPMat::get_row(int i_row, MEDPY_NP_OUTPUT(float** rowv, unsigned long long* rowv_n)) const {
	vector<float> ret;
	o->get_row(i_row, ret);
	vector_to_buf(ret, rowv, rowv_n);
}
void MPMat::get_col(int i_col, MEDPY_NP_OUTPUT(float** colv, unsigned long long* colv_n)) const {
	vector<float> ret;
	o->get_col(i_col, ret);
	vector_to_buf(ret, colv, colv_n);
}

void MPMat::normalize(int norm_type, MEDPY_NP_INPUT(float* wgts, unsigned long long wgts_n)) {
	if (norm_type == Normalize_Cols) {
		if (wgts_n < o->nrows)
			throw runtime_error("normalize failed: wgts should be longer than nrows");
	}else{
		if (wgts_n < o->ncols)
			throw runtime_error("normalize failed: wgts should be longer than ncols");
	}
	o->normalize(norm_type, wgts);
}
void MPMat::normalize(int norm_type) { o->normalize(norm_type); }

void MPMat::normalize(MEDPY_NP_INPUT(float* external_avg, unsigned long long external_avg_n),
	MEDPY_NP_INPUT(float* external_std, unsigned long long external_std_n), int norm_type) {
	vector<float> external_std_vec, external_avg_vec;
	buf_to_vector(external_avg, external_avg_n, external_avg_vec);
	buf_to_vector(external_std, external_std_n, external_std_vec);
}

void MPMat::get_cols_avg_std(MEDPY_NP_OUTPUT(float** buf_avg, unsigned long long* buf_avg_n), MEDPY_NP_OUTPUT(float** buf_std, unsigned long long* buf_std_n)) {
	vector<float> __avg, __std;
	o->get_cols_avg_std(__avg, __std);
	vector_to_buf(__avg, buf_avg, buf_avg_n);
	vector_to_buf(__std, buf_std, buf_std_n);
}

void MPMat::load(float *IN_ARRAY2, unsigned long long DIM1, unsigned long long DIM2) { o->load(IN_ARRAY2, DIM1, DIM2); };
void MPMat::load_transposed(float *IN_ARRAY2, unsigned long long DIM1, unsigned long long DIM2) { o->load_transposed(IN_ARRAY2, DIM1, DIM2); };
void MPMat::load(MPMat& x) { o->load(*(x.o)); };


bool MPMat::is_valid(bool output) { return o->is_valid(output); }

MPSerializableObject MPMat::asSerializable() { return MPSerializableObject(o); }
