#include "MedMat.h"
#include <External/Eigen/Core>
#include <MedTime/MedTime/MedTime.h>
#include <MedUtils/MedUtils/MedGlobalRNG.h>
#include <Logger/Logger/Logger.h>
#include <fstream>
#include <boost/algorithm/string.hpp>

#define LOCAL_SECTION LOG_MEDMAT
#define LOCAL_LEVEL	LOG_DEF_LEVEL
extern MedLogger global_logger;

using namespace Eigen;
#include <thread>

//...........................................................................................
void flags_to_indexes(vector<int> &flags, vector<int> &inds)
{
	inds.clear();
	if (flags.size() == 0)
		return;

	for (int i=0; i<flags.size(); i++) {
		if (flags[i] != 0)
			inds.push_back(i);
	}
}

//...........................................................................................
int get_rand_medmat(MedMat<float> &A) // fills mat with uniform 0-1 numbers
{
	vector<float> &m = A.get_vec();
	for (size_t i=0; i<A.size(); i++) {
		m[i] = rand_1();
	}
	return 0;
}

int fast_multiply_scalar_vector(vector<float> &v, float s) //v = s * v
{
	Map<MatrixXf> x(&v[0], 1, v.size());

	x *= s;

	return 0;
}

int fast_multiply_scalar_vector(vector<float> &v, float s, vector<float> &w) //w = s * v
{
	Map<MatrixXf> x(&v[0], 1, v.size());
	Map<MatrixXf> y(&w[0], 1, w.size());

	y = x * s;

	return 0;
}

int fast_element_dot_vector_vector(vector<float> &v, vector<float> &u, vector<float> &w) //w = v * u elementwise
{
	Map<MatrixXf> x(&v[0], 1, v.size());
	Map<MatrixXf> y(&u[0], 1, u.size());
	Map<MatrixXf> z(&w[0], 1, w.size());

	z = x.array() * y.array();

	return 0;
}

int fast_element_dot_vector_vector(vector<float> &v, vector<float> &u) //v = v * u elementwise
{
	Map<MatrixXf> x(&v[0], 1, v.size());
	Map<MatrixXf> y(&u[0], 1, u.size());

	x = x.array() * y.array();

	return 0;
}

int fast_element_affine_scalar(vector<float> &v, float s, vector<float> &u) // v = v + s*u element wise
{
	Map<MatrixXf> x(&v[0], 1, v.size());
	Map<MatrixXf> y(&u[0], 1, u.size());

	x = x.array() + s*y.array();

	return 0;
}

int fast_element_affine_scalar(float s1, vector<float> &v, float s2, vector<float> &u) // v = s1*v + s2*u element wise
{
	Map<MatrixXf> x(&v[0], 1, v.size());
	Map<MatrixXf> y(&u[0], 1, u.size());

	x = s1*x.array() + s2*y.array();

	return 0;
}

int fast_element_affine_scalar(vector<float> &v, float s, vector<float> &u, vector<float> &w) // w = v + s*u element wise
{
	Map<MatrixXf> x(&v[0], 1, v.size());
	Map<MatrixXf> y(&u[0], 1, u.size());
	Map<MatrixXf> z(&w[0], 1, w.size());

	z = x.array() +  s*y.array();

	return 0;
}

int fast_multiply_medmat_(const MedMat<float> &A, const MedMat<float> &B, MedMat<float> &C) // A:n x m B:m x k --> gets C=s*A*B C:n x k
{
	if (A.ncols != B.nrows) {
		MERR("ERROR: multiply_medmat: Mats dimension don't match: (%d x %d) , (%d x %d) ... \n", A.nrows, A.ncols, B.nrows, B.ncols);
		return -1;
	}

	//int ncores = std::thread::hardware_concurrency();
	//Eigen::setNbThreads(3*ncores/4);

	C.resize(A.nrows, B.ncols);

	Map<const MatrixXf> x(A.data_ptr(), A.ncols, A.nrows);
	Map<const MatrixXf> y(B.data_ptr(), B.ncols, B.nrows);

	Map<MatrixXf> z(C.data_ptr(), C.ncols, C.nrows);

	z = y*x;

	return 0;

}


//...........................................................................................
// multiplying using the (really fast) Eigen library.
int fast_multiply_medmat(const MedMat<float> &A, const MedMat<float> &B, MedMat<float> &C, float s) // A:n x m B:m x k --> gets C=s*A*B C:n x k
{
	if (A.ncols != B.nrows) {
		MERR("ERROR: multiply_medmat: Mats dimension don't match: (%d x %d) , (%d x %d) ... \n", A.nrows, A.ncols, B.nrows, B.ncols);
		return -1;
	}

	//int ncores = std::thread::hardware_concurrency();
	//Eigen::setNbThreads(3*ncores/4);

	C.resize(A.nrows, B.ncols);

	Map<const MatrixXf> x(A.data_ptr(), A.ncols, A.nrows);
	Map<const MatrixXf> y(B.data_ptr(), B.ncols, B.nrows);

	Map<MatrixXf> z(C.data_ptr(), C.ncols, C.nrows);

	z = s*y*x;

	return 0;

}

//...........................................................................................
// multiplying using the (really fast) Eigen library.
int fast_multiply_medmat_transpose(const MedMat<float> &A, const MedMat<float> &B, MedMat<float> &C, int transpose_flag) // A:n x m B:m x k --> gets C=A*B C:n x k , but allows transposing each mat
{
	if ((transpose_flag == 0x0 && A.ncols != B.nrows) ||
		(transpose_flag == 0x1 && A.nrows != B.nrows) ||
		(transpose_flag == 0x2 && A.ncols != B.ncols) ||
		(transpose_flag == 0x3 && A.nrows != B.ncols))
	{
		MERR("ERROR: multiply_medmat: Mats dimension don't match: (%d x %d) , (%d x %d) transpose_flag %d ... \n",A.nrows,A.ncols,B.nrows,B.ncols,transpose_flag);
		return -1;
	}

	//int ncores = std::thread::hardware_concurrency();
	//Eigen::setNbThreads(3*ncores/4);

	Map<const MatrixXf> x(A.data_ptr(),A.ncols,A.nrows);
	Map<const MatrixXf> y(B.data_ptr(),B.ncols,B.nrows);

	int nr=0, nc=0;
	if (transpose_flag == 0x0) {nr=A.nrows; nc=B.ncols;}
	else if (transpose_flag == 0x1) {nr=A.ncols; nc=B.ncols;}
	else if (transpose_flag == 0x2) {nr=A.nrows; nc=B.nrows;}
	else if (transpose_flag == 0x3) {nr=A.ncols; nc=B.nrows;}

	C.resize(nr,nc);
	Map<MatrixXf> z(C.data_ptr(),C.ncols,C.nrows);

	if (transpose_flag == 0x0) {z = y*x;}
	else if (transpose_flag == 0x1) {z = y*x.transpose();}
	else if (transpose_flag == 0x2) {z = y.transpose()*x;}
	else if (transpose_flag == 0x3) {z = y.transpose()*x.transpose();}

	return 0;

}


int fast_multiply_medmat_transpose(const MedMat<float> &A, const MedMat<float> &B, MedMat<float> &C, int transpose_flag, float s) // A:n x m B:m x k --> gets C=A*B C:n x k , but allows transposing each mat
{
	if ((transpose_flag == 0x0 && A.ncols != B.nrows) ||
		(transpose_flag == 0x1 && A.nrows != B.nrows) ||
		(transpose_flag == 0x2 && A.ncols != B.ncols) ||
		(transpose_flag == 0x3 && A.nrows != B.ncols))
	{
		MERR("ERROR: multiply_medmat: Mats dimension don't match: (%d x %d) , (%d x %d) transpose_flag %d ... \n", A.nrows, A.ncols, B.nrows, B.ncols, transpose_flag);
		return -1;
	}

	//int ncores = std::thread::hardware_concurrency();
	//Eigen::setNbThreads(3*ncores/4);

	Map<const MatrixXf> x(A.data_ptr(), A.ncols, A.nrows);
	Map<const MatrixXf> y(B.data_ptr(), B.ncols, B.nrows);

	int nr = 0, nc = 0;
	if (transpose_flag == 0x0) { nr = A.nrows; nc = B.ncols; }
	else if (transpose_flag == 0x1) { nr = A.ncols; nc = B.ncols; }
	else if (transpose_flag == 0x2) { nr = A.nrows; nc = B.nrows; }
	else if (transpose_flag == 0x3) { nr = A.ncols; nc = B.nrows; }

	C.resize(nr, nc);
	Map<MatrixXf> z(C.data_ptr(), C.ncols, C.nrows);

	if (transpose_flag == 0x0) { z = s * y * x; }
	else if (transpose_flag == 0x1) { z = s * y * x.transpose(); }
	else if (transpose_flag == 0x2) { z = s * y.transpose()*x; }
	else if (transpose_flag == 0x3) { z = s * y.transpose()*x.transpose(); }

	return 0;

}


//...........................................................................................
int multiply_medmat(MedMat<float> &A, MedMat<float> &B, MedMat<float> &C) // A:n x m B:m x k --> gets C=A*B C:n x k
{
	if (A.ncols != B.nrows) {
		MERR("ERROR: multiply_medmat: Mats dimension don't match: (%d x %d) , (%d x %d) ... \n",A.nrows,A.ncols,B.nrows,B.ncols);
		return -1;
	}

	C.resize(A.nrows,B.ncols);

	C.zero();

	int i,j,k;

	for (i=0; i<A.nrows; i++)
		for (j=0; j<B.ncols; j++)
			for (k=0; k<A.ncols; k++)
				C(i,j) += A(i,k)*B(k,j);

	return 0;

}

//...........................................................................................
// A: n x m , output: Asum: 1 x m , summing all rows, done with matrix mult with factor (more efficient this way)
int fast_sum_medmat_rows(MedMat<float> &A, MedMat<float> &Asum, float factor)
{
	if (A.ncols == 0 || A.nrows == 0)
		return -1;

	MedMat<float> Ones(1, A.nrows);
	vector<float> &m = Ones.get_vec();
	fill(m.begin(), m.end(), factor);

	return (fast_multiply_medmat_(Ones, A, Asum));
}

//...........................................................................................
// A: n x m , output: Asum: n x 1 , summing all cols, done with matrix mult with factor (more efficient this way)
int fast_sum_medmat_cols(MedMat<float> &A, MedMat<float> &Asum, float factor)
{
	if (A.ncols == 0 || A.nrows == 0)
		return -1;

	MedMat<float> Ones(A.ncols, 1);
	vector<float> &m = Ones.get_vec();
	fill(m.begin(), m.end(), factor);

	return (fast_multiply_medmat_(A, Ones, Asum));
}

//...........................................................................................
//double corr_mats_cols(MedMat &A, int Acol, MedMat &B, int Bcol)
//{
//	double sx,sy,sxy,sxx,syy,n;
//
//	if (A.nrows != B.nrows)
//		return -2;
//	if (A.nrows == 0)
//		return -2;
//
//	sx = sy = sxy = sxx = syy = 0;
//
//	for (int i=0; i<A.nrows; i++) {
//		sx += A(i,Acol);
//		sy += B(i,Bcol);
//		sxx += A(i,Acol)*A(i,Acol);
//		syy += B(i,Bcol)*B(i,Bcol);
//		sxy += A(i,Acol)*B(i,Bcol);
//	}
//
//	n = (double)A.nrows;
//
//	sx /= n;
//	sy /= n;
//	sxx /= n;
//	syy /= n;
//	sxy /= n;
//
//	double c1 = sxy - sx*sy;
//	double c2 = sxx - sx*sx;
//	double c3 = syy - sy*sy;
//
//	double epsilon = 1e-8;
//	if (c2 < epsilon || c3 < epsilon)
//		return 0;
//	return (c1/(sqrt(c2)*sqrt(c3)));
//}

////...........................................................................................
//void get_rand_binary_vec(vector<int> &v, double p, int len)
//{
//	v.resize(len);
//	for (int i=0; i<len; i++) {
//		double r = (double)globalRNG::rand()/(double)globalRNG::max();
//		if (r < p)
//			v[i] = 1;
//		else
//			v[i] = 0;
//	}
//}
//
////...........................................................................................
//void split_mat_by_rows(MedMat &A, vector<int> &flag, MedMat &B, MedMat &C)
//{
//	float *a = A.data_ptr();
//
//	B.resize(0,A.ncols);
//	C.resize(0,A.ncols);
//
//	for (int i=0; i<A.nrows; i++)
//		if (flag[i])
//			B.add_rows(&a[i*A.ncols],1);
//		else
//			C.add_rows(&a[i*A.ncols],1);
//}
//
////...........................................................................................
//void split_vector_by_flag(vector<int> &A, vector<int>flag, vector<int> &B, vector<int> &C)
//{
//	B.clear();
//	C.clear();
//
//	for (int i=0; i<A.size(); i++) {
//		if (flag[i])
//			B.push_back(A[i]);
//		else
//			C.push_back(A[i]);
//	}
//}
//
////...........................................................................................
//void get_mat_by_rows(MedMat &A, vector<int> &flag, MedMat &B)
//{
//	float *a = A.data_ptr();
//
//	B.resize(0,A.ncols);
//
//	for (int i=0; i<A.nrows; i++)
//		if (flag[i])
//			B.add_rows(&a[i*A.ncols],1);
//}
//
////...........................................................................................
//void medmat_sum_rows(MedMat &A, MedMat &B)
//{
//	B.resize(1,A.ncols);
//	B.zero();
//	for (int i=0; i<A.nrows; i++)
//		for (int j=0; j<A.ncols; j++)
//			B(0,j) += A(i,j);
//}
//
////...........................................................................................
//void medmat_sum_cols(MedMat &A, MedMat &B)
//{
//	B.resize(A.nrows,1);
//	B.zero();
//	for (int i=0; i<A.nrows; i++)
//		for (int j=0; j<A.ncols; j++)
//			B(i,0) += A(i,j);
//}
//
////...........................................................................................
//void medmat_avg_rows(MedMat &A, MedMat &B)
//{
//	medmat_sum_rows(A,B);
//	if (A.nrows > 0)
//		medmat_scalar_mult(B,(float)1/(float)A.nrows);
//}
//
////...........................................................................................
//void medmat_avg_cols(MedMat &A, MedMat &B)
//{
//	medmat_sum_cols(A,B);
//	if (A.ncols > 0)
//		medmat_scalar_mult(B,(float)1/(float)A.ncols);
//}
//
////...........................................................................................
//void medmat_scalar_mult(MedMat &A, float s)
//{
//	for (int i=0; i<A.nrows; i++)
//		for (int j=0; j<A.ncols; j++)
//			A(i,j) *= s;	
//}
