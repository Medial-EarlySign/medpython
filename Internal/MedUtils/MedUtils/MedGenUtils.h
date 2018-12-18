//
// MedGenUtils - general needed utilities (time conversions, vector ops, randomizing, etc)
//

#ifndef __MED_GEN_UTILS_H__
#define __MED_GEN_UTILS_H__

#include <stdlib.h>
#include <vector>
#include <map>
#include <MedUtils/MedUtils/MedGlobalRNG.h>


#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <Logger/Logger/Logger.h>

using namespace std;

//================================================================
//  Time / Date related
//================================================================
// days from 01011900 - function is in h file since it is inline.
//inline 
inline int date_to_days(int date);
inline int year_to_days(int year) {return date_to_days(year*10000+0701);}
inline float get_age(int date, int byear) {return ((float)(date_to_days(date) - year_to_days(byear))/365.0f);}

//================================================================
// vector/array operations
//================================================================

#define VEC_DATA(v) ((v).size()>0 ? &((v)[0]) : NULL)
//#define VEC_DATA(v) ((v).size()>0 ? (v).data() : 0)

// get the a vector of all indexes such that v[i] == 0
template <class T> void get_zero_inds(T *v, int len, vector<int> &inds);
// get the a vector of all indexes such that v[i] != 0
template <class T> void get_nonzero_inds(T *v, int len, vector<int> &inds);
// get the a vector of all indexes such that v[i] == 0 , c++ ver
template <class T> void get_zero_inds(vector<T> &v, vector<int> &inds) {get_zero_inds(&v[0],(int)v.size(),inds);}
// get the a vector of all indexes such that v[i] != 0 , c++ ver
template <class T> void get_nonzero_inds(vector<T> &v, vector<int> &inds) {get_nonzero_inds(&v[0],(int)v.size(),inds);}
// gets all elements in a vector<vector<T>> to a single vector<T>
template <class T> void get_vec_from_vecvec(vector<vector<T>> &v_in, vector<T> &v_out);
// gets number of different values in a vector
template <class T> int get_vec_ndiff_vals(vector<T> &v);

// generate an arithmetic sequence start:step:finish. if step <= 0 returns -1 otherwise 0,
// and fills seq with the sequence
// isForward - go from start to finish or backward from finish
template<typename T> int sequence(T start, T finish, T step, vector<T>& seq, bool isForward = true);

// gets a (sorted) partition on the values in a vector and gets the matching categorized vec
void categorize_vec(vector<float> &in, vector<float> &bounds, vector<float> &out);

void get_probs_vec(vector<float> &v);

//================================================================
// Random Numbers
//================================================================

// wrappers for globalRNG::rand()
//....................

// set the random seed, if called with -1 will take seed from current time
void set_rand_seed(int seed);

void get_rand_splits(vector<int> &split, int nsplits, int size);

// gets first size elements of a permutaion of 0..N-1
void get_rand_vector_no_repetitions(vector<int> &v, int N, int size);

// gets size random values in range 0...N-1
void get_rand_vector_with_repetitions(vector<int> &v, int N, int size);


//==========================================================
// Inline functions
//==========================================================
// conversions
// days from 01011900 - function is in h file since it is inline.
inline int date_to_days(int date)
{
	int y,d,m;

	y = date/10000;
	m = (date - y*10000)/100;
	d = date - y*10000 - m*100;

	return (y*365 + (m-1)*30 + d);
	
} 

int get_day(int val) ;
int get_day_approximate(int val);
int get_date(int days);

//================================================================
// Initialization Utilities
//================================================================
// OLD CODE
//int initialization_text_to_map(const string& text, map<string, string>& init_map);

//===========================
// Opertaing System Utilities
//===========================

//is the OS Windows
bool is_windows_os(void);

// templated functions in _imp file
#include "MedGenUtils_imp.h"

#endif
