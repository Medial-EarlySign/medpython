// An assortment of comparison functions.

#ifndef __MED_COMPARE_H_
#define __MED_COMPARE_H_


#include "stdlib.h"
#include "stdio.h"
#include "string.h"

#include <vector> 


using namespace std ;

struct compare_indexed_float {
    bool operator()(const pair<int, float>& left, const pair<int, float>& right) {
		return (left.second < right.second) ;
    }
} ;

#endif