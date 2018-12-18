//
// MedDataStructures.h :: several needed data structures and handling them
//

#ifndef __MED_DATA_STRUCTURES_H__
#define __MED_DATA_STRUCTURES_H__

#include "Logger/Logger/Logger.h"
#include <vector>
#include <algorithm>

using namespace std;

// ================================================================
//   HEAP
// ================================================================
struct heap_elem {
	float score;
	int index;
};
inline bool operator< (const heap_elem &l, const heap_elem &r) {return (l.score < r.score);}

// simple heap to efficiently get lowest N scores and their indexes
// if max needed - simply use -score
class MedHeap {
	public:
		vector<heap_elem> h;
		int max_size;
		int size() {return (int)h.size();};
		int set_max_size(int s) {max_size = s;}
		
		MedHeap(int s) {max_size = s; h.clear(); make_heap(h.begin(),h.end());}
		void clear(){h.clear(); make_heap(h.begin(),h.end());}

		void push(float score, int index);
};




#endif
