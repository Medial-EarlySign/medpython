#include "MedDataStructures.h"

#define LOCAL_SECTION LOG_MEDMAT
#define LOCAL_LEVEL	LOG_DEF_LEVEL
extern MedLogger global_logger;

//-----------------------------------------------------------------------------
void MedHeap::push(float score, int index)
{
	heap_elem he;

	he.score = score;
	he.index = index;

	h.push_back(he);
	push_heap(h.begin(),h.end());

	while (size() > max_size) {
		pop_heap(h.begin(),h.end());
		h.pop_back();
	}
}
