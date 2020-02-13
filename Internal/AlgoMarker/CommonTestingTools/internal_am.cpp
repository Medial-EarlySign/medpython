#include "internal_am.h"
#include <AlgoMarker/AlgoMarker/AlgoMarker.h>

void set_am_matrix(void *am, const std::string& am_csv_file) {
	MedialInfraAlgoMarker *m_am = (MedialInfraAlgoMarker *)am;
	m_am->set_am_matrix(am_csv_file);
}


