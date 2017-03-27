#include "DoCalcFeatGenerator.h"

void DoCalcFeatGenerator::set_names() {
	if (names.empty()) {
		string name = "FTR_" + int_to_string_digits(serial_id, 6) + "." + signalName + ".";
		names.push_back(name);
	}
}

DEF_FEATURE_GENERATOR(DoCalcFeatGenerator);