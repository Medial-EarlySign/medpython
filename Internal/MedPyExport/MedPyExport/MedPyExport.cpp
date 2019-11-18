#include "MedPyExport.h"
static const std::vector<std::string> public_objects = { PUBLIC_OBJECTS };

std::vector<std::string> get_public_objects() {
	return public_objects;
}