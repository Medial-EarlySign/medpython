#include "FeatureSelector.h"
#include "Logger/Logger/Logger.h"

#define LOCAL_SECTION LOG_FEAT_SELECTOR
#define LOCAL_LEVEL	LOG_DEF_LEVEL

//=======================================================================================
// FeatsCleaner
//=======================================================================================
// Cleaner types
FeatureSelectorTypes feature_selector_name_to_type(const string& selector_name) {

	if (selector_name == "all")
		return FTR_SLCTR_ALL;
	else
		return FTR_SLCTR_LAST;
}

// Initialization
//.......................................................................................
FeatureSelector* FeatureSelector::make_selector(string selector_name) {

	return make_selector(feature_selector_name_to_type(selector_name));
}

//.......................................................................................
FeatureSelector * FeatureSelector::make_selector(string selector_name, string init_string) {

	return make_selector(feature_selector_name_to_type(selector_name), init_string);
}

//.......................................................................................
FeatureSelector * FeatureSelector::make_selector(FeatureSelectorTypes selector_type) {

	if (selector_type == FTR_SLCTR_ALL)
		return new DummyFeatsSelector;
	else
		return NULL;

}

//.......................................................................................
FeatureSelector * FeatureSelector::make_selector(FeatureSelectorTypes selector_type, string init_string) {

	FeatureSelector *newSelector = make_selector(selector_type);
	newSelector->init_from_string(init_string);
	return newSelector;
}

// (De)Serialize
//.......................................................................................
size_t FeatureSelector::get_selector_size() {
	return sizeof(selector_type) + get_size();
}

//.......................................................................................
size_t FeatureSelector::selector_serialize(unsigned char *blob) {

	size_t ptr = 0;
	memcpy(blob + ptr, &selector_type, sizeof(FeatureSelectorTypes)); ptr += sizeof(FeatureSelectorTypes);
	ptr += serialize(blob + ptr);

	return ptr;
}
