#include "MedModel.h"
#include "MedProcessUtils.h"
#include <omp.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/find.hpp>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string/regex.hpp>
#include <string>
#include "StripComments.h"

#define LOCAL_SECTION LOG_MED_MODEL
#define LOCAL_LEVEL	LOG_DEF_LEVEL

#define CHECK_CRC 0

using namespace boost::property_tree;

template <typename T>
std::vector<T> as_vector(ptree const& pt, ptree::key_type const& key)
{
	std::vector<T> r;
	for (auto& item : pt.get_child(key))
		r.push_back(item.second.get_value<T>());
	return r;
}

void MedModel::init_from_json_file_with_alterations(const string &fname, vector<string>& alterations) {
	string json_contents = file_to_string(0, fname, alterations);
	istringstream no_comments_stream(json_contents);

	MLOG("init model from json file [%s], stripping comments and displaying first 5 lines:\n", fname.c_str());
	int i = 5; string my_line;
	while (i-- > 0 && getline(no_comments_stream, my_line))
		MLOG("%s\n", my_line.c_str());
	no_comments_stream.clear();
	no_comments_stream.seekg(0);

	ptree pt;
	read_json(no_comments_stream, pt);
	string model_json_version = pt.get<string>("model_json_version");
	string ser = pt.get<string>("serialize_learning_set", to_string(this->serialize_learning_set).c_str());
	this->serialize_learning_set = stoi(ser);

	for (ptree::value_type &p : pt.get_child("model_actions"))
	{
		MLOG("[%s]\n", p.first);
		vector<vector<string>> all_attr_values;
		for (ptree::value_type &attr : p.second) {
			MLOG("\t[%s]\n", p.first);
		}
	}
	if (pt.count("predictor") > 0) {
		auto my_pred = pt.get_child("predictor");
		auto my_pred_params = pt.get_child("predictor_params");
		set_predictor(my_pred.data(), my_pred_params.data());
	}
	else MWARN("NOTE: no [predictor] node found in file\n");

}
