//
// MedDictionary.c
//
#define __INFRAMED_DLL

#ifndef _SCL_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS
#endif

#include "InfraMed.h"
#include "MedDictionary.h"
#include "Logger/Logger/Logger.h"
#include <fstream>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
//#define _SCL_SECURE_NO_WARNINGS

#include <queue>
#include <unordered_set>
#include <algorithm>

#define LOCAL_SECTION LOG_DICT
#define LOCAL_LEVEL LOG_DEF_LEVEL
extern MedLogger global_logger;

using namespace boost;
using namespace std;

mutex lock_dict_changes;

//-----------------------------------------------------------------------------------------------
int MedDictionary::read(vector<string> &dfnames)
{
	int rc = 0;
	for (int i = 0; i < dfnames.size(); i++) {
		rc += read(dfnames[i]);
	}
	return rc;
}

//-----------------------------------------------------------------------------------------------
int MedDictionary::read(string path, vector<string> &dfnames)
{
	int rc = 0;
	for (int i = 0; i < dfnames.size(); i++) {
		string fname = (path == "") ? dfnames[i] : path + "/" + dfnames[i];
		rc += read(fname);
	}
	return rc;
}

//-----------------------------------------------------------------------------------------------
// read:: gets a dictionary file name, opens it, parses it, and loads internal data structures
int MedDictionary::read(const string &fname)
{
	ifstream inf(fname);

	if (!inf) {
		MERR("MedDictionary: read: Can't open file %s\n", fname.c_str());
		return -1;
	}

	MLOG_D("MedDictinary: read: reading dictionary file %s\n", fname.c_str());
	fnames.push_back(fname); // TD : check that we didn't already load this file
	string curr_line;
	while (getline(inf, curr_line)) {
		trim(curr_line);
		if ((curr_line.size() > 1) && (curr_line[0] != '#')) {

			vector<string> fields;
			split(fields, curr_line, boost::is_any_of("\t"));

			if (fields.size() >= 3) {
				if (fields[0].compare(0, 3, "DEF") == 0) {
					int n_id = stoi(fields[1]);
					if (Name2Id.find(fields[2]) != Name2Id.end())
						MWARN("Warning::MedDictionary.read read id %d with existing key name \"%s\" ("
							"with maybe diffrent id) in dictionary %s\n", n_id, fields[2].c_str(),
							fname.c_str());
					Name2Id[fields[2]] = n_id;
					Id2Name[n_id] = fields[2];
					if (used.find(fields[2]) == used.end()) {
						Id2Names[n_id].push_back(fields[2]);
						used[fields[2]] = 1;
					}
					MLOG_D("n_id = %d Name =@@@ %s @@@, name2id %d , id2name %s\n", n_id, fields[2].c_str(), Name2Id[fields[2]], Id2Name[n_id].c_str());
				}
				else if (fields[0].compare(0, 6, "SIGNAL") == 0) {
					int n_id = stoi(fields[2]);
					if (Name2Id.find(fields[2]) != Name2Id.end())
						MWARN("Warning::Signal \"%s\" is defined again. now with id %d ("
							"with maybe diffrent id) in file %s\n", fields[2].c_str(), n_id,
							fname.c_str());
					Name2Id[fields[1]] = n_id;
					Id2Name[n_id] = fields[1];
					if (used.find(fields[1]) == used.end()) {
						Id2Names[n_id].push_back(fields[1]);
						used[fields[1]] = 1;
					}
					MLOG_D("SIG n_id = %d Name =@@@ %s @@@, name2id %d , id2name %s\n", n_id, fields[1].c_str(), Name2Id[fields[1]], Id2Name[n_id].c_str());
				}
				else if (fields[0].compare(0, 3, "SET") == 0) {
					if (Name2Id.find(fields[1]) == Name2Id.end() || Name2Id.find(fields[2]) == Name2Id.end()) {
						MERR("MedDictionary: read: SET line with undefined elements: [%s] [%s]\n", fields[1].c_str(), fields[2].c_str());
					}
					else {
						int member_n = Name2Id[fields[2]];
						int set_n = Name2Id[fields[1]];

						pair<int, int> p;
						p.first = set_n;
						p.second = member_n;

						int first_def = 1;
						if (MemberInSet.find(p) != MemberInSet.end()) first_def = 0;

						MemberInSet[p] = 1;

						MLOG_D("Added SET : %s : ---> : %s : %d ---> %d\n", fields[2].c_str(), fields[1].c_str(), member_n, set_n);

						if (first_def) {
							if (Set2Members.find(set_n) == Set2Members.end())
								Set2Members[set_n] = vector<int>();
							Set2Members[set_n].push_back(member_n);
							if (Member2Sets.find(member_n) == Member2Sets.end())
								Member2Sets[member_n] = vector<int>();
							Member2Sets[member_n].push_back(set_n);
						}
					}
				}

			}
			else if (fields[0] != "SECTION")
				MWARN("Warning in MedDictionary::read in file %s - got line with less than 3 tokens: \"%s\"\n", 
					fname.c_str(), curr_line.c_str());
		}

	}
	inf.close();
	return 0;
}

//-----------------------------------------------------------------------------------------------
int MedDictionary::id_list(vector<string> &names, vector<int> &ids)
{
	int rc = 0;

	ids.resize(names.size());
	for (int i = 0; i < names.size(); i++) {
		rc += (ids[i] = id(names[i]));
	}
	return rc;
}
//-----------------------------------------------------------------------------------------------
int MedDictionary::id(const string &name)
{
	if (Name2Id.find(name) == Name2Id.end())
		return -1;
	return Name2Id[name];
}

//-----------------------------------------------------------------------------------------------
string MedDictionary::name(int id)
{
	if (Id2Name.find(id) == Id2Name.end())
		return string("");
	return Id2Name[id];
}

//-----------------------------------------------------------------------------------------------
int MedDictionary::is_in_set(const string& member, int set_id)
{
	if (Name2Id.find(member) == Name2Id.end())
		return 0;

	return is_in_set(set_id, Name2Id[member]);
}

//-----------------------------------------------------------------------------------------------
int MedDictionary::is_in_set(int member_id, const string &set_name)
{
	if (Name2Id.find(set_name) == Name2Id.end())
		return 0;

	return is_in_set(member_id, Name2Id[set_name]);
}

//-----------------------------------------------------------------------------------------------
int MedDictionary::is_in_set(const string& member, const string& set_name)
{
	if (Name2Id.find(member) == Name2Id.end())
		return 0;
	if (Name2Id.find(set_name) == Name2Id.end())
		return 0;
	return is_in_set(Name2Id[member], Name2Id[set_name]);
}

//-----------------------------------------------------------------------------------------------
int MedDictionary::is_in_set(int member_id, int set_id)
{
	pair<int, int> p;
	p.first = set_id;
	p.second = member_id;

	// first (backwward compatible) testing a quick test for direct inclusion.
	if ((member_id == set_id) || (MemberInSet.find(p) != MemberInSet.end()))
		return 1;

	if (Member2Sets.find(member_id) == Member2Sets.end())
		return 0;
	// now we are at a case in which the member is not a direct member of the set, 
	// yet it may be possible that it is a memeber of a set included in the set (at any depth).
	// in order to test that we have to search all possible paths from our member

	map<int, int> used;
	queue<int> set_q;

	// init first level (we know there's no direct inclusion)
	for (int i = 0; i < Member2Sets[member_id].size(); i++) {
		set_q.push(Member2Sets[member_id][i]);
		used[Member2Sets[member_id][i]] = 1;
	}

	// pop , compare and search
	int *member2set = NULL;
	while (set_q.size() > 0) {

		int id_set = set_q.front();
		set_q.pop();

		if (Member2Sets.find(id_set) != Member2Sets.end()) {

			int len = (int)Member2Sets[id_set].size();
			if (len > 0)
				member2set = &Member2Sets[id_set][0];

			for (int i = 0; i < len; i++) {

				int new_set = member2set[i];
				if (used.find(new_set) == used.end()) {

					if (new_set == set_id)
						return 1; // we found it !!

					set_q.push(new_set);
					used[new_set] = 1;
				}

			}

		}

	}

	return 0;
}

//-----------------------------------------------------------------------------------------------------------
// next function calculates for each member ALL the sets it is contained in , applying all set transitivity
void MedDictionary::get_members_to_all_sets(vector<int> &members, unordered_map<int, vector<int>> &Member2AllSets)
{
	vector <int> dummy;
	return get_members_to_all_sets(members, dummy, Member2AllSets);

/*
	if (members.size() == 0) {
		for (auto &i : Id2Name) members.push_back(i.first);
	}

	Member2AllSets.clear();
#pragma omp parallel for
	for (int i=0; i<members.size(); i++) {
		int member = members[i];
		unordered_set<int> _used;
		queue<int> q;
		q.push(member);
		vector<int> v_sets;
		while (q.size() > 0) {
			int set_n = q.front();
			q.pop();
			if (_used.find(set_n) == _used.end()) {
				v_sets.push_back(set_n);
				_used.insert(set_n);
				if (Member2Sets.find(set_n) != Member2Sets.end())
					for (auto n : Member2Sets[set_n])
						q.push(n);
			}

		}
#pragma omp critical
		Member2AllSets[member] = v_sets;
	}
	*/
}

//-----------------------------------------------------------------------------------------------------------
// next function calculates for each member ALL the sets it is contained in , applying all set transitivity
void MedDictionary::get_members_to_all_sets(vector<int> &members, vector<int> &sets, unordered_map<int, vector<int>> &Member2AllSets)
{
	if (members.size() == 0) {
		for (auto &i : Id2Name) members.push_back(i.first);
	}

	set<int> use_sets(sets.begin(), sets.end());
	MLOG("use_sets size %d\n", use_sets.size());


	Member2AllSets.clear();
#pragma omp parallel for
	for (int i = 0; i<members.size(); i++) {
		int member = members[i];
		unordered_set<int> _used;
		queue<int> q;
		q.push(member);
		while (q.size() > 0) {
			int set_n = q.front();
			q.pop();
			if (_used.find(set_n) == _used.end()) {
				_used.insert(set_n);
				if (Member2Sets.find(set_n) != Member2Sets.end())
					for (auto n : Member2Sets[set_n])
						q.push(n);
			}

		}

		if (!use_sets.empty()) {
			vector<int> to_rm;
			for (auto &v : _used) if (use_sets.find(v) == use_sets.end()) to_rm.push_back(v);
			for (auto &v : to_rm) _used.erase(v);
		}

		vector<int> v_sets(_used.begin(), _used.end());

#pragma omp critical
		Member2AllSets[member] = v_sets;
	}
}


//-----------------------------------------------------------------------------------------------
void MedDictionary::get_regex_names(string regex_s, vector<string> &names)
{
	boost::regex regf(regex_s);
	for (auto &e : Id2Names) {
		for (auto &v : e.second) {
			if (boost::regex_match(v, regf)) {
				names.push_back(v);
				break;
			}
		}
	}
}

//-----------------------------------------------------------------------------------------------
void MedDictionary::get_regex_ids(string regex_s, vector<int> &ids)
{
	boost::regex regf(regex_s);
	for (auto &e : Id2Names) {
		for (auto &v : e.second) {
			if (boost::regex_match(v, regf)) {
				ids.push_back(e.first);
				break;
			}
		}
	}
}

//-----------------------------------------------------------------------------------------------
void MedDictionary::get_set_members(const string &set, vector<int> &members)
{
	members.clear();
	if (Name2Id.find(set) == Name2Id.end())
		return;
	int id = Name2Id[set];
	return get_set_members(id, members);
}

//-----------------------------------------------------------------------------------------------
void MedDictionary::get_set_members(int set_id, vector<int> &members)
{
	members.clear();
	if (Set2Members.find(set_id) == Set2Members.end())
		return;
	members = Set2Members[set_id];
}
//-----------------------------------------------------------------------------------------------
void MedDictionary::get_member_sets(const string &member, vector<int> &sets)
{
	sets.clear();
	if (Name2Id.find(member) == Name2Id.end())
		return;
	int id = Name2Id[member];
	return get_member_sets(id, sets);
}

//-----------------------------------------------------------------------------------------------
void MedDictionary::get_member_sets(int member_id, vector<int> &sets)
{
	sets.clear();
	if (Member2Sets.find(member_id) == Member2Sets.end())
		return;
	sets = Member2Sets[member_id];
}

//-----------------------------------------------------------------------------------------------
int MedDictionary::add_def(const string &fname, const string &name, int id)
{
	ofstream of;

	of.open(fname, ofstream::out | ofstream::app);

	if (!of) {
		MERR("MedDictionary: add_def: Can't open file %s\n", fname.c_str());
		return -1;
	}

	of << "DEF\t" << id << "\t" << name << "\n";

	of.close();

	return 0;
}

//-----------------------------------------------------------------------------------------------
int MedDictionary::add_set(const string &fname, const string &member_name, const string &set_name)
{
	ofstream of;

	of.open(fname, ofstream::out | ofstream::app);

	if (!of) {
		MERR("MedDictionary: add_set: Can't open file %s\n", fname.c_str());
		return -1;
	}

	of << "SET\t" << set_name << "\t" << member_name << "\n";
	of.close();

	return 0;
}


//-----------------------------------------------------------------------------------------------
int MedDictionary::prep_sets_lookup_table(const vector<string> &set_names, vector<char> &lut)
{
	// convert names to ids
	vector<int> sig_ids;
	for (auto &name : set_names) {
		int myid = id(name);
		if (myid >= 0)
			sig_ids.push_back(myid);
		else {
			MERR("ERROR in section %d : ", dict_id);
			for (auto s : section_name) MERR("%s,", s.c_str());
			MERR("\n");
			MTHROW_AND_ERR("prep_sets_lookup_table() : Found bad name [%s] :: not found in dictionary()\n", name.c_str());
		}
	}

	int max_id = 1;
	if (Id2Name.size() > 0)
		max_id = Id2Name.rbegin()->first;
	else
		MTHROW_AND_ERR("prep_sets_lookup_table() : Got an empty Id2Name...\n");

	lut.clear();
	lut.resize(max_id + 1, (char)0);

	/*
		for (int j=0; j<sig_ids.size(); j++) {
			MLOG("lut j=%d sig %s %d\n", j, set_names[j].c_str(), sig_ids[j]);
			for (int i=min_id; i<=max_id; i++)
				if (lut[sig_ids[j]] == 0)
					if (is_in_set(i, sig_ids[j]))
						lut[i] = 1;
		}
	*/
	// below MUCH faster version than the previous using depth search and a queue
	for (int j = 0; j < sig_ids.size(); j++) {
		//		MLOG("lut j=%d sig %s %d\n", j, set_names[j].c_str(), sig_ids[j]);
		queue<int> q;
		q.push(sig_ids[j]);

		while (q.size() > 0) {
			int s = q.front();
			q.pop();
			lut[s] = 1;
			for (auto elem : Set2Members[s])
				if (lut[elem] == 0)
					q.push(elem);

		}

	}
	return 0;
}

//-----------------------------------------------------------------------------------------------
int MedDictionary::prep_sets_indexed_lookup_table(const vector<string> &set_names, vector<unsigned char> &lut)
{
	// convert names to ids
	vector<int> sig_ids;
	for (auto &name : set_names) {
		int myid = id(name);
		if (myid > 0)
			sig_ids.push_back(myid);
		else {
			MERR("ERROR in section %d : ", dict_id);
			for (auto s : section_name) MERR("%s,", s.c_str());
			MERR("\n");
			MTHROW_AND_ERR("prep_sets_indexed_lookup_table() : Found bad name %s :: not found in dictionary()\n", name.c_str());
		}
	}

	int max_id = Id2Name.rbegin()->first;

	lut.clear();
	lut.resize(max_id + 1, 0);

	/*
	for (int j=0; j<sig_ids.size(); j++) {
	MLOG("lut j=%d sig %s %d\n", j, set_names[j].c_str(), sig_ids[j]);
	for (int i=min_id; i<=max_id; i++)
	if (lut[sig_ids[j]] == 0)
	if (is_in_set(i, sig_ids[j]))
	lut[i] = 1;
	}
	*/
	// below MUCH faster version than the previous using depth search and a queue
	for (int j = 0; j < sig_ids.size(); j++) {
		//		MLOG("lut j=%d sig %s %d\n", j, set_names[j].c_str(), sig_ids[j]);
		queue<int> q;
		q.push(sig_ids[j]);

		while (q.size() > 0) {
			int s = q.front();
			q.pop();
			lut[s] = j + 1;
			for (auto elem : Set2Members[s])
				if (lut[elem] == 0)
					q.push(elem);

		}

	}
	return 0;
}


//-----------------------------------------------------------------------------------------------
void MedDictionary::push_new_def(string name, int id)
{ 
	lock_guard<mutex> guard(lock_dict_changes);
	Name2Id[name] = id;
	Id2Name[id] = name; 
	Id2Names[id].push_back(name); 
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------------------------
int MedDictionarySections::read(const string &fname)
{
	// first we open the file and search for a SECTION statement within the first 100 lines
	// if there isn't one the section is "DEFAULT"

	string section_name = "DEFAULT";

	ifstream inf(fname);

	if (!inf) {
		MTHROW_AND_ERR("MedDictionarySections: read: Can't open file %s\n", fname.c_str());
	}

	MLOG_D("MedDictinarySections: read: reading dictionary file %s\n", fname.c_str());
	string curr_line;
	int n_line = 0;
	while (n_line < 100 && getline(inf, curr_line)) {
		n_line++;
		if ((curr_line.size() > 1) && (curr_line[0] != '#')) {

			if (curr_line[curr_line.size() - 1] == '\r')
				curr_line.erase(curr_line.size() - 1);

			vector<string> fields;
			split(fields, curr_line, boost::is_any_of("\t"));

			if (fields.size() >= 2) {
				if (fields[0].compare(0, 7, "SECTION") == 0) {
					section_name = fields[1];
					break;
				}
			}
		}
	}

	inf.close();

	vector<string> snames;
	split(snames, section_name, boost::is_any_of(" ,;:/"));
	int is_in = -1;
	for (int i = 0; i < snames.size(); i++)
		if (SectionName2Id.find(snames[i]) != SectionName2Id.end()) {
			is_in = SectionName2Id[snames[i]];
		}

	if (is_in < 0) {
		int section_id = (int)sections_names.size();
		sections_names.push_back(snames[0]);
		is_in = section_id;
		SectionName2Id[snames[0]] = section_id;
		MedDictionary dummy;
		dicts.push_back(dummy);
		section_fnames.push_back(vector<string>());
	}

	for (int i = 0; i < snames.size(); i++)
		SectionName2Id[snames[i]] = is_in;

	int section_id = SectionName2Id[snames[0]];
	section_fnames[section_id].push_back(fname);

	dicts[section_id].dict_id = section_id;
	for (auto s : snames)
		dicts[section_id].section_name.insert(s);

	return (dicts[section_id].read(fname));
}

//-----------------------------------------------------------------------------------------------
int MedDictionarySections::read(vector<string> &dfnames)
{
	int rc = 0;
	for (int i = 0; i < dfnames.size(); i++) {
		rc += read(dfnames[i]);
	}
	return rc;
}

//-----------------------------------------------------------------------------------------------
int MedDictionarySections::read(string path, vector<string> &dfnames)
{
	int rc = 0;
	for (int i = 0; i < dfnames.size(); i++) {
		string fname = (path == "" || (dfnames[i].length() > 0 && dfnames[i][0] == '/')) ? dfnames[i] : path + "/" + dfnames[i];
		try {
			rc += read(fname);
		}
		catch (...) {
			MERR("Error in reading Dictionary %s\n", fname.c_str());
			throw;
		}
	}
	return rc;
}

//------------------------------------------------------------------------------------------------------
void MedDictionarySections::add_section(string new_section_name) 
{ 
	lock_guard<mutex> guard(lock_dict_changes);
	MedDictionary dummy;
	dicts.push_back(dummy); 
	SectionName2Id[new_section_name] = (int)dicts.size() - 1; 
}

//------------------------------------------------------------------------------------------------------
void MedDictionarySections::connect_to_section(string new_section_name, int section_id) 
{ 
	lock_guard<mutex> guard(lock_dict_changes);
	SectionName2Id[new_section_name] = section_id;
}