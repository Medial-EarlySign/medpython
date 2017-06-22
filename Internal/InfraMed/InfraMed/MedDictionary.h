//
// MedDictionary.h - classes for dictionary and sets handling
//

#ifndef __MED__DICTIONARY__H__
#define __MED__DICTIONARY__H__

#include <string>
#include <vector>
#include <map>

using namespace std;

//#define	MAX_DICT_ID_NUM	10000000

class DLLEXTERN  MedDictionary {
	
	public:
		vector<string> fnames;
		map<string,int> Name2Id;			// we allow several names for each Id, the first one read will be the "official" name
		map<int,string> Id2Name;			// from an id to the "official" name
		map<int, vector<string>> Id2Names;	// a list of all the names given to an id.
		map<pair<int,int>,int> MemberInSet; // pair.first = set number, pair second = member number
		map<int,vector<int>> Set2Members;	// for each set a vector of all members of that set
		map<int,vector<int>> Member2Sets;	// for each id - a vector of all sets the id is a member of
		
		void clear() {fnames.clear(); Name2Id.clear(); MemberInSet.clear(); Set2Members.clear();}
		int read(const string &fname);
		int read(vector<string> &dfnames);
		int read(string path, vector<string> &dfnames) ;
		int id(const string &name);
		int id_list(vector<string> &names, vector<int> &ids);
		string name(int id);

		int is_in_set(int member_id, int set_id);
		int is_in_set(const string& member, const string& set_name);
		int is_in_set(int member_id, const string &set_name);
		int is_in_set(const string& member, int set_id);

		int prep_sets_lookup_table(const vector<string> &set_names, vector<char> &lut);
		int prep_sets_indexed_lookup_table(const vector<string> &set_names, vector<unsigned char> &lut);

		void get_set_members(const string &set, vector<int> &members);
		void get_set_members(int set_id, vector<int> &members);
		void get_member_sets(const string &member, vector<int> &sets);
		void get_member_sets(int member_id, vector<int> &sets);

		int add_def(const string &fname, const string &name, int id);
		int add_set(const string &fname, const string &member_name, const string &set_name);

	private:
		map<string, int> used;

};

//
// next class is a natural expansion of the dictionary class
// There are several "namespaces" or sections of dictionaries.
// The main one is the "DEFAULT" section, which always holds the signal ids (at least)
// Each dictionary file can be assigned to a section using the SECTION <name> statement in its first line (can be in one of the 100 first lines to make it easier)
// Each section has its OWN coding system, and hence can use values used in other sections.
// When using a dictionary, using it as done up to now (rep.dict.id()) will simply use the DEFAULT section.
// When we want to use a section we can use an API that has also the section name as input (better call section names with the name of the signal used)
// We can also set a certain section as the default section (and reset back to "DEFAULT")
//

class DLLEXTERN MedDictionarySections {

//class MedDictionarySections {

	public:
	vector<string> sections_names;	// [0] is always "DEFAULT"
	vector<vector<string>> section_fnames;
	vector<MedDictionary> dicts;	
	map<string, int>	SectionName2Id;
	int default_section;
	int curr_section;
	int read_state;

	// clearing
	void clear() { sections_names.clear(); section_fnames.clear(); dicts.clear(); SectionName2Id.clear(); default_section = 0;  curr_section = 0; init(); }

	// set curr/default
	int set_curr_dict(int id) { curr_section = id; }
	int set_curr_dict(const string name) { curr_section = SectionName2Id[name]; }
	int reset() { curr_section = default_section; }

	// initializations
	void init() { sections_names.push_back("DEFAULT");  section_fnames.push_back(vector<string>()); 
							  dicts.resize(1); SectionName2Id["DEFAULT"] = 0; 
							  default_section = 0; curr_section = 0; read_state = 0;
	}
	MedDictionarySections() { init(); }
	int read(const string &fname);
	int read(vector<string> &dfnames);
	int read(string path, vector<string> &dfnames);


	// in the following calls note that WE DO NOT CHECK IF section_id or section_name exist and appear in the map and correct range.
	// this is done for efficiency (saving an if for each call).

	// get a pointer to relevant dictionary
	MedDictionary *curr_dict() { return &dicts[curr_section]; }
	MedDictionary *dict(int section_id) { return &dicts[section_id]; }

	int section_id(const string name) { if (SectionName2Id.find(name) == SectionName2Id.end()) return 0; else return SectionName2Id[name]; }

	// use MedDictionary functions for a specific section id
	int id(int section_id, const string &name)							{ return dicts[section_id].id(name); }
	int id_list(int section_id,vector<string> &names, vector<int> &ids) { return dicts[section_id].id_list(names, ids); }
	string name(int section_id, int id)									{ return dicts[section_id].name(id); }

	int is_in_set(int section_id, int member_id, int set_id)						{ return dicts[section_id].is_in_set(member_id, set_id);}
	int is_in_set(int section_id, const string& member, const string& set_name)		{ return dicts[section_id].is_in_set(member, set_name); }
	int is_in_set(int section_id, int member_id, const string &set_name)			{ return dicts[section_id].is_in_set(member_id, set_name); }
	int is_in_set(int section_id, const string& member, int set_id)					{ return dicts[section_id].is_in_set(member, set_id); }

	void get_set_members(int section_id, const string &set, vector<int> &members)	{ return dicts[section_id].get_set_members(set, members); }
	void get_set_members(int section_id, int set_id, vector<int> &members)			{ return dicts[section_id].get_set_members(set_id, members); }
	void get_member_sets(int section_id, const string &member, vector<int> &sets)	{ return dicts[section_id].get_member_sets(member, sets); }
	void get_member_sets(int section_id, int member_id, vector<int> &sets)			{ return dicts[section_id].get_member_sets(member_id, sets); }

	// use curr_dict (typically default one)
	int id(const string &name)									{ return id(curr_section, name); }
	int id_list(vector<string> &names, vector<int> &ids)		{ return id_list(curr_section, names, ids); }
	string name(int id)											{ return name(curr_section, id); }

	int is_in_set(int member_id, int set_id)					{ return is_in_set(curr_section, member_id, set_id);}
	int is_in_set(const string& member, const string& set_name) { return is_in_set(curr_section, member, set_name); }
	int is_in_set(int member_id, const string &set_name)		{ return is_in_set(curr_section, member_id, set_name); }
	int is_in_set(const string& member, int set_id)				{ return is_in_set(curr_section, member, set_id); }

	void get_set_members(const string &set, vector<int> &members)	{ return get_set_members(curr_section, set, members); }
	void get_set_members(int set_id, vector<int> &members)			{ return get_set_members(curr_section, set_id, members); }
	void get_member_sets(const string &member, vector<int> &sets)	{ return get_member_sets(curr_section, member, sets); }
	void get_member_sets(int member_id, vector<int> &sets)			{ return get_member_sets(curr_section, member_id, sets); }


	int prep_sets_lookup_table(int section_id, const vector<string> &set_names, vector<char> &lut) { return dicts[section_id].prep_sets_lookup_table(set_names, lut); }
	int prep_sets_indexed_lookup_table(int section_id, const vector<string> &set_names, vector<unsigned char> &lut) { return dicts[section_id].prep_sets_indexed_lookup_table(set_names, lut); }
};


#endif
