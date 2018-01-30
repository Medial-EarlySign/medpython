#pragma once
//
// MedCrypt contains tools to sign files and test them, and to encrypt/decrypt files
// Main usage: to enable a library to build basic anti-tamperting mechanisms.
//
#include <string>
#include <vector>
#include <sstream>

using namespace std;

class MedSignFiles
{
	public:
		int get_signature(const string &fname, const string& signature_type, string &out_signature);

		static int calc_md5(const vector<char> &in, string &md5_out);

		// helpers
		static int read_file_to_char_vec(const string &fname, vector<char> &data);

		static string int_to_hex(int i)
		{
			std::stringstream stream;
			//stream << std::setfill('0') << std::setw(sizeof(int)*2) << std::hex << i;
			stream << std::hex << i;
			return stream.str();
		}
	


};