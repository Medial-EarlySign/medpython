//
// MedArchiver.h :: A model for archiving general data for learn/predict
//

#ifndef __MED_ARCHIVER_H__
#define __MED_ARCHIVER_H__

#include <vector>
#include <algorithm>
#include <cstring>
#include "MedAlgo/MedAlgo/MedAlgo.h"

using namespace std;

class MedArchiver {	
	private:
		vector<unsigned char> data ;
		size_t index ;
		size_t last_write_index ;
		size_t last_read_index;
		
	public:
		MedArchiver() {index = last_write_index = last_read_index =  0; data.clear(); };
		
		// Control
		void reset() {index = 0;}

		// Add
		template <class T> size_t add (T& newData) ;
		template <class T> size_t add (vector<T>& newData) ;
		template <class T, class S> size_t add (map <T,S> & newData);
		void add (const unsigned char *newData, const size_t& size) ;
		size_t add(string& info) ;
		size_t add(int& value) ;
		size_t add(float& value) ;
		size_t add(double& value) ;
		size_t add(char& value) ;
		size_t add(MedPredictor *predictor) ;

		// Extract
		template <class T> size_t get (T& newData) ; 
		template <class T> size_t get (vector<T>& newData) ;
		template <class T, class S> size_t get (map <T,S> & newData);
		void get(unsigned char *newData, const size_t& size) ;
		size_t get(string &info) ;
		size_t get(int & value) ;
		size_t get(float & value) ;
		size_t get(double & value) ;
		size_t get(char& value) ;
		size_t get(MedPredictor* &predictor) ;
	
		// Read/Write
		int read(const string& file_name) ;
		int write(const string& file_name) ;

} ;

#include "MedArchiver_imp.h"

#endif