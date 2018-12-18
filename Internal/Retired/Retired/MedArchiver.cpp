#define _CRT_SECURE_NO_WARNINGS

//
// MedArchiver
//

#include "MedArchiver.h"

//Add
void MedArchiver::add (const unsigned char *newData, const size_t& size) {
	
	data.insert(data.end(),newData,newData+size) ;
	index += size ;
}

size_t MedArchiver::add(string& info) {
	size_t size = info.length() ;
	add((unsigned char *) &size, sizeof(size_t)) ;
	add((unsigned char *) info.c_str(),size+1) ;

	size_t tot_size = sizeof(size_t) + size + 1 ;
	return tot_size ;
}

size_t MedArchiver::add(int& value) {
	size_t size = sizeof(int) ;
	data.insert(data.end(),(unsigned char *) &value,((unsigned char *) &value)+size ) ;
	index += size ;
	return size ;
}


size_t MedArchiver::add(char& value) {
	size_t size = sizeof(char) ;
	data.insert(data.end(),(unsigned char *) &value,((unsigned char *) &value)+size ) ;
	index += size ;
	return size ;
}


size_t MedArchiver::add(float& value) {
	size_t size = sizeof(float) ;
	data.insert(data.end(),(unsigned char *) &value,((unsigned char *) &value)+size  ) ;
	index += size ;
	return size ;
}



size_t MedArchiver::add(double& value) {
	size_t size = sizeof(double) ;
	data.insert(data.end(),(unsigned char *) &value,((unsigned char *) &value)+size ) ;
	index += size ;
	return size ;
}


// Add classifier_type + predictor
size_t MedArchiver::add(MedPredictor *predictor) {

	MedPredictorTypes type = predictor->classifier_type ;
	add((unsigned char *) &type, sizeof(MedPredictorTypes)) ;
	size_t size =  sizeof(MedPredictorTypes) ;

	size += add(*predictor) ;
	
	return size ;
}


// Extract
void MedArchiver::get(unsigned char *newData, const size_t& size) {
	
	memcpy(newData,&(data[index]),size) ;
	index += size ;

}

size_t MedArchiver::get(string& info) {
	size_t size ;
	get((unsigned char *)&size, sizeof(size_t)) ;

	vector<unsigned char> temp(size + 1) ;
	get(&(temp[0]),size+1) ;
	info = string((char *) &(temp[0])) ;

	size_t tot_size = sizeof(size_t) + size + 1 ;
	return tot_size ;
}

size_t MedArchiver::get(int & value) {
	size_t size = sizeof(int) ;
	memcpy(&value,&(data[index]),size) ;
	index += size ;
	return size ;
}


size_t MedArchiver::get(char & value) {

	size_t size = sizeof(char) ;
	memcpy(&value,&(data[index]),size) ;
	index += size ;
	return size ;
}

size_t MedArchiver::get(float & value) {
	size_t size = sizeof(float) ;
	memcpy(&value,&(data[index]),size) ;
	index += size ;
	return size ;
}

size_t MedArchiver::get(double & value) {
	size_t size = sizeof(double) ;
	memcpy(&value,&(data[index]),size) ;
	index += size ;
	return size ;
}



size_t MedArchiver::get(MedPredictor* &predictor) {

	size_t tot_size = 0 ; 
	MedPredictorTypes type ;
	get((unsigned char *)&type, sizeof(MedPredictorTypes)) ;
	size_t size = sizeof(MedPredictorTypes) ;

	tot_size += size ;

	predictor = MedPredictor::make_predictor(type) ;
	if (predictor == NULL)
		return -1 ;

	size = get(*predictor) ;
	tot_size += size ;

	return tot_size ;
}

// Read/Write
int MedArchiver::read(const string& file_name) {
	
	FILE *fp = fopen(file_name.c_str(),"rb") ;
	if (fp == NULL) {
		fprintf(stderr,"Cannot open %s for binary reading\n",file_name.c_str()) ;
		return -1 ;
	}

	int model_size ;
	if (fread(&model_size,sizeof(int),1,fp) != 1) {
		fprintf(stderr,"Cannot read size from %s\n",file_name.c_str()) ;
		return -1 ;
	}

	data.resize(model_size) ;

	size_t rc ;
	if ((rc = fread(&(data[last_read_index]),(size_t)model_size,1,fp)) != 1) {
		fprintf(stderr,"Cannot read model from %s (rc = %zd)\n",file_name.c_str(),rc) ;
		return -1 ;
	}

	last_read_index += model_size ;

	fclose(fp) ;

	return 0 ;
}

int MedArchiver::write(const string& file_name) {

	FILE *fp = fopen(file_name.c_str(),"wb") ;
	if (fp == NULL) {
		fprintf(stderr,"Cannot open %s for binary writing\n",file_name.c_str()) ;
		return -1 ;
	}

	int model_size = (int) (index - last_write_index) ;
	if (fwrite(&model_size,sizeof(int),1,fp) != 1) {
		fprintf(stderr,"Cannot write to %s\n",file_name.c_str()) ;
		return -1 ;
	}


	size_t rc;
	if (model_size) {
		if ((rc = (int)fwrite(&(data[last_write_index]), model_size, 1, fp)) != 1) {
			fprintf(stderr, "Cannot write data of size %d to %s model (rc = %lld) \n", model_size, file_name.c_str(), rc);
			return -1;
		}
	}
	fclose(fp) ;

	last_write_index = index ;
	return 0 ;
}