#define _CRT_SECURE_NO_WARNINGS

// Input/Output
#include "medial_utilities.h"

int read_blob(const char *file_name, unsigned char **data) {
	
	FILE *fp = safe_fopen(file_name, "rb", false) ;
	if (fp == NULL) {
		fprintf(stderr,"Cannot open %s for reading\n",file_name) ;
		return -1 ;
	}

	int size ;
	if (fread(&size,sizeof(int),1,fp) != 1) {
		fprintf(stderr,"Reading from %s failed\n",file_name) ;
		return -1 ;
	}

	if (((*data) = (unsigned char *) malloc (size)) == NULL) {
		fprintf(stderr,"Allocation of size %d failed\n",size) ;
		return -1 ;
	}

	if (fread(*data,1,size,fp) != size) {
		fprintf(stderr,"Reading from %s failed\n",file_name) ;
		return -1 ;
	}

	fclose(fp) ;

	return size ;
}