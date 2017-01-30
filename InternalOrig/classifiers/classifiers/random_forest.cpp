// Random Forest functions

#define _CRT_SECURE_NO_WARNINGS

#include "classifiers.h"

// Single tree prediction
int tree_pred(rf_tree tr, double *values) {

	int inode = 0 ;
	while (tr.nodes[inode].status == 1) {
		if (values[tr.nodes[inode].svar - 1] >= tr.nodes[inode].sval)
			inode = tr.nodes[inode].right - 1 ;
		else
			inode = tr.nodes[inode].left - 1 ;
	}

	return tr.nodes[inode].pred - 1 ;
}

// Forest prediction
void forest_pred(random_forest *rf, double *values, int *counts, int nclass) {

	for (int i=0; i<nclass; i++)
		counts[i] = 0 ;

	for (int i=0; i<rf->ntrees; i++)
		counts[tree_pred((rf->trees)[i],values)]++ ;

	return ;
}

// Read from file
int read_forest(char *file_name, random_forest *rf) {

	unsigned char *data ;
	if (read_blob(file_name,&data) == -1) {
		fprintf(stderr,"Reading RF model from \'%s\' failed\n",file_name) ;
		return -1 ;
	}

	if (rf_deserialize(data, rf) == -1) {
		fprintf(stderr,"DeSerialization of RF model failed\n") ;
		return -1 ;
	}

	free(data) ;
	return 0 ;
}

int read_text_forest(char *file_name, random_forest *rf) {

	FILE *fin = safe_fopen(file_name,"r", false) ;
	if (fin == NULL) {
		fprintf(stderr,"Cannot open %s for reading\n",file_name) ;
		return -1 ;
	}

	int ntrees ;
	char dummy[256] ;
	if (fscanf(fin,"%s %d\n",dummy,&ntrees) != 2) {
		fprintf(stderr,"Cannot read ntrees\n") ;
		return -1 ;
	}

	rf->ntrees = ntrees ;
	rf->trees = (rf_tree *) malloc(ntrees*sizeof(rf_tree)) ;
	if (rf->trees == NULL) {
		fprintf(stderr,"Allocation of forest failed\n") ;
		return -1 ;
	}

	for (int i=0; i<ntrees; i++) {
		int nnodes ;
		if (fscanf(fin,"%s %d\n",dummy,&nnodes) != 2) {
			fprintf(stderr,"Cannot read number of nodes for tree %d\n",i) ;
			return -1 ;
		}

		(rf->trees)[i].nnodes = nnodes ;
		(rf->trees)[i].nodes = (rf_node *) malloc(nnodes*sizeof(rf_node)) ;
		if ((rf->trees)[i].nodes == NULL) {
			fprintf(stderr,"Allocation of tree %d failed\n",i) ;
			return -1 ;
		}

		// Header
		fscanf(fin,"left daughter right daughter split var split point status prediction\n") ;

		for (int j=0; j<nnodes; j++) {
			int id ;
			int left,right,svar,status,pred ;
			double sval ;

			if (fscanf(fin,"%d %d %d %d %lf %d %d\n",&id,&left,&right,&svar,&sval,&status,&pred) != 7) {
				fprintf(stderr,"Cannot read info for node %d at tree %d\n",j+1,i) ;
				return -1 ;
			}

			if (id != j+1) {
				fprintf(stderr,"inconssistent line for tree %d node %d\n",i,j+1) ;
				return -1 ;
			}

			(rf->trees)[i].nodes[j].left = left ;
			(rf->trees)[i].nodes[j].right = right ;
			(rf->trees)[i].nodes[j].svar = svar ;
			(rf->trees)[i].nodes[j].sval = sval ;
			(rf->trees)[i].nodes[j].status = status ;
			(rf->trees)[i].nodes[j].pred = pred ;
		}
	}

	return 0 ;
}

// (De)Serialization
size_t get_rf_size(random_forest *rf) {

	size_t size = 2 * sizeof(int) ; // nTrees + nrNodes
	int ntrees = rf->ntrees ;

	int nrnodes = 0 ;
	for (int i=0 ; i<ntrees; i++) {
		if ((rf->trees)[i].nnodes > nrnodes)
			nrnodes = (rf->trees)[i].nnodes ;
	}
	int vec_size = nrnodes * ntrees ;

	size += ntrees * sizeof(int) ; // nNodes 
	size += vec_size * sizeof(int) ; // Status
	size += 2 * vec_size * sizeof(int) ; // tree (doughter)
	size += vec_size * sizeof(int) ; // pred 
	size += vec_size * sizeof(int) ; // var
	size += vec_size * sizeof(double) ; // value ;

	return size ;
}

int rf_serialize(random_forest *forest, unsigned char *rf_data) {
	
	size_t idx = 0 ;
	memcpy(rf_data,&(forest->ntrees),sizeof(int)) ; idx += sizeof(int) ;
	
	int nrnodes = 0 ;
	for (int i=0 ; i<forest->ntrees; i++) {
		if ((forest->trees)[i].nnodes > nrnodes)
			nrnodes = (forest->trees)[i].nnodes ;
	}

	memcpy(rf_data+idx,&nrnodes,sizeof(int)) ; idx += sizeof(int) ;

	int size = nrnodes * forest->ntrees ;
	
	// Point to data
	size_t nnodes_ptr = idx ; idx += forest->ntrees*sizeof(int) ;
	size_t status_ptr = idx ; idx += size*sizeof(int) ;
	size_t tree_ptr = idx   ; idx += 2*size*sizeof(int) ;
	size_t pred_ptr = idx   ; idx += size*sizeof(int) ;
	size_t var_ptr = idx    ; idx += size*sizeof(int) ;
	size_t value_ptr = idx  ; idx += size*sizeof(double) ;

	// Collect from trees
	for (int i=0; i<forest->ntrees; i++) {
		memcpy(rf_data + nnodes_ptr + i*sizeof(int),&(forest->trees)[i].nnodes,sizeof(int)) ;

		for (int j=0; j<(forest->trees)[i].nnodes; j++) {
			memcpy(rf_data + tree_ptr + (2*nrnodes*i + j)*sizeof(int),&((forest->trees)[i].nodes[j].left),sizeof(int)) ;
			memcpy(rf_data + tree_ptr + (2*nrnodes*i + nrnodes + j)*sizeof(int),&((forest->trees)[i].nodes[j].right),sizeof(int)) ;
			memcpy(rf_data + var_ptr + (nrnodes*i+j)*sizeof(int),&((forest->trees)[i].nodes[j].svar),sizeof(int)) ;
			memcpy(rf_data + value_ptr + (nrnodes*i+j)*sizeof(double),&((forest->trees)[i].nodes[j].sval),sizeof(double)) ;
			memcpy(rf_data + status_ptr + (nrnodes*i+j)*sizeof(int),&((forest->trees)[i].nodes[j].status),sizeof(int)) ;
			memcpy(rf_data + pred_ptr + (nrnodes*i+j)*sizeof(int),&((forest->trees)[i].nodes[j].pred),sizeof(int)) ;
		}
	}

	return (int) idx ;

}

int rf_deserialize(unsigned char *rf_data, random_forest *forest) {

	size_t idx = 0 ;
	memcpy(&(forest->ntrees),rf_data,sizeof(int)) ; idx += sizeof(int) ;

	int nrnodes ;
	memcpy(&nrnodes,rf_data+idx,sizeof(int)) ; idx += sizeof(int) ;
	int size = nrnodes * forest->ntrees ;

	// Point to data
	size_t nnodes_ptr = idx ; idx += forest->ntrees*sizeof(int) ;
	size_t status_ptr = idx ; idx += size*sizeof(int) ;
	size_t tree_ptr = idx   ; idx += 2*size*sizeof(int) ;
	size_t pred_ptr = idx   ; idx += size*sizeof(int) ;
	size_t var_ptr = idx    ; idx += size*sizeof(int) ;
	size_t value_ptr = idx  ; idx += size*sizeof(double) ;

	// ReArrange in tree
	(forest->trees) = (rf_tree *) malloc(forest->ntrees*sizeof(rf_tree)) ;
	if (forest->trees == NULL) {
		fprintf(stderr,"Allocation of forest failed\n") ;
		return -1 ;
	}

	for (int i=0; i<forest->ntrees; i++) {
		memcpy(&((forest->trees)[i].nnodes),rf_data + nnodes_ptr + i*sizeof(int),sizeof(int)) ;
		(forest->trees)[i].nodes = (rf_node *) malloc((forest->trees)[i].nnodes*sizeof(rf_node)) ;
		if ((forest->trees)[i].nodes == NULL) {
			fprintf(stderr,"Allocation of tree %d failed\n",i) ;
			return -1 ;
		}

		for (int j=0; j<(forest->trees)[i].nnodes; j++) {
			memcpy(&((forest->trees)[i].nodes[j].left),rf_data + tree_ptr + (2*nrnodes*i + j)*sizeof(int),sizeof(int)) ;
			memcpy(&((forest->trees)[i].nodes[j].right),rf_data + tree_ptr + (2*nrnodes*i + nrnodes + j)*sizeof(int),sizeof(int)) ;
			memcpy(&((forest->trees)[i].nodes[j].svar),rf_data + var_ptr + (nrnodes*i+j)*sizeof(int),sizeof(int)) ;
			memcpy(&((forest->trees)[i].nodes[j].sval),rf_data + value_ptr + (nrnodes*i+j)*sizeof(double),sizeof(double)) ;
			memcpy(&((forest->trees)[i].nodes[j].status),rf_data + status_ptr + (nrnodes*i+j)*sizeof(int),sizeof(int)) ;
			memcpy(&((forest->trees)[i].nodes[j].pred),rf_data + pred_ptr + (nrnodes*i+j)*sizeof(int),sizeof(int)) ;
		}
	}

	return (int) idx ;
}

int write_forest(random_forest *forest, char *file_name) {

	FILE *fp = safe_fopen(file_name, "wb", false) ;
	if (fp==NULL) {
		fprintf(stderr,"Cannot open %s for writing\n",file_name) ;
		return -1 ;
	}

	int size = (int) get_rf_size(forest) ;
	unsigned char *rf_data  = (unsigned char *) malloc(size) ;
	if (rf_data == NULL) {
		fprintf(stderr,"Allocation of RandomForest failed\n") ;
		return -1 ;
	}

	if (rf_serialize(forest,rf_data) != size) {
		fprintf(stderr,"Serialization of RandomForest failed\n") ;
		return -1 ;
	}

	if (fwrite(rf_data,1,size,fp) != size) {
		fprintf(stderr,"Cannot write to %s\n",file_name) ;
		return -1 ;
	}

	return 0 ;
}

void print_forest(random_forest *forest, FILE *fp) {
	
	fprintf(fp,"Ntrees = %d\n",forest->ntrees) ;

	for (int i=0; i<forest->ntrees; i++) {
		fprintf(fp,"Tree %d : Nnodes = %d\n",i,forest->trees[i].nnodes) ;

		for (int j=0; j<forest->trees[i].nnodes; j++) {
			fprintf(fp,"Tree %d Node %d : %d %d %d %.16g %d %d\n",i,j,forest->trees[i].nodes[j].left,forest->trees[i].nodes[j].right,
					forest->trees[i].nodes[j].svar,forest->trees[i].nodes[j].sval,forest->trees[i].nodes[j].status,forest->trees[i].nodes[j].pred) ;
		}
	}

	return ;
}

// Cleaning
void clear_random_forest(random_forest *forest) {

	for (int i=0; i<forest->ntrees; i++)
		free((forest->trees)[i].nodes) ;

	free(forest->trees) ;
}