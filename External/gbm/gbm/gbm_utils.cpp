#include "stdafx.h"
#include "gbm.h"
#include "gbm_utils.h"
#include "gbm_structs.h"
#include "gbm_predict.h"

#define XIDX(i,j,ncol) ((i)*(ncol) + (j))

struct compare_pair {
    bool operator()(const pair<double,int>& left, const pair<double,int>& right) {
		return (left.first < right.first) ;
    }
} ;


// Learn 
int get_gbm_predictor(double *x, double *y, double *w, int nrows, int ncols, double shrinkage, double bag_p, bool take_all_pos, int ntrees, int depth,
						int min_obs_in_node, gbm_tree *trees, double *initF, int ***rcSplits, int **rcSplitSizes, int *ncSplits,
						GBM_LossFunctions lossFunction, double alpha) {

	// Take care of bag_p, take_all_pos, and pos/neg ratio.
	int npos = 0 ;
	for (int i=0; i<nrows; i++) {
		if (y[i] > 0)
			npos ++ ;
	}

	double pos_p = (npos + 0.0)/nrows ;
	if (take_all_pos && bag_p < 2 * pos_p) {
		fprintf(stderr,"Prevalance (%f) is more than half of bag_p (%f) => Turning off take_all_pos\n",pos_p,bag_p) ;
		take_all_pos = false ;
//		bag_p = 1.0 ;
	}

	// Prepare
	int *order ;
	if((order = get_order(x,nrows,ncols))==NULL) {
		fprintf(stderr,"Ordering failed\n") ;
		return -1 ;
	}

	int *monotone = (int *) malloc(ncols*sizeof(int)) ;
	int *var_types = (int *) malloc(ncols*sizeof(int)) ;
	if (monotone == NULL || var_types == NULL) {
		fprintf(stderr,"Allocation failed\n") ;
		return -1 ;
	}
	memset (monotone,0,ncols*sizeof(int)) ;
	memset (var_types,0,ncols*sizeof(int)) ;

	CDataset learningSet ;
	if (learningSet.SetData(x,order,y,NULL,w,NULL,nrows,ncols,var_types,monotone) != GBM_OK) {
		fprintf(stderr,"DataSet initialization failed\n") ;
		return -1 ;
	}

	CGBM *pGBM ;
	pGBM = new CGBM() ;

	VEC_VEC_CATEGORIES vecSplitCodes;
	
	PCDistribution pDist;
	switch(lossFunction){
	case GBM_Loss_AdaBoost:
		pDist = new CAdaBoost();
		break;
	case GBM_Loss_Bernulli:
		pDist = new CBernoulli();
		break;
	case GBM_Loss_Laplace:
		pDist = new CLaplace();
		break;
	case GBM_Loss_Poisson:
		pDist = new CPoisson();
		break;
	case GBM_Loss_Quantile:
		pDist = new CQuantile(alpha);
		break;
	case GBM_Loss_Gaussian:
		pDist = new CGaussian();
		break;
	default:
		fprintf(stderr, "The chosen distribution for gbm is not valid!!\n");
		return(-1);
	}

	if (GBM_FAILED(pGBM->Initialize(&learningSet,pDist,shrinkage,nrows,bag_p,depth,min_obs_in_node,take_all_pos))) {
		fprintf(stderr,"Initialization failed\n") ;
		return -1 ;
	}

	double init_f;
	double *f = (double *) malloc(nrows*sizeof(double)) ;
	if (f==NULL) {
		fprintf(stderr,"Predictions allocation failed\n") ;
		return -1 ;
	}

	if (GBM_FAILED(pDist->InitF(learningSet.adY,learningSet.adMisc,learningSet.adOffset,learningSet.adWeight,init_f,nrows))) {
		fprintf(stderr,"Initialization of predictions failed\n") ;
		return -1 ;
	}
     
	for(int i=0; i < nrows; i++)
		f[i] = init_f ;


	double *train_err = (double *) malloc(ntrees*sizeof(double)) ;
	double *valid_err = (double *) malloc(ntrees*sizeof(double)) ;
	double *oob_improve = (double *) malloc(ntrees*sizeof(double)) ;
	if (train_err == NULL || valid_err == NULL || oob_improve == NULL) {
		fprintf(stderr,"Allocation of errors failed\n") ;
		return -1 ;
	}

	int c_nodes = 0 ;

	for (int i=0; i<ntrees; i++) {
		if (i%15 == 0)
			fprintf(stderr,"Starting Tree %d ... ",i) ;

		if (GBM_FAILED(pGBM->iterate(f,train_err[i],valid_err[i],oob_improve[i],c_nodes))) {	
			fprintf(stderr,"Tree %d failed\n",i) ;
			return -1 ;
		}

		if (GBM_FAILED(set_tree(pGBM,c_nodes,vecSplitCodes,trees,i))) {
			fprintf(stderr,"Trees management failed at iteration %d\n",i) ;
			return -1 ;
		}
	}
	fprintf(stderr,"\n") ;

	*initF = init_f ;

	*ncSplits = (int) vecSplitCodes.size() ;
	if (((*rcSplits) = (int **) malloc (vecSplitCodes.size()*sizeof(int *)))==NULL) {
		fprintf(stderr,"Allocation failed\n") ;
		return -1 ;
	}
	if (((*rcSplitSizes) = (int *) malloc (vecSplitCodes.size()*sizeof(int)))==NULL) {
		fprintf(stderr,"Allocation failed\n") ;
		return -1 ;
	}
	
    for(unsigned int i=0; i<vecSplitCodes.size(); i++) {
		(*rcSplitSizes)[i] = (int) vecSplitCodes[i].size() ;
		if (((*rcSplits)[i] = (int *) malloc(vecSplitCodes[i].size() * sizeof(int)))==NULL) {
			fprintf(stderr,"Allocation failed at %d\n",i) ;
			return -1 ;
		}

        if (GBM_FAILED(gbm_transfer_catsplits_to_R(i,vecSplitCodes,(*rcSplits)[i]))) {
			fprintf(stderr,"Splits management failed\n") ;
			return -1 ;
		}
    }

	free(order) ;
	free(monotone) ;
	free(var_types) ;
	free(train_err) ; 
	free(valid_err) ;
	free(oob_improve) ;

	return 0 ;
}

// Envelope - 
int get_gbm_predictor(double *x, double *y, double *w, int nrows, int ncols, double shrinkage, 
					  double bag_p, bool take_all_pos, int ntrees, int depth, int min_obs_in_node, 
					  full_gbm_learn_info_t *gbm_info, GBM_LossFunctions lossFunction, double alpha) {
	
	gbm_info->ntrees = ntrees ;
	gbm_info->trees = allocate_gbm_trees(ntrees) ;
	if (gbm_info->trees==NULL) {
		fprintf(stderr,"GBM allocation failed\n") ;
		return -1 ;
	}

	return get_gbm_predictor(x,y,w,nrows,ncols,shrinkage,bag_p,take_all_pos,ntrees,depth,min_obs_in_node, gbm_info->trees,&(gbm_info->initF),&(gbm_info->rcSplits),
		&(gbm_info->rcSplitSizes),&(gbm_info->ncSplits), lossFunction) ;
}

// Learn 
int get_gbm_regressor(double *x, double *y, double *w, int nrows, int ncols, double shrinkage, double bag_p, int ntrees, int depth, int min_obs_in_node,
						gbm_tree *trees, double *initF, int ***rcSplits, int **rcSplitSizes, int *ncSplits) {

	// Prepare
	int *order ;
	if((order = get_order(x,nrows,ncols))==NULL) {
		fprintf(stderr,"Ordering failed\n") ;
		return -1 ;
	}

	int *monotone = (int *) malloc(ncols*sizeof(int)) ;
	int *var_types = (int *) malloc(ncols*sizeof(int)) ;
	if (monotone == NULL || var_types == NULL) {
		fprintf(stderr,"Allocation failed\n") ;
		return -1 ;
	}
	memset (monotone,0,ncols*sizeof(int)) ;
	memset (var_types,0,ncols*sizeof(int)) ;

	CDataset learningSet ;
	if (learningSet.SetData(x,order,y,NULL,w,NULL,nrows,ncols,var_types,monotone) != GBM_OK) {
		fprintf(stderr,"DataSet initialization failed\n") ;
		return -1 ;
	}

	CGBM *pGBM ;
	pGBM = new CGBM() ;

	VEC_VEC_CATEGORIES vecSplitCodes;
	PCDistribution pDist = new CGaussian();

	if (GBM_FAILED(pGBM->Initialize(&learningSet,pDist,shrinkage,nrows,bag_p,depth,min_obs_in_node,false))) {
		fprintf(stderr,"Initialization failed\n") ;
		return -1 ;
	}

	double init_f;
	double *f = (double *) malloc(nrows*sizeof(double)) ;
	if (f==NULL) {
		fprintf(stderr,"Predictions allocation failed\n") ;
		return -1 ;
	}

	if (GBM_FAILED(pDist->InitF(learningSet.adY,learningSet.adMisc,learningSet.adOffset,learningSet.adWeight,init_f,nrows))) {
		fprintf(stderr,"Initialization of predictions failed\n") ;
		return -1 ;
	}
     
	for(int i=0; i < nrows; i++)
		f[i] = init_f ;


	double *train_err = (double *) malloc(ntrees*sizeof(double)) ;
	double *valid_err = (double *) malloc(ntrees*sizeof(double)) ;
	double *oob_improve = (double *) malloc(ntrees*sizeof(double)) ;
	if (train_err == NULL || valid_err == NULL || oob_improve == NULL) {
		fprintf(stderr,"Allocation of errors failed\n") ;
		return -1 ;
	}

	int c_nodes = 0 ;

	for (int i=0; i<ntrees; i++) {
		if (i%15 == 0)
			fprintf(stderr,"Starting Tree %d ... ",i) ;

		if (GBM_FAILED(pGBM->iterate(f,train_err[i],valid_err[i],oob_improve[i],c_nodes))) {	
			fprintf(stderr,"Tree %d failed\n",i) ;
			return -1 ;
		}

		if (GBM_FAILED(set_tree(pGBM,c_nodes,vecSplitCodes,trees,i))) {
			fprintf(stderr,"Trees management failed at iteration %d\n",i) ;
			return -1 ;
		}
	}
	fprintf(stderr,"\n") ;

	*initF = init_f ;

	*ncSplits = (int) vecSplitCodes.size() ;
	if (((*rcSplits) = (int **) malloc (vecSplitCodes.size()*sizeof(int *)))==NULL) {
		fprintf(stderr,"Allocation failed\n") ;
		return -1 ;
	}
	if (((*rcSplitSizes) = (int *) malloc (vecSplitCodes.size()*sizeof(int)))==NULL) {
		fprintf(stderr,"Allocation failed\n") ;
		return -1 ;
	}
	
    for(unsigned int i=0; i<vecSplitCodes.size(); i++) {
		(*rcSplitSizes)[i] = (int) vecSplitCodes[i].size() ;
		if (((*rcSplits)[i] = (int *) malloc(vecSplitCodes[i].size() * sizeof(int)))==NULL) {
			fprintf(stderr,"Allocation failed at %d\n",i) ;
			return -1 ;
		}

        if (GBM_FAILED(gbm_transfer_catsplits_to_R(i,vecSplitCodes,(*rcSplits)[i]))) {
			fprintf(stderr,"Splits management failed\n") ;
			return -1 ;
		}
    }

	free(order) ;
	free(monotone) ;
	free(var_types) ;
	free(train_err) ; 
	free(valid_err) ;
	free(oob_improve) ;

	return 0 ;
}

// Envelope - 
int get_gbm_regressor(double *x, double *y, double *w, int nrows, int ncols, double shrinkage, double bag_p, int ntrees, int depth,int min_obs_in_node, full_gbm_learn_info_t *gbm_info) {

	gbm_info->ntrees = ntrees ;
	gbm_info->trees = allocate_gbm_trees(ntrees) ;
	if (gbm_info->trees==NULL) {
		fprintf(stderr,"GBM allocation failed\n") ;
		return -1 ;
	}

	return get_gbm_regressor(x,y,w,nrows,ncols,shrinkage,bag_p,ntrees,depth,min_obs_in_node,gbm_info->trees,&(gbm_info->initF),&(gbm_info->rcSplits),
		&(gbm_info->rcSplitSizes),&(gbm_info->ncSplits)) ;
}

// Predict
int gbm_predict(double *x, int nrows, int ncols, const gbm_tree *trees, int ntrees, double init_f, int **rcSplits, double *preds) {

	int *var_types = (int *) malloc(ncols*sizeof(int)) ;
	if (var_types == NULL) {
		fprintf(stderr,"Allocation failed\n") ;
		return -1 ;
	}
	memset (var_types,0,ncols*sizeof(int)) ;
	
	if (GBM_FAILED(gbm_pred(x,nrows,ncols,ntrees,init_f,trees,rcSplits,var_types,preds))) {
		fprintf(stderr,"Prediction failed\n") ; ;
		return -1 ;
	}

	free(var_types) ;
	return 0 ;
}

int gbm_predict(double *x, int nrows, int ncols, int ntrees, const full_gbm_learn_info_t *gbm_info, double *preds) {

	if (ntrees > gbm_info->ntrees) {
		fprintf(stderr,"Required ntrees larger than availble\n") ;
		return -1;
	}

	return gbm_predict(x,nrows,ncols,gbm_info->trees,ntrees,gbm_info->initF,gbm_info->rcSplits,preds) ;
}

// Ordering
int *get_order(double *x, int nrows, int ncols) {

	int *order = (int *) malloc(nrows*ncols*sizeof(int)) ;
	if (order == NULL) {
		fprintf(stderr,"Allocation failed\n") ;
		return NULL ;
	}

	for (int icol=0; icol<ncols; icol++) {
		vector<pair<double,int> > vec ;
		for (int irow=0; irow<nrows; irow++) {
			pair<double,int> value(x[XIDX(icol,irow,nrows)],irow) ;
			vec.push_back(value) ;
		}

		sort(vec.begin(),vec.end(),compare_pair()) ;

		for (int irow=0; irow<nrows; irow++)
			order[XIDX(icol,irow,nrows)] = vec[irow].second ;
	}


	return order ;
}

// Clear
void clear_gbm_info(full_gbm_learn_info_t *gbm_info) {

	for (int i=0; i<gbm_info->ncSplits; i++)
		 free((gbm_info->rcSplits)[i]) ;

	if (gbm_info->ncSplits) {
		free(gbm_info->rcSplits) ;
		free(gbm_info->rcSplitSizes) ;
	}

	free_gbm_trees(gbm_info->trees) ;
}

// Read from File
int read_full_gbm_info(full_gbm_learn_info_t *gbm_info, char *fname) {

	FILE *fp ;
	if ((fp = fopen(fname,"rb"))==NULL) {
		fprintf(stderr,"Cannot open %s for reading\n",fname) ;
		return -1 ;
	}

	int rc = read_full_gbm_info(gbm_info,fp) ;
	fclose(fp) ;
	return rc ;
}

int read_full_gbm_info(full_gbm_learn_info_t *gbm_info, FILE *fp) {

	int ntrees ;
	if (fscanf(fp,"Ntrees = %d\n",&ntrees) != 1) {
		fprintf(stderr,"Reading of ntrees failed\n") ;
		return -1 ;
	}
	gbm_info->ntrees = ntrees ;

	gbm_info->trees = allocate_gbm_trees(ntrees) ;
	if (gbm_info->trees==NULL) {
		fprintf(stderr,"Trees allocation failed\n") ;
		return -1 ;
	}

	int itree ;
	for (int i=0; i<ntrees; i++) {
		if (fscanf(fp,"Tree %d\n",&itree) != 1 || itree != i) {
			fprintf(stderr,"Tree %d header cannot be read\n",i) ;
			return -1 ;
		}

		if (read_gbm_tree(fp,gbm_info->trees,i)==-1) {
			fprintf(stderr,"Reading of tree %d failed\n",i) ;
			return -1 ;
		}
	}

	double initF ;
	if (fscanf(fp,"initF = %lf\n",&initF) != 1) {
		fprintf(stderr,"Reading of initF failed\n") ;
		return -1 ;
	}
	gbm_info->initF = initF ;

	int ncSplits ;
	if (fscanf(fp,"ncSplits = %d\n",&ncSplits) != 1) {
		fprintf(stderr,"Reading of ncSplits failed\n") ;
		return -1 ;
	}
	gbm_info->ncSplits = ncSplits ;
	
	gbm_info->rcSplits = (int **) malloc(ncSplits*sizeof(int *)) ;
	gbm_info->rcSplitSizes = (int *) malloc(ncSplits*sizeof(int)) ;
	if (gbm_info->rcSplits==NULL || gbm_info->rcSplitSizes==NULL) {
		fprintf(stderr,"Splits Allocation failed\n") ;
		return -1 ;
	}

	int ival,jval,val; 
	for (int i=0; i<ncSplits; i++) {
		if (fscanf(fp,"SplitSize %d = %d\n",&ival,&val)!=2 || ival != i) {
			fprintf(stderr,"Reading of SplitSize %d failed\n",i) ;
			return -1 ;
		}
		gbm_info->rcSplitSizes[i] = val ;
		if(((gbm_info->rcSplits)[i] = (int *) malloc(val*sizeof(int))) == NULL) {
			fprintf(stderr,"Splits %d allocation failed\n",i) ;
			return -1 ;
		}

		for (int j=0; j<(gbm_info->rcSplitSizes)[i]; j++) {
			if (fscanf(fp,"Split %d %d = %d\n",&ival,&jval,&val)!=3 || ival != i || jval != j) {
				fprintf(stderr,"Reading of Split %d,%d failed\n",i,j) ;
				return -1 ;
			}
			(gbm_info->rcSplits)[i][j] = val ;
		}
	}

	return 0 ;
}

// Write to File
int write_full_gbm_info(full_gbm_learn_info_t *gbm_info, char *fname) {

	FILE *fp ;
	if ((fp = fopen(fname,"wb"))==NULL) {
		fprintf(stderr,"Cannot open %s for writing\n",fname) ;
		return -1 ;
	}

	write_full_gbm_info(gbm_info,fp) ;
	fclose(fp) ;
	return 0 ;
}

void write_full_gbm_info(full_gbm_learn_info_t *gbm_info, FILE *fp) {

	fprintf(fp,"Ntrees = %d\n",gbm_info->ntrees) ;
	
	for (int i=0; i<gbm_info->ntrees; i++) {
		fprintf(fp,"Tree %d\n",i) ;
		print_gbm_tree(fp,gbm_info->trees,i) ;
	}
	
	fprintf(fp,"initF = %f\n",gbm_info->initF) ;
	fprintf(fp,"ncSplits = %d\n",gbm_info->ncSplits) ;
	
	for (int i=0; i<gbm_info->ncSplits; i++) {
		fprintf(fp,"SplitSize %d = %d\n",i,(gbm_info->rcSplitSizes)[i]) ;
		for (int j=0; j<(gbm_info->rcSplitSizes)[i]; j++)
			fprintf(fp,"Split %d %d = %d\n",i,j,(gbm_info->rcSplits)[i][j]) ;
	}

	return ;
}

void write_full_gbm_info(const string& prefix, full_gbm_learn_info_t *gbm_info, FILE *fp) {

	fprintf(fp,"%s: Ntrees = %d\n",prefix.c_str(),gbm_info->ntrees) ;
	
	for (int i=0; i<gbm_info->ntrees; i++) {
		fprintf(fp,"%s: Tree %d\n",prefix.c_str(),i) ;
		print_gbm_tree(prefix,fp,gbm_info->trees,i) ;
	}
	
	fprintf(fp,"%s: initF = %f\n",prefix.c_str(),gbm_info->initF) ;
	fprintf(fp,"%s: ncSplits = %d\n",prefix.c_str(),gbm_info->ncSplits) ;
	
	for (int i=0; i<gbm_info->ncSplits; i++) {
		fprintf(fp,"%s: SplitSize %d = %d\n",prefix.c_str(),i,(gbm_info->rcSplitSizes)[i]) ;
		for (int j=0; j<(gbm_info->rcSplitSizes)[i]; j++)
			fprintf(fp,"%s: Split %d %d = %d\n",prefix.c_str(),i,j,(gbm_info->rcSplits)[i][j]) ;
	}

	return ;
}

// (De)Serialize gbm_info
size_t get_gbm_info_size (full_gbm_learn_info_t *gbm_info) {

	size_t blob_size = sizeof(int) + sizeof(double) + sizeof(int) ; // nTrees + initF + ncSplits
	for (int i=0; i<gbm_info->ntrees; i++) {
		blob_size += sizeof(int) ; // Per Tree : Nnodes
		blob_size += gbm_info->trees[i].nnodes * (sizeof(int) + sizeof(double) + sizeof(int) + sizeof(int)) ; // Per Tree, Node : SplitVar + SplitPoint + LeftNode + RightNode
	}

	for (int i=0; i<gbm_info->ncSplits; i++) {
		blob_size += sizeof(int) ; // rcSplitSize
		blob_size += gbm_info->rcSplitSizes[i] * sizeof(int)	; // rcSplits
	}

	return blob_size ;
}
	
int gbm_serialize (full_gbm_learn_info_t *gbm_info, unsigned char *gbm_data) {

	// Insert structure to blob

	size_t idx = 0 ;
	memcpy(gbm_data+idx,&(gbm_info->ntrees),sizeof(int)) ; idx += sizeof(int) ;
	memcpy(gbm_data+idx,&(gbm_info->initF),sizeof(double)) ; idx += sizeof(double) ;
	memcpy(gbm_data+idx,&(gbm_info->ncSplits),sizeof(int)) ; idx += sizeof(int) ;

	for (int i=0; i<gbm_info->ntrees; i++) {
		memcpy(gbm_data+idx,&(gbm_info->trees[i].nnodes),sizeof(int)) ; idx += sizeof(int) ;
		memcpy(gbm_data+idx,gbm_info->trees[i].split_var,gbm_info->trees[i].nnodes * sizeof(int)) ; idx += gbm_info->trees[i].nnodes * sizeof(int) ;
		memcpy(gbm_data+idx,gbm_info->trees[i].split_point,gbm_info->trees[i].nnodes * sizeof(double)) ; idx += gbm_info->trees[i].nnodes * sizeof(double) ;
		memcpy(gbm_data+idx,gbm_info->trees[i].left_node,gbm_info->trees[i].nnodes * sizeof(int)) ; idx += gbm_info->trees[i].nnodes * sizeof(int) ;
		memcpy(gbm_data+idx,gbm_info->trees[i].right_node,gbm_info->trees[i].nnodes * sizeof(int)) ; idx += gbm_info->trees[i].nnodes * sizeof(int) ;
	}

	for (int i=0; i<gbm_info->ncSplits; i++) {
		memcpy(gbm_data+idx,&(gbm_info->rcSplitSizes[i]),sizeof(int)) ; idx += sizeof(int) ;
		memcpy(gbm_data+idx,gbm_info->rcSplits[i],gbm_info->rcSplitSizes[i] * sizeof(int)) ; idx += gbm_info->rcSplitSizes[i] * sizeof(int) ;
	}

	return (int) idx ;
}

int gbm_deserialize (unsigned char *gbm_data, full_gbm_learn_info_t *gbm_info) {

	size_t idx = 0 ;
	memcpy(&(gbm_info->ntrees),gbm_data+idx,sizeof(int)) ; idx += sizeof(int) ;
	memcpy(&(gbm_info->initF),gbm_data+idx,sizeof(double)) ; idx += sizeof(double) ;
	memcpy(&(gbm_info->ncSplits),gbm_data+idx,sizeof(int)) ; idx += sizeof(int) ;

	gbm_info->trees = (gbm_tree *) malloc(gbm_info->ntrees * sizeof(gbm_tree)) ;
	if (gbm_info->trees == NULL) {
		fprintf(stderr,"Allocation of gbm structure failed\n") ;
		return -1 ;
	}

	for (int i=0; i<gbm_info->ntrees; i++) {
		memcpy(&(gbm_info->trees[i].nnodes),gbm_data+idx,sizeof(int)) ; idx += sizeof(int) ;
		
		gbm_info->trees[i].split_var = (int *) malloc(gbm_info->trees[i].nnodes * sizeof(int)) ;
		gbm_info->trees[i].split_point = (double *) malloc(gbm_info->trees[i].nnodes * sizeof(double)) ;
		gbm_info->trees[i].left_node = (int *) malloc(gbm_info->trees[i].nnodes * sizeof(int)) ;
		gbm_info->trees[i].right_node = (int *) malloc(gbm_info->trees[i].nnodes * sizeof(int)) ;
		if (gbm_info->trees[i].split_var == NULL || gbm_info->trees[i].split_point == NULL || gbm_info->trees[i].left_node == NULL || gbm_info->trees[i].right_node == NULL) {
			fprintf(stderr,"Allocation of gbm structure failed\n") ;
			return -1 ;
		}

		memcpy(gbm_info->trees[i].split_var,gbm_data+idx, gbm_info->trees[i].nnodes * sizeof(int)); idx += gbm_info->trees[i].nnodes * sizeof(int) ;
		memcpy(gbm_info->trees[i].split_point,gbm_data+idx, gbm_info->trees[i].nnodes * sizeof(double)); idx += gbm_info->trees[i].nnodes * sizeof(double) ;
		memcpy(gbm_info->trees[i].left_node,gbm_data+idx, gbm_info->trees[i].nnodes * sizeof(int)); idx += gbm_info->trees[i].nnodes * sizeof(int) ;
		memcpy(gbm_info->trees[i].right_node,gbm_data+idx, gbm_info->trees[i].nnodes * sizeof(int)); idx += gbm_info->trees[i].nnodes * sizeof(int) ;
	}

	gbm_info->rcSplitSizes = (int *) malloc(gbm_info->ncSplits * sizeof(int)) ;
	if (gbm_info->rcSplitSizes == NULL) {
		fprintf(stderr,"Allocation of gbm structure failed\n") ;
		return -1 ;
	}

	for (int i=0; i<gbm_info->ncSplits; i++) {
		memcpy(&(gbm_info->rcSplitSizes[i]),gbm_data+idx,sizeof(int)) ; idx += sizeof(int) ;

		gbm_info->rcSplits[i] = (int *) malloc(gbm_info->rcSplitSizes[i] * sizeof(int)) ;
		if (gbm_info->rcSplits[i] == NULL) {
			fprintf(stderr,"Allocation of gbm structure failed\n") ;
			return -1 ;
		}
		memcpy(gbm_info->rcSplits[i],gbm_data+idx, gbm_info->rcSplitSizes[i] * sizeof(int)) ; idx += gbm_info->rcSplitSizes[i] * sizeof(int) ;
	}

	return (int) idx ;
}