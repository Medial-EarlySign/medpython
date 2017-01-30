#define _CRT_SECURE_NO_WARNINGS

#define _CRT_RAND_S

#ifdef __linux
#include <boost/random/mersenne_twister.hpp>
#endif

#include "stdlib.h"

#include "gbm/gbm/gbm_utils.h"
#include "classifiers.h" 

double *nfc;
int auccompare( const void *arg1, const void *arg2 )
{
	if (nfc[*(int *)arg1]<nfc[*(int *)arg2]) return(1);
	if (nfc[*(int *)arg1]>nfc[*(int *)arg2]) return(-1);
	return(0);
}


//int get_gbm_predictor(double *x, double *y, double *w, int nrows, int ncols, double shrinkage, double bag_p, bool take_all_pos, int ntrees, int depth,
//						gbm_tree *trees, double *initF, int ***rcSplits, int **rcSplitSizes, int *ncSplits) {

int get_gbm_predictor_ens(double *x1, double *y1, int train_size, int nftrs, int ens_size, full_gbm_ens_learn_info_t *full_info_ens, gbm_parameters *gbm_params)
{
	double *w1 = (double *) malloc(train_size*sizeof(double)) ;
	if (w1==NULL) {
		fprintf(stderr,"Allocation Failed\n") ;
		return -1 ;
	}

	int rc = get_gbm_predictor_ens(x1,y1,w1,train_size,nftrs,ens_size,full_info_ens,gbm_params) ;

	free(w1) ;
	return rc ;
}

int get_gbm_predictor_ens(double *x1, double *y1, double *w1, int train_size, int nftrs, int ens_size, full_gbm_ens_learn_info_t *full_info_ens, gbm_parameters *gbm_params)
{
#ifdef __linux
	boost::random::mt19937 mt_rng;
#endif

	full_info_ens->ens_size = ens_size ;

	int ens_ntrees = gbm_params->ntrees;
	double ens_shrinkage = gbm_params->shrinkage;
	double ens_bag_p = gbm_params->bag_p;
	int ens_depth = gbm_params->depth ;
	bool ens_take_all = gbm_params->take_all_pos ;
	int min_obs_in_node = gbm_params->min_obs_in_node ;

	double *newx1 = (double *) calloc(train_size*3*nftrs, sizeof(double));
	double *newx1half1 = (double *) calloc(train_size/2*nftrs , sizeof(double));
	double *newx1half2 = (double *) calloc(train_size/2*nftrs , sizeof(double)); 			
			
	double *y1half1 = (double *)calloc(train_size/2, sizeof(double));
	double *y1half2 = (double *)calloc(train_size/2, sizeof(double));
//	double *newx2 = (double *) calloc(test_size*nftrs , sizeof(double)); 

	if (newx1==NULL || newx1half1==NULL || newx1half2==NULL || y1half1==NULL || y1half2==NULL) {
		fprintf(stderr,"Allocation failed\n") ;
		return -1 ;
	}

	full_info_ens->alpha = (double *) malloc (ens_size*sizeof(double)) ;
	full_info_ens->full_info_tree = (full_gbm_learn_info_t *) malloc (ens_size*sizeof(full_gbm_learn_info_t)) ;

	double **useyhalf1 = (double **) malloc(ens_size*sizeof(double *)) ;
	double **useyhalf2 = (double **) malloc(ens_size*sizeof(double *)) ;
	//double **useytest = (double **) malloc(ens_size*sizeof(double *));

	if (full_info_ens->alpha==NULL || full_info_ens->full_info_tree==NULL || useyhalf1==NULL || useyhalf2==NULL) {
		fprintf(stderr,"Allocation failed\n") ;
		return -1 ;
	}
	
	// Do scramble....

	for(int it=0; it<train_size; it++)
	{
		int i1 = it;
		unsigned int randnum;
#ifdef __linux
		randnum = mt_rng();
#else
		rand_s(&randnum);
#endif
		int i2 = randnum%train_size;
		for(int j=0; j<nftrs; j++)
		{
			double tmp = x1[j*train_size+i1];
			x1[j*train_size+i1] = x1[j*train_size+i2];
			x1[j*train_size+i2] = tmp;
		}
		double tmp =y1[i1];
		y1[i1] = y1[i2];
		y1[i2] = tmp;
	}





	for(int j=0; j<ens_size; j++)
	{
		useyhalf1[j] = (double *)calloc(train_size/2, sizeof(double));
		useyhalf2[j] = (double *)calloc(train_size/2, sizeof(double));
	//	useytest[j] = (double *)calloc(test_size, sizeof(double));

		if (useyhalf1[j]==NULL || useyhalf2[j]==NULL) {	
			fprintf(stderr,"Allocation failed\n") ;
			return -1 ;
		}
	}

	for(int noftr=0; noftr<nftrs ; noftr++)
	{
		 
		for(int i=0;i<train_size/2; i++)
		{
			newx1half1[noftr*train_size/2+i] = x1[noftr*train_size+i];
			newx1half2[noftr*train_size/2+i] = x1[noftr*train_size+i+train_size/2];
			y1half1[i] = y1[i];
			y1half2[i] = y1[i+train_size/2];

		}

		/*for(int i=0; i<test_size; i++)
		{
			newx2[noftr*test_size+i] = x2[noftr*test_size+i];
		}*/
	}
 
	double *temp_w1 = (double *)( calloc(train_size/2, sizeof(double)));

	for(int tnum = 0; tnum<ens_size; tnum++)
	{	
		for(int i=0; i<train_size/2; i++)
		{
			unsigned int num=0;
#ifdef __linux
                        num = mt_rng();
#else
			rand_s(&num);
#endif
			if (num%3) temp_w1[i]=0; else temp_w1[i]=w1[i];
			 
		}


//		full_gbm_learn_info_t  full_info;
		if (get_gbm_predictor(newx1half1,y1half1,temp_w1,train_size/2,nftrs,ens_shrinkage,ens_bag_p,ens_take_all,ens_ntrees,ens_depth,min_obs_in_node,
			(full_info_ens->full_info_tree)+tnum)==-1) {							 
			fprintf(stderr,"GBM learning failed\n") ;
			return -1 ;
		}

		//		memcpy(&(full_info_ens->full_info_tree[tnum]), &full_info, sizeof(full_info));
		//		memcpy(save_trees[tnum], trees, sizeof(trees));
	
		// Predict
		if (gbm_predict(newx1half1,train_size/2,nftrs,ens_ntrees,(full_info_ens->full_info_tree)+tnum,useyhalf1[tnum])==-1) {
			fprintf(stderr,"GBM prediction failed\n") ;			 
			return -1 ;
		}

		// Predict
		if (gbm_predict(newx1half2,train_size/2,nftrs,ens_ntrees,(full_info_ens->full_info_tree)+tnum,useyhalf2[tnum])==-1) {
			fprintf(stderr,"GBM prediction failed\n") ;
			return -1 ;
		}
	}

	/*	for(int tn=0; tn<ens_size; tn++)
		{
			for(int i=0; i<test_size; i++)
			ipreds[i]+=useytest[tn][i];
	}*/
				 
	double *treeb = (double *) calloc(ens_size, sizeof(double));	
	int *tosort = (int *) calloc(train_size/2, sizeof(int));

	nfc = (double *) calloc(train_size/2, sizeof(double));
	double *ofc = (double *) calloc(train_size/2, sizeof(double));

	if (treeb==NULL || tosort==NULL || nfc==NULL || ofc==NULL) {
		fprintf(stderr,"Allocation failed\n") ;
		return -1 ;
	}

	for(int i=0; i<train_size/2; i++)
		tosort[i] = i;
	
	double cstep = 1.0;
	double besterr=0.9;

	int total0=0;
	int total1=0;
	for(int i1=train_size/2; i1<train_size; i1++)
		if (y1[i1]<0.1) total0++; else total1++;

	for(int it=0; it<150; it++)	
	{
		if (it>1) cstep*=0.98; 
		if (it==1) cstep*=0.1;

		for(int sign=-1; sign<2; sign+=2)
		{		 
			for(int j=0; j<ens_size; j++)
			{				 			
				if (treeb[j]+cstep*sign<0.1) continue;;
				if (treeb[j]+cstep*sign>10.0) continue;;
						 	
				for(int i=0; i<train_size/2; i++)
				{
					double ux = useyhalf2[j][i];
					ofc[i] = nfc[i];
					nfc[i]+=ux*cstep*sign;
					 
				}
				 
				qsort(tosort , train_size/2, sizeof(int), auccompare );

				int current0=total0;
				int sum=0;
				for(int i1=0; i1<train_size/2; i1++)
				{
					if (y1[tosort[i1]+train_size/2]==1) sum+=current0; else current0--;
				}

				double auc = (double)sum/(double)(total0*total1);


				sum=0;
				for(int i1=0; i1<(train_size/2)*0.1; i1++)
				{
					if (y1[tosort[i1]+train_size/2]==1) sum+=1;  

				}

				double sens90 = (double)sum/(double)(total1);
			
				double err = -auc;//-sens90;

				if ((err<besterr) || (it==0))
				{
					besterr=err;
					treeb[j]+=cstep*sign;						
					printf("Auc:%.4f Sens90:%.4f Err:%.4f\n", auc, sens90, err);
				} else
				{		 
					for(int i=0; i<train_size/2; i++)
							nfc[i] = ofc[i];
				}				
			}	
		}
	}

	double sumtreeb=0;
	for(int j=0; j<ens_size; j++)
	{
		fprintf(stderr,"Treeb %d  = %.6f\n", j, treeb[j]);
		sumtreeb += treeb[j];	
	}

	for(int j=0; j<ens_size; j++)
	{

		treeb[j]/=sumtreeb;
		full_info_ens->alpha[j] = treeb[j];
	}

	// Free
	free(temp_w1) ;
	free(newx1) ;
	free(newx1half1) ;
	free(newx1half2) ;
	free(y1half1) ;
	free(y1half2) ;

	for(int j=0; j<ens_size; j++)
	{
		free(useyhalf1[j]) ;
		free(useyhalf2[j]) ;
	}
	free(useyhalf1) ;
	free(useyhalf2) ;
	free(treeb) ;
	free(tosort) ;
	free(nfc) ;
	free(ofc) ;

	return 0 ;
}


int gbm_predict_ens(double *x, int nrows, int ncols, full_gbm_ens_learn_info_t *full_info_ens, double *preds) {

	int ens_size = full_info_ens->ens_size ;

	double **useytest = (double **) malloc(ens_size*sizeof(double *)) ;
	if (useytest==NULL) {
		fprintf(stderr,"Allocation failed\n") ;
		return -1 ;
	}

	for(int j=0; j<ens_size; j++)
	{
		useytest[j] = (double *) calloc(nrows, sizeof(double));

		 if (gbm_predict(x,nrows,ncols,(full_info_ens->full_info_tree)[j].ntrees,&(full_info_ens->full_info_tree[j]), useytest[j])==-1) {
				 fprintf(stderr,"GBM ENS prediction failed\n") ;
				  return -1 ;
		 }					
	}
	
	for(int i=0; i<nrows; i++)
	{
	  double fc=0;

	  for(int treeno =0; treeno<ens_size; treeno++)
	  {
			 
			 double ux = useytest[treeno][i];
	  		 fc +=ux*full_info_ens->alpha[treeno];
			 
	  }
		 
	  //fc+=avry;		  
	  preds[i] = fc;///sumtreeb;
		 
	  //printf("%.6f\n", ipreds[i]);
	}

	return 0 ;
}

// Clear
void clear_gbm_predictor_ens(full_gbm_ens_learn_info_t *full_info_ens) {

	free(full_info_ens->alpha) ;

	for (int i=0; i<full_info_ens->ens_size; i++)
		clear_gbm_info(full_info_ens->full_info_tree+i) ;

	free(full_info_ens->full_info_tree) ;
}

// Read and Write
int write_full_gbm_ens_info(full_gbm_ens_learn_info_t *ens_info, char *fname) {

	FILE *fp ;
	if ((fp = safe_fopen(fname, "w", false))==NULL) {
		fprintf(stderr,"Cannot open %s for reading\n",fname) ;
		return -1 ;
	}

	write_full_gbm_ens_info(ens_info,fp) ;
	fclose(fp) ;
	return 0 ;
}

void write_full_gbm_ens_info(full_gbm_ens_learn_info_t *ens_info, FILE *fp)  {

	fprintf(fp,"EnsSize = %d\n",ens_info->ens_size) ;
	for (int i=0; i<ens_info->ens_size; i++)
		fprintf(fp,"Alpha %d = %f\n",i,ens_info->alpha[i]) ;

	for (int i=0; i<ens_info->ens_size; i++) {
		fprintf(fp,"GBM %d\n",i) ;
		write_full_gbm_info((ens_info->full_info_tree)+i,fp) ;
	}
}

int read_full_gbm_ens_info(full_gbm_ens_learn_info_t *ens_info, char *fname) {

	FILE *fp ;
	if ((fp = safe_fopen(fname, "r", false))==NULL) {
		fprintf(stderr,"Cannot open %s for reading\n",fname) ;
		return -1 ;
	}

	int ens_size ;
	if (fscanf(fp,"EnsSize = %d\n",&ens_size) != 1) {
		fprintf(stderr,"Reading EnsSize failed\n") ;
		return -1 ;
	}
	ens_info->ens_size = ens_size ;

	ens_info->alpha = (double *) malloc(ens_size*sizeof(double)) ;
	ens_info->full_info_tree = (full_gbm_learn_info_t *) calloc (1,ens_size*sizeof(full_gbm_learn_info_t)) ;
	if (ens_info->alpha==NULL || ens_info->full_info_tree==NULL) {
		fprintf(stderr,"Allocation failed\n") ;
		return -1 ;
	}

	double alpha ;
	int ival ;

	for (int i=0; i<ens_size; i++) {
		if (fscanf(fp,"Alpha %d = %lf\n",&ival,&alpha)!=2 || ival != i) {
			fprintf(stderr,"Reading Alpha %d failed\n",i) ;
			return -1 ;
		}

		 
		ens_info->alpha[i] = alpha ;

	}

	int rc ;
	for (int i=0; i<ens_size; i++) {
		if (((rc=fscanf(fp,"GBM %d\n",&ival)) != 1) || ival != i) {

			fprintf(stderr,"Reading failed at %d. Rc = %d\n",i,rc) ;
			return -1 ;
		}

		if (read_full_gbm_info( (ens_info->full_info_tree)+i,fp) == -1) {
			fprintf(stderr,"Reading data for GBM %d failed\n",i) ;
			return -1 ;
		}

	
	}

 

	return 0 ;
}

// (De)Serialize gbm_ens_info
// Serialization
size_t get_gbm_ens_info_size(full_gbm_ens_learn_info_t *ens_info) {
	
	size_t ens_info_size = sizeof(int) ; // EnsSize
	ens_info_size += ens_info->ens_size * sizeof(double) ; //  Alphas 

	for (int i=0; i<ens_info->ens_size; i++)
		ens_info_size += get_gbm_info_size(ens_info->full_info_tree + i) ;
		
	return ens_info_size ;
}

int gbm_ens_serialize(full_gbm_ens_learn_info_t *ens_info, unsigned char *ens_data) {
	size_t idx = 0 ;
	memcpy(ens_data+idx,&(ens_info->ens_size),sizeof(int)) ; idx += sizeof(int) ;
	memcpy(ens_data+idx,ens_info->alpha,ens_info->ens_size * sizeof(double)) ; idx += ens_info->ens_size * sizeof(double) ;

	for (int i=0; i<ens_info->ens_size; i++)
		idx += gbm_serialize(ens_info->full_info_tree + i, ens_data + idx) ;

	return (int) idx ;
}

int gbm_ens_deserialize(unsigned char *ens_data, full_gbm_ens_learn_info_t *ens_info) {
	size_t idx = 0 ;
	memcpy(&(ens_info->ens_size),ens_data+idx,sizeof(int)) ; idx += sizeof(int) ;
	memcpy(ens_info->alpha,ens_data+idx, ens_info->ens_size * sizeof(double)); idx += ens_info->ens_size * sizeof(double) ;

	for (int i=0; i<ens_info->ens_size; i++)
		idx += gbm_deserialize(ens_data + idx, ens_info->full_info_tree + i) ;

	return (int) idx ;
}


