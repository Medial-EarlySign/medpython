#define _CRT_SECURE_NO_WARNINGS

#include "classifiers.h"
#include "medial_utilities/medial_utilities/globalRNG.h"

#define STREE_IT 1000
#define STREENUM 250
#define NPARTS 3
#define FOREST_PART 1.0 // 0.125
#define TOP_PART 0.01
#define MINV 12

double ***stree_saf,***stree_val ;
char ***stree_sign,***stree_leaf ;
int ***stree_x,***stree_next ;
double **mytesty,**streeb ;
int *treeones,*treeones2 ;

int allocate_nir_trees() {
	stree_saf = (double ***) malloc (NPARTS*sizeof(double **)) ;
	stree_sign = (char ***) malloc (NPARTS*sizeof(char **)) ;
	stree_x = (int ***) malloc(NPARTS*sizeof(int **)) ;
	stree_next = (int ***) malloc(NPARTS*sizeof(int **)) ;
	stree_leaf = (char ***) malloc (NPARTS*sizeof(char **)) ;
	stree_val = (double ***) malloc (NPARTS*sizeof(double **)) ;
	mytesty = (double **) malloc(NPARTS*sizeof(double *)) ;
	streeb = (double **) malloc(NPARTS*sizeof(double *)) ;

	if (stree_saf==NULL || stree_sign==NULL || stree_x==NULL || stree_next==NULL || stree_leaf==NULL || stree_val==NULL || mytesty==NULL || streeb==NULL) {
		fprintf(stderr,"Nir-Tree Allocation failed\n") ;
		return -1 ;
	}

	for (int i=0; i<NPARTS; i++) {
		stree_saf[i] = (double **) malloc(STREENUM*sizeof(double *)) ;
		stree_sign[i] = (char **) malloc(STREENUM*sizeof(char *)) ;
		stree_x[i] = (int **) malloc(STREENUM*sizeof(int *)) ;
		stree_next[i] = (int **) malloc(STREENUM*sizeof(int *)) ;
		stree_leaf[i] = (char **) malloc(STREENUM*sizeof(char *)) ;
		stree_val[i] = (double **) malloc(STREENUM*sizeof(double *)) ;
		streeb[i] = (double *) malloc(STREENUM*sizeof(double)) ;

		if (stree_saf[i]==NULL || stree_sign[i]==NULL || stree_x[i]==NULL || stree_next[i]==NULL || stree_leaf[i]==NULL || stree_val[i]==NULL || streeb[i]==NULL) {
			fprintf(stderr,"Nir-Tree Allocation failed\n") ;
			return -1 ;
		}

		for (int j=0; j<STREENUM; j++) {
			stree_saf[i][j] = (double *) malloc(STREENUM*sizeof(double)) ;
			stree_sign[i][j] = (char *) malloc(STREENUM*sizeof(char)) ;
			stree_x[i][j] = (int *) malloc(STREENUM*sizeof(int)) ;
			stree_next[i][j] = (int *) malloc(STREENUM*sizeof(int)) ;
			stree_leaf[i][j] = (char *) malloc(STREENUM*sizeof(char)) ;
			stree_val[i][j] = (double *) malloc(STREENUM*sizeof(double)) ;

			if (stree_saf[i][j]==NULL || stree_sign[i][j]==NULL || stree_x[i][j]==NULL || stree_next[i][j]==NULL || stree_leaf[i][j]==NULL || stree_val[i][j]==NULL) {
				fprintf(stderr,"Nir-Tree Allocation failed\n") ;
				return -1 ;
			}
		}
	}

	treeones = (int *) malloc(STREENUM*sizeof(int)) ;
	treeones2 = (int *) malloc(STREENUM*sizeof(int)) ;

	if (treeones==NULL || treeones2==NULL) {
		fprintf(stderr,"Nir-Tree Allocation failed\n") ;
		return -1 ;
	}

	return 0 ;
}

//
// x = predictors matrix (transposed)
// y = response vector
// (nrows,ncols) = dimensions of x
// training = flags for internal-training set
//
int make_trees(double *x, double *y, int nrows, int ncols, int *training, int gloop) {
	
	// Allocation
	int *dontusei = (int *) malloc(nrows*sizeof(int)) ;
	int *group = (int *) malloc(nrows*sizeof(int)) ;
	int *dumax = (int *) malloc(nrows*sizeof(int)) ;
	int *tgroup = (int *) malloc(nrows*sizeof(int)) ;
	int *workonlist = (int *) malloc(nrows*sizeof(int)) ;
	double *usexu = (double *) malloc(nrows*sizeof(double)) ;
	double *usey = (double *) malloc(nrows*sizeof(double)) ;
	int *du = (int *) malloc(ncols*sizeof(int)) ;
	double *streefc = (double *) malloc (STREENUM*nrows*sizeof(double)) ;
	double *fc = (double *) malloc(nrows*sizeof(double)) ;

	if (dontusei==NULL || du==NULL || dumax==NULL || group==NULL || tgroup==NULL || workonlist==NULL || usexu==NULL || usey==NULL || streefc==NULL || fc==NULL) {
		fprintf(stderr,"Allocation failed\n") ;
		return -1 ;
	}
	
	// Build trees		
	for(int itree=0; itree<STREENUM; itree++) {
		clock_t start = clock() ;
		fprintf(stderr,"Building tree %d/%d ...",itree+1,STREENUM) ;

		// Choose subset of samples
		for(int i=0; i<nrows; i++) {
					 
			if ((y[i]==0) && (globalRNG::rand()%5)) 
				dontusei[i] = 1; 
			else 
				dontusei[i] = 0;		
		}

		// Choose subset of features
		for(int j=0; j<ncols; j++) {
			if (globalRNG::rand()%3) 
				du[j] = 0; 
			else 
				du[j] = 1;				 		
		}

		for(int i=0; i<nrows;i++) {
			group[i]=0;
			if (!training[i]) group[i]=-1;
			if ( dontusei[i]) group[i]=-1;
		}

		// Initialize 
		for(int i=0; i<STREE_IT; i++)
			stree_leaf[gloop][itree][i]=0;

		int workon=0;
		int lastworkon=1;
		int workonsize=0;

		// Split nodes		
		for(int git=0; git<STREE_IT; git++) {
			workonsize=0;
			for(int i=0; i<nrows;i++) {
				if ((group[i]==workon))
					workonlist[workonsize++]=i;
			}
												
			int cnt=0;
			double sum=0;
			for(int i1=0; i1<workonsize; i1++) {	
				int i=workonlist[i1];
				cnt++;						 
				sum+=y[i];		
			}

			if (cnt==0) 
				break ;

			if ((cnt<(MINV+1)*2) || (sum==0) || (stree_leaf[gloop][itree][workon]==1)) {
				stree_leaf[gloop][itree][workon] = 1;
				workon++;
				continue;
			}

			stree_leaf[gloop][itree][workon] = 0;
			double maxp=0;
			double max_saf ;	
			int max_tryj,max_sign ;
			double cnt1,sum1,cnt2,sum2 ; 

			// Find optimal split
			for(int tryj=0; tryj<ncols; tryj++) {		
				if (du[tryj])
					continue;

				double maxv=-1e13;
				double minv=1e13;
				 for(int i1=0; i1<workonsize; i1++) {	
						int i=workonlist[i1];
						if (x[XIDX(tryj,i,nrows)]>maxv)
							maxv= x[XIDX(tryj,i,nrows)];

						if (x[XIDX(tryj,i,nrows)]<minv)
							minv= x[XIDX(tryj,i,nrows)];
				 }
				 
				 double start = minv+0.0001;
				 double end = maxv;
				 double step = (maxv-minv)/2.001    ;
				
				 if (step==0)
					 continue;

				 for(double saf=start; saf<=end; saf+= step) {
					 if (saf>end)
						 break;
			
					int sign=0;
					for(int sign=0; sign<2; sign++) {
						cnt1=sum1=cnt2=sum2=0 ;
					
						for(int i11=0; i11<workonsize; i11++) {
							int i=workonlist[i11];
							usexu[i] = x[XIDX(tryj,i,nrows)]; 
						}

						if (sign==0) {
						  for(int i1=0; i1<workonsize; i1++) {					
							  int i=workonlist[i1];								  
							  if (usexu[i] >saf) {
								  cnt1++;
								  sum1+=y[i];

							  } else {
								  cnt2++;
								  sum2+=y[i];							
							  }
						  }
						} else {
							for(int i1=0; i1<workonsize; i1++) {	
								int i=workonlist[i1]; 
								if (usexu[i]<saf) {
									cnt1++;
									sum1+=y[i] ;
								   
								} else {
									cnt2++;
									sum2+=y[i] ;								
								}
							}   
						}	
					
						double parta=0;
						double partb=0;

						sum1+=(globalRNG::rand()%1000)/1000000000.0;
						sum2+=(globalRNG::rand()%1000)/1000000000.0;
						if ((cnt1>MINV)  && (cnt2>MINV)  /* &&  (cnt1>cnt*0.1)  && (cnt2>cnt*0.1) */ ) {
							parta = sum1/cnt1;
							partb = sum2/cnt2;
						} else
							continue;				
								
							
						double val  = (parta -partb)+(double)(cnt1-cnt2+1)/1000000.0;
								 
						if (val>maxp) {						
							maxp=val;
							max_sign = sign ;
							max_saf = saf ;
							max_tryj = tryj ;
						}
					}
				 }
			}
	
			if (maxp > 0) {
				cnt1=sum1=cnt2=sum2=0 ;
				
				for(int i11=0; i11<workonsize; i11++) {
					int i=workonlist[i11];
					usexu[i] = x[XIDX(max_tryj,i,nrows)]; 
				}

				if (max_sign==0) {
				  for(int i1=0; i1<workonsize; i1++) {					
					  int i=workonlist[i1];								  
					  if (usexu[i] > max_saf) {
						  cnt1++;
						  sum1+=y[i];

					  } else {
						  cnt2++;
						  sum2+=y[i];							
					  }
				  }
				} else {
					for(int i1=0; i1<workonsize; i1++) {	
						int i=workonlist[i1]; 
						if (usexu[i]<max_saf) {
							cnt1++;
							sum1+=y[i] ;
						   
						} else {
							cnt2++;
							sum2+=y[i] ;								
						}
					}   
				}

				sum1+=(globalRNG::rand()%1000)/1000000000.0;
				sum2+=(globalRNG::rand()%1000)/1000000000.0;
              
				double parta = sum1/cnt1;
				double partb = sum2/cnt2;
						
				for(int i=0; i<nrows ;i++) {										
					stree_leaf[gloop][itree][lastworkon]=0;				
					stree_leaf[gloop][itree][lastworkon+1]=0;
					
					stree_val[gloop][itree][lastworkon] = parta;					
					stree_val[gloop][itree][lastworkon+1] = partb;

					if ((cnt1<(MINV+1)*2))										
						stree_leaf[gloop][itree][lastworkon]=1;
										 
					if ((cnt2<(MINV+1)*2)   || (sum1<=1.001) )
						stree_leaf[gloop][itree][lastworkon+1]=1;

					if (max_sign==0) {
						 if (x[XIDX(max_tryj,i,nrows)]>max_saf)
							tgroup[i]=lastworkon;
						 else
						   tgroup[i]=lastworkon+1;
					} else {
						 if (x[XIDX(max_tryj,i,nrows)]<max_saf)
							tgroup[i]=lastworkon;
						 else
							 tgroup[i]=lastworkon+1;
					}
				}
		
				stree_sign[gloop][itree][git] = max_sign;		
				stree_saf[gloop][itree][git] = max_saf;
				stree_x[gloop][itree][git]=max_tryj;
				stree_next[gloop][itree][workon] = lastworkon;				

			} else {
	
				stree_leaf[gloop][itree][workon] = 1;
				workon++;
				continue;		
			}

					
			for(int i=0; i<nrows ;i++) {			
				if (group[i]==workon)
					group[i]=tgroup[i];	
			}

			lastworkon+=2;
			if (lastworkon >= STREE_IT) {
				fprintf(stderr,"Tree Exceeded nodes limit %d\n",STREE_IT) ;
				return -1 ;
			}

			workon++;
		}

		fprintf(stderr," (%d) took %f seconds\n",workon, (0.0 + clock() - start)/CLOCKS_PER_SEC) ;
	}
				 
	// Calc weights of trees by least square on internal test-set
	int sumtfc=0;
			
	// Create matrix
	for (int i=0; i<nrows ;i++) {			
		for(int itree=0; itree<STREENUM; itree++) {
			streefc[XIDX(itree,i,nrows)] = 0 ;
			int workon=0;

			for(;;) {
				if (stree_leaf[gloop][itree][workon]) {
					streefc[XIDX(itree,i,nrows)]= stree_val[gloop][itree][workon] ;						
					break;		
				}

				int tryj= stree_x[gloop][itree][workon];
				int sign = stree_sign[gloop][itree][workon];
				double saf = stree_saf[gloop][itree][workon];

				if (sign==0) {
					 if (x[XIDX(tryj,i,nrows)]>saf)
						 workon=stree_next[gloop][itree][workon]; 
					 else 
						 workon=stree_next[gloop][itree][workon]+1 ;
				} else {
					 if (x[XIDX(tryj,i,nrows)]<saf) 
						 workon=stree_next[gloop][itree][workon];
					 else
						 workon=stree_next[gloop][itree][workon]+1;
				}		
			}					
		}		
	}
	
	// Find subset of trees by performance
	for(int itree=0; itree<STREENUM; itree++) {
 
		for(int j=0; j<nrows; j++)
			dumax[j] = training[j];

		int onesctr=0;

		for(int maxones=0; maxones<TOP_PART * nrows/(NPARTS*1.0) ; maxones++) {
			double max=-0.1;
			int ptr=-1;

			for(int i=0; i<nrows; i++) {	
				if (dumax[i]) continue;
				if (streefc[XIDX(itree,i,nrows)]>max) {
					max = streefc[XIDX(itree,i,nrows)];
					ptr = i;
				}
			}

			dumax[ptr] = 1;
			if (y[ptr]>0) 
				onesctr++;
		}

		treeones[itree] = onesctr;
		treeones2[itree] = onesctr ;
			
	}
			
	int maxusetree=0;
	for(int l1=0; l1<STREENUM; l1++) {
		for(int l2=0; l2<l1; l2++) {
			if (treeones[l1]>treeones[l2]) {
				int temp = treeones[l1];
				treeones[l1] = treeones[l2];
				treeones[l2]= temp;
			}
		}
	}

	maxusetree= treeones[(int)(STREENUM*FOREST_PART)];
	fprintf(stderr,"FIRST:!!!! %d    MAXUSE:%d\n", treeones[0], maxusetree);
			
	memcpy( usey, y, nrows*sizeof(double));
	for(int wn =0; wn<STREENUM; wn++)
		streeb[gloop][wn] = 0;

	// Linear regression
	for(int it=0; it<100; it++){
		double err=0;
		for(int i=0; i<nrows ;i++) {
			if (training[i]) 
				continue;
			err+=usey[i]*usey[i];
		}

		for(int wn =0; wn<STREENUM; wn++) {
			if (treeones2[wn]<maxusetree)
				continue;

			double sumxx=0;
			double sumxy=0;
			for(int i=0; i<nrows ;i++) {
				if (training[i])
					continue;
		
				double x= streefc[XIDX(wn,i,nrows)] ;
				double y=usey[i];
				sumxx+=x*x;
				sumxy+=x*y;
			}

						
			double alpha = sumxy/sumxx;
			double oldb=streeb[gloop][wn];
			streeb[gloop][wn]+=alpha;

			if (streeb[gloop][wn]<0) 
				streeb[gloop][wn]=0;

			streeb[gloop][wn] *= 0.5 ;
			alpha = streeb[gloop][wn]-oldb;

			for(int i=0; i<nrows ;i++) {
				if (training[i])
					continue;
				double x= streefc[XIDX(wn,i,nrows)];
				usey[i]-=x*alpha;
			}
		}
	}


	double avrfc=0;
	double avry=0;
	double cnt1=0;
	for(int i=0; i<nrows ;i++) {
		if (training[i])
			continue;
		fc[i] = y[i]-usey[i];

		cnt1++;
		avrfc+=fc[i];
		avry+=y[i];
	}
	avry/=cnt1;
	avrfc/=cnt1;


	double sumxx=0;
	double sumyy=0;
	double sumxy=0;
	for(int i=0; i<nrows ;i++)
	{
		if (training[i]) 
			continue;
		double xval = fc[i]-avrfc;
		double yval =  y[i] - avry;
		sumxx+=xval*xval;
		sumxy+=xval*yval;
		sumyy+=yval*yval;
	}

	double correl = sumxy/sqrt(sumxx*sumyy);
	fprintf(stderr,"TOTAL CORREL: %.4f\n", correl);


	// Clearn
	free(dontusei); free(group) ; free(dumax) ;
	free(tgroup); free(workonlist); free(usexu); free(usey) ;
	free(du); free(streefc); free(fc) ;

	return 0 ;
}

int runontest(double *x, int nrows, int ncols, int gloop) {
	
	// Allocation
	double *streefc = (double *) malloc (STREENUM*nrows*sizeof(double)) ;
	mytesty[gloop] = (double *) malloc(nrows*sizeof(double)) ;

	if (streefc==NULL || mytesty[gloop]==NULL) {
		fprintf(stderr,"Allocation failed\n") ;
		return -1 ;
	}
	int sumtfc=0;
	
	for(int i=0; i<nrows ;i++) {		
		for(int itree=0; itree<STREENUM; itree++) {

			streefc[XIDX(itree,i,nrows)]=0;
			int workon=0;

			for(;;) {			 
				if (stree_leaf[gloop][itree][workon]) {
					streefc[XIDX(itree,i,nrows)]= stree_val[gloop][itree][workon] ;
					break;
				}

				int tryj= stree_x[gloop][itree][workon];
				int sign = stree_sign[gloop][itree][workon];
				double saf = stree_saf[gloop][itree][workon];

				if (sign==0) {
					 if (x[XIDX(tryj,i,nrows)]>saf) 
						 workon=stree_next[gloop][itree][workon];
					 else 
						 workon=stree_next[gloop][itree][workon]+1;
				} else {
					 if (x[XIDX(tryj,i,nrows)]<saf) 
						 workon=stree_next[gloop][itree][workon]; 
					 else 
						 workon=stree_next[gloop][itree][workon]+1;
				}
			}				 
		}
	}

	for(int i=0; i<nrows; i++) {
		mytesty[gloop][i] = 0;
		for(int l=0; l<STREENUM; l++)
			mytesty[gloop][i]+=streefc[XIDX(l,i,nrows)]*streeb[gloop][l];
	}

	// Clean
	free(streefc) ;	

	return 0 ;
}

void finalresult(double *x, int nrows, int ncols, double *predictions) {

		for(int i=0; i<nrows ;i++) {
			for(int j1=0; j1<NPARTS; j1++) {
				for(int j2=0; j2<j1; j2++) {
					if (mytesty[j1][i]>mytesty[j2][i]) {
						double temp  =mytesty[j1][i];
						mytesty[j1][i] = mytesty[j2][i];
						mytesty[j2][i] = temp;
					}
				}
			}

			predictions[i] = mytesty[1][i]+mytesty[0][i]+mytesty[2][i];
		}
}

int print_nir_trees(char *file_name) {
	
	FILE *fp = safe_fopen(file_name, "w", false) ;
	if (fp==NULL) {
		fprintf(stderr,"Cannot open %s for writing\n",file_name) ;
		return -1 ;
	}

	fprintf(fp,"%d %d %d\n",NPARTS,STREENUM,STREE_IT) ;
	for (int gloop=0; gloop<NPARTS; gloop++) {
		for (int itree=0; itree<STREENUM; itree++) {
			fprintf(fp,"%lf\n",streeb[gloop][itree]) ;

			for (int i=0; i<STREE_IT; i++) {
				fprintf(fp,"%d %d %d %d %lf %lf\n",stree_leaf[gloop][itree][i],stree_x[gloop][itree][i],stree_sign[gloop][itree][i],stree_next[gloop][itree][i],
					stree_val[gloop][itree][i],stree_saf[gloop][itree][i]) ;
			}
		}
	}

	return 0 ;
}

int read_nir_trees(char *file_name) {
	
	if (allocate_nir_trees() == -1)
		return -1 ;

	FILE *fp = safe_fopen(file_name, "r", false) ;
	if (fp==NULL) {
		fprintf(stderr,"Cannot open %s for writing\n",file_name) ;
		return -1 ;
	}

	int nparts,streenum,stree_it ;
	fscanf(fp,"%d %d %d\n",&nparts,&streenum,&stree_it) ;
	if (nparts != NPARTS || streenum != STREENUM || stree_it != STREE_IT) {
		fprintf(stderr,"Wrong parameters of nir trees in file : %d %d %d\n",nparts,streenum,stree_it) ;
		return -1 ;
	}
	for (int gloop=0; gloop<NPARTS; gloop++) {
		for (int itree=0; itree<STREENUM; itree++) {
			if (fscanf(fp,"%lf\n",&(streeb[gloop][itree])) != 1) {
				fprintf(stderr,"Problems reading b for %d %d\n",gloop,itree) ;
				return -1 ;
			}
//			streeb[gloop][itree] = 1 ;

			for (int i=0; i<STREE_IT; i++) {
				if (fscanf(fp,"%hhd %d %hhd %d %lf %lf\n",&(stree_leaf[gloop][itree][i]),&(stree_x[gloop][itree][i]),&(stree_sign[gloop][itree][i]),&(stree_next[gloop][itree][i]),
					&(stree_val[gloop][itree][i]),&(stree_saf[gloop][itree][i])) != 6) {
						fprintf(stderr,"Problem reading data for %d %d %d\n",gloop,itree,i) ;
						return -1 ;
				}
			}
		}
	}

	return 0 ;
}

int learn_nir_trees(double *x, double *y, int nrows, int ncols, char *file_name) {

	if (allocate_nir_trees() == -1)
		return -1 ;

	int *training = (int *) malloc(nrows*sizeof(int)) ;
	if (training == NULL) {
		fprintf(stderr,"Allocation failed\n") ;
		return -1 ;
	}

	for(int gloop=0; gloop<NPARTS; gloop++) {
		fprintf(stderr,"Loop %d/%d\n",gloop+1,NPARTS) ;

		for(int i=0; i<nrows; i++) {
        if (i%3==gloop)
			training[i]=0; 
		else 
			training[i] = 1;
		}

		make_trees(x,y,nrows,ncols,training,gloop);
	}

	return print_nir_trees(file_name) ;
}

int nir_trees_predict(double *x, double *preds, int nrows, int ncols, char *file_name) {

	if (read_nir_trees(file_name) == -1) 
		return -1 ;

	for (int gloop=0; gloop<NPARTS; gloop++) {
		fprintf(stderr,"Predicting loop %d/%d\n",gloop+1,NPARTS) ;
		if (runontest(x,nrows,ncols,gloop) == -1)
			return -1 ;
	}

	finalresult(x,nrows,ncols,preds);

	return 0;
}
