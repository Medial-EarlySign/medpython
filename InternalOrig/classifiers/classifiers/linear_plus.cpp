// Linear+ predictor (Nir)
#define _CRT_SECURE_NO_WARNINGS

#include "classifiers.h"
#include "medial_utilities/medial_utilities/globalRNG.h"

double *nfc_lp;
int min_auc_ftr =-1;
double exparray[100000];

inline int lp_auccompare( const void *arg1, const void *arg2 )
{
	return((int)(100*(nfc_lp[*(int *)arg2]-nfc_lp[*(int *)arg1])+0.5));
	//if (nfc_lp[*(int *)arg1]<nfc_lp[*(int *)arg2]) return(1); 
	//if (nfc_lp[*(int *)arg1]>nfc_lp[*(int *)arg2]) return(-1);
	return(0);
}

inline void rand_func(unsigned int *num)
{
	int rand1 = globalRNG::rand();
	int rand2 = globalRNG::rand();
	*num = (rand1 <<15)^rand2 ;
}

int get_linear_plus_predictions(double *x1, double *y1, int nrows1, double *x2, int nrows2, int ncols, double *preds) {

	// Learn
	full_linear_plus_info_t linear_plus_model ;

	get_linear_plus(x1,y1,nrows1,ncols,&linear_plus_model) ;

	// Predict
	if (predict_linear_plus(x2,nrows2,ncols,&linear_plus_model,preds)==-1) {
		fprintf(stderr,"Prediction with Linear-Plus failed\n") ;
		return -1 ;
	}

	clear_linear_plus(&linear_plus_model) ;

	return 0 ;
}

int learn_linear_plus_predictor(double *x, double *y, int nrows, int ncols, unsigned char **model) {
	
	// Learn
	full_linear_plus_info_t linear_plus_model ;

	get_linear_plus(x,y,nrows,ncols,&linear_plus_model) ;

	// Allocate
	size_t linear_plus_size = get_linear_plus_size(&linear_plus_model) ;
	if (((*model) = (unsigned char *) malloc(linear_plus_size))==NULL) {
		fprintf(stderr,"Allocation of LinearPlus model failed\n") ;
		return -1 ;
	}

	// Serialize
	if (linear_plus_serialize(&linear_plus_model, *model) == -1) {
		fprintf(stderr,"Serialization of GBM model failed\n") ;
		return -1 ;
	}

	 // Clear
	 clear_linear_plus(&linear_plus_model) ;

	 return (int) linear_plus_size ;
}

void get_linear_plus(double *x1, double *y1, int train_size, int nftrs, full_linear_plus_info_t *full_info_lin_plus )
{

	int cntpos=0;
	for(int i=0; i<train_size; i++)
		if (y1[i]>0.1) cntpos++;
	fprintf(stderr,"In get linear plus,  negative: %d  positive: %d\n", train_size-cntpos, cntpos);


	double *usevec  = (double *) calloc(sizeof(double), train_size);
	double *usey = (double *) calloc(sizeof(double), train_size);
	double *usey2 = (double *) calloc(sizeof(double), train_size);


	double *correl = (double *) calloc(sizeof(double), nftrs);





	double *newx1 = (double *) calloc(sizeof(double), train_size*nftrs);
	double *newy1 =  (double *) calloc(sizeof(double), train_size);
	int *tosort = (int *) calloc(train_size, sizeof(int));

	int *whatplace = (int *) calloc(sizeof(int), train_size);
	int *bestwhatplace = (int *) calloc(sizeof(int), train_size);


	 



	time_t start, end;

	time(&start);
	fprintf(stderr, "calc exp");

	for(double num = -500; num<500; num+=0.01)
	{
		exparray[50000+(int)(num*100)] =  exp(num);
	}

		 



	fprintf(stderr, "finished to calc exp");


	full_info_lin_plus->nftrs = nftrs;
	full_info_lin_plus->avrgs = (double *) calloc(nftrs, sizeof(double));
	full_info_lin_plus->stdevs = (double *) calloc(nftrs, sizeof(double));

	full_linear_plus_info_t best_full_info_lin_plus;

	full_info_lin_plus->linesno = 0;
	full_info_lin_plus->lines = (full_linear_plus_info_line_t *) calloc(LP_MAXLINES, sizeof(full_linear_plus_info_line_t));
	best_full_info_lin_plus.lines = (full_linear_plus_info_line_t *) calloc(LP_MAXLINES, sizeof(full_linear_plus_info_line_t));

	int bestlineno=0;


	nfc_lp = (double *) calloc(train_size, sizeof(double));
	double *ofc = (double *) calloc(train_size, sizeof(double));
	double *savenfc_lp = (double *) calloc(train_size, sizeof(double));

	for(int j=0; j<nftrs; j++)
	{
		double avrx=0;
		double cntx=0;
		double avrxx=0;

		for(int i=0; i<train_size; i++)
		{
			double val= x1[j*train_size+i];

			avrx+=val;
			avrxx+=val*val;
			cntx+=1.0;



			newx1[j*train_size+i] = val;
		 

		}

		avrx/=cntx;
		avrxx = sqrt(avrxx/cntx-avrx*avrx);

		 
		full_info_lin_plus->avrgs[j] = avrx;
		full_info_lin_plus->stdevs[j] = avrxx;




		printf("Avr: %.6f   Stdev: %.6f\n", avrx,avrxx);


		for(int i=0; i<train_size; i++)
		{
			newx1[j*train_size+i]-=avrx;
			newx1[j*train_size+i]/=avrxx;

		}


	}

	double avry=0;
	for(int i=0; i<train_size; i++)
		avry+= y1[i];

	avry/=train_size;


	int total0=0;
	int total1=0;
	
	for(int i=0; i<train_size; i++)
	{
		usey[i] = y1[i];
		usey2[i] = usey[i];
		tosort[i] = i;
		if (y1[i]<0.1) total0++; else total1++;

	}

	for(int j=0; j<nftrs; j++)
	{
		for(int i=0; i<train_size; i++)
			nfc_lp[i] = newx1[j*train_size+i];


		qsort(tosort , train_size , sizeof(int), lp_auccompare );

		int current0=total0;
		int sum=0;
		for(int i1=0; i1<train_size ; i1++)
		{
			if (y1[tosort[i1]]==1) sum+=current0; else current0--;
		}

		double auc = (double)sum/(double)(total0*total1);

		 
		correl[j] = auc;
		 
		printf("-----> %d %.6f  %d %d %de\n", j, correl[j], sum, total0, total1);
	 

		correl[j] = fabs(0.5-auc);
	}


	double min = 1000;
	
	for(int j=0; j<nftrs; j++)
	{
		if (correl[j]<min)
		{
			min = correl[j];
			min_auc_ftr = j;
		}
	}

	printf("Min auc ftr: %d\n", min_auc_ftr);



	int  BESTO=10;
	int CHOOSE_FROM = 35;
	if (CHOOSE_FROM>nftrs) CHOOSE_FROM = nftrs;
	if (BESTO>nftrs) BESTO = nftrs/2;


	memset(nfc_lp, 0, train_size*sizeof(double));
	int *bestar = (int *) calloc(nftrs,sizeof(int));
	 
	for(int it=0; it<nftrs; it++)
	{
		double bcorrel=0;
		int bj=-1;
		for(int j=0; j<nftrs; j++)
		{
			if (fabs(correl[j])>bcorrel)
			{
				bcorrel = fabs(correl[j]);
				bj = j;
			}
		}
		bestar[it] = bj;
		printf("%d)  ftr: %d  Correl: %.6f\n", it, bj, bcorrel);
		correl[bj] = -0.001;
	}


	

	 



	 

    

	for(int i=0; i<train_size; i++)
	{
		 
			nfc_lp[i] = 0 ;
			ofc[i] = 0 ;

			
		 

	}


	 
		 




	double bestpauc[AUCPARTS];
	double bestsens90[AUCPARTS];


 
	 
	double bestrauc = 0;
	memset(bestpauc,0, sizeof(bestpauc));
	memset(bestsens90,0, sizeof(bestsens90));

	 
 

 


/*	
#define ITERATION_FOR_ONE_CHOICE 3000
#define NUMBER_OF_FEATURES_PER_RUN  20
#define ENS_SIZE 1
*/


	double RIDGE_FACTOR = 1.0;



	int lastrecordit = 0;
	 
	int lastzerofeature = 0;
	int totalsize = ITERATION_FOR_ONE_CHOICE*NUMBER_OF_FEATURES_PER_RUN*ENS_SIZE;



	for(int it=0; it<ITERATION_FOR_ONE_CHOICE*NUMBER_OF_FEATURES_PER_RUN*ENS_SIZE; it++)
	{
		if (it%1000==0) printf("[IT=%d]\n", it);
		if (it%ITERATION_FOR_ONE_CHOICE==(ITERATION_FOR_ONE_CHOICE-1))
			{

				int lineno = full_info_lin_plus->linesno;

				int rand_ftr=full_info_lin_plus->lines[lineno]. l_rand_ftr;
		 
				int rand_ftrb=full_info_lin_plus->lines[lineno].l_rand_ftrb;
				int rand_ftrc=full_info_lin_plus->lines[lineno].l_rand_ftrc;
				int rand_ftrd=full_info_lin_plus->lines[lineno].l_rand_ftrd;
		 
				
				full_info_lin_plus->lines[lineno].l_rand_size *= RIDGE_FACTOR;		 


				printf("now adding  %d %d %d %d\n", rand_ftr, rand_ftrb, rand_ftrc, rand_ftrd );
				full_info_lin_plus->linesno = (full_info_lin_plus->linesno)+1;

				RIDGE_FACTOR-=0.01;
				if (RIDGE_FACTOR<0.6) RIDGE_FACTOR=0.6;

				bestlineno=0;
				int reshufflecounter=0;
				double minstdx = 1e20;
				double stdxsaf = 0.0005;

				// Do the calcs from beginning
reshuffle:
				memset(bestpauc,0, sizeof(bestpauc));
				memset(bestsens90,0, sizeof(bestsens90));
			 
				bestrauc = 0;
				 


				for(int it2=0; it2<5; it2++)
				for(int it3=0; it3<train_size; it3++)
				{
					int i1 = it3;
					int i2 = globalRNG::rand30()%train_size;
					for(int j=0; j<nftrs; j++)
					{
						double temp = newx1[j*train_size+i1];
						newx1[j*train_size+i1]=newx1[j*train_size+i2];
						newx1[j*train_size+i2] = temp;
					}
					double temp = y1[i1];
					y1[i1]= y1[i2];
					y1[i2] = temp;
				}




				for(int i=0; i<train_size; i++)
				{
					nfc_lp[i] = 0;
				}
				
				for(int j1=lastzerofeature; j1<full_info_lin_plus->linesno; j1++)
				{
					int rand_act=full_info_lin_plus->lines[j1].l_rand_act;  
					int rand_dir=full_info_lin_plus->lines[j1].l_rand_dir;
					int rand_ftr=full_info_lin_plus->lines[j1].l_rand_ftr;
					int rand_ftrb=full_info_lin_plus->lines[j1].l_rand_ftrb;
					int rand_ftrc=full_info_lin_plus->lines[j1].l_rand_ftrc;
					int rand_ftrd=full_info_lin_plus->lines[j1].l_rand_ftrd;
					int rand_ab=full_info_lin_plus->lines[j1].l_rand_ab;
					int rand_abc=full_info_lin_plus->lines[j1].l_rand_abc;
					int rand_abcd=full_info_lin_plus->lines[j1].l_rand_abcd;
					double first_saf=full_info_lin_plus->lines[j1].l_first_saf;
					double second_saf=full_info_lin_plus->lines[j1].l_second_saf;
					double first_safb=full_info_lin_plus->lines[j1].l_first_safb;
					double second_safb=full_info_lin_plus->lines[j1].l_second_safb;
					double first_safc=full_info_lin_plus->lines[j1].l_first_safc;
					double second_safc=full_info_lin_plus->lines[j1].l_second_safc;
					double first_safd=full_info_lin_plus->lines[j1].l_first_safd;
					double second_safd=full_info_lin_plus->lines[j1].l_second_safd;
					double rand_size=full_info_lin_plus->lines[j1].l_rand_size;
				

				 

					 
					for(int i=0; i<train_size; i++)
					{
						double vala = newx1[rand_ftr*train_size+i];
						double valb = newx1[rand_ftrb*train_size+i];
						double valc = newx1[rand_ftrc*train_size+i];
						double vald = newx1[rand_ftrd*train_size+i];

					 

						double useit=0;
						double useitb=0;
						double useitc=0;
						double useitd=0;

						 
						if (vala>second_saf) useit = 1; else
						if (vala<first_saf)  useit = 0; else
						useit = (vala-first_saf)/(second_saf-first_saf);

						
						if (valb>second_safb) useitb = 1; else
						if (valb<first_safb)  useitb = 0; else
						useitb = (valb-first_safb)/(second_safb-first_safb);
						if (rand_ab) useitb = 1- useitb;


						if (valc>second_safc) useitc = 1; else
						if (valc<first_safc)  useitc = 0;   else
						useitc = (valc-first_safc)/(second_safc-first_safc);
						if (rand_abc) useitc = 1- useitc;

						if (vald>second_safd) useitd = 1; else
						if (vald<first_safd)  useitd = 0;   else
						useitd = (vald-first_safd)/(second_safd-first_safd);
						if (rand_abcd) useitd = 1- useitd;


						 useit*=useitb; 
						 useit*=useitc; 
						 useit*=useitd; 
						 useit *=rand_size;


					
						nfc_lp[i]+=useit;					 	 
						savenfc_lp[i]=nfc_lp[i];
						ofc[i]=nfc_lp[i];
					}


				}
				
				
				{
					double total0[AUCPARTS];
					double total1[AUCPARTS];
					memset(total0, 0, sizeof(total0));
					memset(total1, 0, sizeof(total0));

					for(int j3=0; j3<AUCPARTS; j3++)
					{

						for(int i=j3*train_size/AUCPARTS; i<(j3+1)*train_size/AUCPARTS; i++)
						{					
							tosort[i] = i;
							if (y1[i]<0.1) total0[j3]++; else total1[j3]++;

						}

						qsort(tosort+j3*train_size/AUCPARTS , train_size/AUCPARTS , sizeof(int), lp_auccompare );

						int current0=(int)total0[j3];
						int sum=0;
						for(int i=j3*train_size/AUCPARTS; i<(j3+1)*train_size/AUCPARTS; i++)
						{
							if (y1[tosort[i]]==1) sum+=current0; else current0--;
						}


						bestpauc[j3] = (double)sum/(double)(total0[j3]*total1[j3]) * 1.0001;

						double part1=0;
						for(int i=j3*train_size/AUCPARTS; i<(j3+0.1)*train_size/AUCPARTS; i++)
						{
							if (y1[tosort[i]]==1) part1+=1.0;
						}

						bestsens90[j3]= part1/total1[j3]-0.001;
						 
						 
					}


					double avg=0;
					double stdx=0;
					for(int j3=0; j3<AUCPARTS; j3++)
					{
						avg+=bestpauc[j3];
					}
					avg/=AUCPARTS;


					for(int j3=0; j3<AUCPARTS; j3++)
					{
						stdx+=(bestpauc[j3]-avg)*(bestpauc[j3]-avg);
					}
					

					if (stdx<minstdx)
					{
						minstdx=stdx;
					}

					stdx = stdx/AUCPARTS;
					stdx = sqrt(stdx);

					printf("Stdx: %.4f (%.4f)  Auc: ", stdx, stdxsaf);

					if (reshufflecounter%10==9) stdxsaf+=0.0002;


					for(int j3=0; j3<AUCPARTS; j3++)
					{
						printf("%.4f ", bestpauc[j3]);
					}

					if (0)
					//if (stdx>stdxsaf)
					{
						printf("... Reshuffle\n");
						reshufflecounter++;
						goto reshuffle;
					}


					printf("Use it!!!!\n");
				}


			 
					

			 

				if (it%(ITERATION_FOR_ONE_CHOICE*NUMBER_OF_FEATURES_PER_RUN)==(ITERATION_FOR_ONE_CHOICE*NUMBER_OF_FEATURES_PER_RUN-1))
				{
					printf("Zero all..............\n");
					for(int i=0; i<train_size; i++)
					{
						 
							nfc_lp[i]=0 ;
							savenfc_lp[i]=0 ;
							ofc[i] = 0 ;
						 
					}
					 RIDGE_FACTOR = 1.0;
					bestrauc = 0;
					memset(bestpauc,0, sizeof(bestpauc));
					memset(bestsens90,0, sizeof(bestsens90));
					lastzerofeature = full_info_lin_plus->linesno;

					


				}



			}

 
		unsigned int num;
		rand_func(&num);
		int rand_ftr = bestar[num%CHOOSE_FROM];

		if (globalRNG::rand()%3==0)
		{
			rand_ftr = bestar[num%BESTO];
		}


		rand_func(&num);
		int rand_ftrb =bestar[num%CHOOSE_FROM];


		if (globalRNG::rand()%3==0)
		{
			rand_ftrb = bestar[num%BESTO];
		}


		rand_func(&num);
		int rand_ftrc =bestar[num%CHOOSE_FROM];


		if (globalRNG::rand()%3==0)
		{
			rand_ftrc = bestar[num%BESTO];
		}
 

		rand_func(&num);
		int rand_ftrd =bestar[num%CHOOSE_FROM];


		if (globalRNG::rand()%3==0)
		{
			rand_ftrd = bestar[num%BESTO];
		}
 


	


	





		rand_func(&num);
		int rand_ab = num%2;
		rand_func(&num);
		int rand_abc = num%2;
		rand_func(&num);
		int rand_abcd = num%2;


		

	

// Parameter 
		 
		 
 


 
		rand_func(&num);
		double first_saf  = ((double)(num%8000)/1000.0)-4.0;		
		rand_func(&num);
		double second_saf = ((double)(num%8000)/1000.0);				
		second_saf+=first_saf+0.01;

		
		rand_func(&num);
		double first_safb  = ((double)(num%8000)/1000.0)-4.0;		
		rand_func(&num);
		double second_safb = ((double)(num%8000)/1000.0);				
		second_safb+=first_safb+0.01;

		
		rand_func(&num);
		double first_safc  = ((double)(num%8000)/1000.0)-4.0;		
		rand_func(&num);
		double second_safc = ((double)(num%8000)/1000.0);				
		second_safc+=first_safc+0.01;

		rand_func(&num);
		double first_safd  = ((double)(num%8000)/1000.0)-4.0;		
		rand_func(&num);
		double second_safd = ((double)(num%8000)/1000.0);				
		second_safd+=first_safd+0.01;
		


		 


	 
		 
		  

	 
		rand_func(&num);
		double rand_size = (double)(num%25000) /5000.0;

		 
		if (globalRNG::rand()%2) rand_size=-rand_size;

		



		if (globalRNG::rand()%5==0)
		{
			first_safc=-100;
			second_safc=-99;
			rand_abc=0;
		}


		if (globalRNG::rand()%5==0)
		{
			first_safb=-100;
			second_safb=-99;
			rand_ab=0;
		}


		 


		 
		
		 




		if (((it%ITERATION_FOR_ONE_CHOICE)>(ITERATION_FOR_ONE_CHOICE*0.3)) && (globalRNG::rand()%2))
		{
			 int lineno = full_info_lin_plus->linesno;
			 rand_ftr=full_info_lin_plus->lines[lineno].l_rand_ftr;
			 rand_ftrb=full_info_lin_plus->lines[lineno].l_rand_ftrb;
			 rand_ftrc=full_info_lin_plus->lines[lineno].l_rand_ftrc;
			 rand_ftrd=full_info_lin_plus->lines[lineno].l_rand_ftrd;
			 rand_ab=full_info_lin_plus->lines[lineno].l_rand_ab;
			 rand_abc=full_info_lin_plus->lines[lineno].l_rand_abc;
			 rand_abcd=full_info_lin_plus->lines[lineno].l_rand_abcd;
			 first_saf=full_info_lin_plus->lines[lineno].l_first_saf;
			 second_saf=full_info_lin_plus->lines[lineno].l_second_saf;
			 first_safb=full_info_lin_plus->lines[lineno].l_first_safb;
			 second_safb=full_info_lin_plus->lines[lineno].l_second_safb;
			 first_safc=full_info_lin_plus->lines[lineno].l_first_safc;
			 second_safc=full_info_lin_plus->lines[lineno].l_second_safc;
			  first_safd=full_info_lin_plus->lines[lineno].l_first_safd;
			 second_safd=full_info_lin_plus->lines[lineno].l_second_safd;
			 rand_size=full_info_lin_plus->lines[lineno].l_rand_size;


			 switch(globalRNG::rand()%3)
			 {
			 case 0 :
				 rand_func(&num);
				   rand_ftr = num% nftrs;
				   first_saf  = ((double)(num%8000)/1000.0)-4.0;		
				 rand_func(&num);
				   second_saf = ((double)(num%8000)/1000.0);				
				 second_saf+=first_saf+0.01;
				 break;
			 case 1:
				 rand_func(&num);
				   rand_ftrb = num% nftrs;
				   first_safb  = ((double)(num%8000)/1000.0)-4.0;		
				 rand_func(&num);
				   second_safb = ((double)(num%8000)/1000.0);				
				 second_saf+=first_saf+0.01;
				 rand_func(&num);
				  rand_ab = num%2;
				 break;
			 case 2:
				 rand_func(&num);
				   rand_ftr = num% nftrs;
				   first_saf  = ((double)(num%8000)/1000.0)-4.0;		
				 rand_func(&num);
				   second_saf = ((double)(num%8000)/1000.0);				
				 second_saf+=first_saf+0.01;
				  rand_func(&num);
				  rand_abc = num%2;
				 break;
			 }

			 if (globalRNG::rand()%2)
			 {
				 rand_func(&num);
				double rand_size = (double)(num%25000) /5000.0;		 
				if (globalRNG::rand()%2) rand_size=-rand_size;
			 }

		}









	  	if (((it%ITERATION_FOR_ONE_CHOICE)>(ITERATION_FOR_ONE_CHOICE*0.70))   && (bestlineno>0))
		{

			 if (globalRNG::rand()%2)
			 {
				 int lineno = globalRNG::rand()%bestlineno;//full_info_lin_plus->linesno;
			 
				


				 rand_ftr=best_full_info_lin_plus.lines[lineno].l_rand_ftr;
				 rand_ftrb=best_full_info_lin_plus.lines[lineno].l_rand_ftrb;
				 rand_ftrc=best_full_info_lin_plus.lines[lineno].l_rand_ftrc;
				 rand_ftrd=best_full_info_lin_plus.lines[lineno].l_rand_ftrd;
				 rand_ab=best_full_info_lin_plus.lines[lineno].l_rand_ab;
				 rand_abc=best_full_info_lin_plus.lines[lineno].l_rand_abc;
				 rand_abcd=best_full_info_lin_plus.lines[lineno].l_rand_abcd;
				 first_saf=best_full_info_lin_plus.lines[lineno].l_first_saf;
				 second_saf=best_full_info_lin_plus.lines[lineno].l_second_saf;
				 first_safb=best_full_info_lin_plus.lines[lineno].l_first_safb;
				 second_safb=best_full_info_lin_plus.lines[lineno].l_second_safb;
				 first_safc=best_full_info_lin_plus.lines[lineno].l_first_safc;
				 second_safc=best_full_info_lin_plus.lines[lineno].l_second_safc;
				 first_safd=best_full_info_lin_plus.lines[lineno].l_first_safd;
				 second_safd=best_full_info_lin_plus.lines[lineno].l_second_safd;
				 rand_size=best_full_info_lin_plus.lines[lineno].l_rand_size;
			 } else
			 {
				  int lineno = full_info_lin_plus->linesno;
		 
				 rand_ftr=full_info_lin_plus->lines[lineno].l_rand_ftr;
				 rand_ftrb=full_info_lin_plus->lines[lineno].l_rand_ftrb;
				 rand_ftrc=full_info_lin_plus->lines[lineno].l_rand_ftrc;
				 rand_ftrd=full_info_lin_plus->lines[lineno].l_rand_ftrd;
				 rand_ab=full_info_lin_plus->lines[lineno].l_rand_ab;
				 rand_abc=full_info_lin_plus->lines[lineno].l_rand_abc;
				 rand_abcd=full_info_lin_plus->lines[lineno].l_rand_abcd;
				 first_saf=full_info_lin_plus->lines[lineno].l_first_saf;
				 second_saf=full_info_lin_plus->lines[lineno].l_second_saf;
				 first_safb=full_info_lin_plus->lines[lineno].l_first_safb;
				 second_safb=full_info_lin_plus->lines[lineno].l_second_safb;
				 first_safc=full_info_lin_plus->lines[lineno].l_first_safc;
				 second_safc=full_info_lin_plus->lines[lineno].l_second_safc;
				 first_safd=full_info_lin_plus->lines[lineno].l_first_safd;
				 second_safd=full_info_lin_plus->lines[lineno].l_second_safd;
				 rand_size=full_info_lin_plus->lines[lineno].l_rand_size;

			 }


			
			 for(int it44=0; it44<2; it44++)
			 {
			  switch(globalRNG::rand()%11)
			  {
			  case 0:						 
							first_saf *= (double)((globalRNG::rand()%400)+800.0)/1000.0;
							break;
						 
			  case 1:						 
							second_saf *= (double)((globalRNG::rand()%400)+800.0)/1000.0;
							break;

			  case 2:		first_safb *= (double)((globalRNG::rand()%400)+800.0)/1000.0;
							break;

			  case 3:		second_safb *= (double)((globalRNG::rand()%400)+800.0)/1000.0;
							break;
						 
			  case 4:       first_safc *= (double)((globalRNG::rand()%400)+800.0)/1000.0;
							break;
							 
			  case 5:		second_safc *= (double)((globalRNG::rand()%400)+800.0)/1000.0;
							break;

			  case 6:       first_safd *= (double)((globalRNG::rand()%400)+800.0)/1000.0;
							break;
							 
			  case 7:		second_safd *= (double)((globalRNG::rand()%400)+800.0)/1000.0;
							break;
				
			  case 10:
			  case 9:
			  case 8:       rand_size *= (double)((globalRNG::rand()%400)+800.0)/1000.0;
							break;
			  }
							 
			 }

			 
		}  

		 





		if (it%ITERATION_FOR_ONE_CHOICE==(ITERATION_FOR_ONE_CHOICE-1))
		{
			 int lineno = (full_info_lin_plus->linesno)-1;
		 
			 rand_ftr=full_info_lin_plus->lines[lineno].l_rand_ftr;
			 rand_ftrb=full_info_lin_plus->lines[lineno].l_rand_ftrb;
			 rand_ftrc=full_info_lin_plus->lines[lineno].l_rand_ftrc;
			 rand_ftrd=full_info_lin_plus->lines[lineno].l_rand_ftrd;
			 rand_ab=full_info_lin_plus->lines[lineno].l_rand_ab;
			 rand_abc=full_info_lin_plus->lines[lineno].l_rand_abc;
			  rand_abcd=full_info_lin_plus->lines[lineno].l_rand_abcd;
			 first_saf=full_info_lin_plus->lines[lineno].l_first_saf;
			 second_saf=full_info_lin_plus->lines[lineno].l_second_saf;
			 first_safb=full_info_lin_plus->lines[lineno].l_first_safb;
			 second_safb=full_info_lin_plus->lines[lineno].l_second_safb;
			 first_safc=full_info_lin_plus->lines[lineno].l_first_safc;
			 second_safc=full_info_lin_plus->lines[lineno].l_second_safc;
			  first_safd=full_info_lin_plus->lines[lineno].l_first_safd;
			 second_safd=full_info_lin_plus->lines[lineno].l_second_safd;
			 rand_size=(full_info_lin_plus->lines[lineno].l_rand_size/RIDGE_FACTOR)*(1-RIDGE_FACTOR);
		}


		if (lastzerofeature==full_info_lin_plus->linesno)
		{
			if (globalRNG::rand()%2) rand_size = 5; else rand_size=-5;
		}

			 

		
		if (rand_size>5) rand_size =5;
		if (rand_size<-5) rand_size = -5;
		if (first_saf>5) first_saf=5;
		if (first_safb>5) first_safb=5;
		if (first_safc>5) first_safc=5;
		if (second_saf>5) second_saf=5.01;
		if (second_safb>5) second_safb=5.01;
		if (second_safc>5) second_safc=5.01;
		if (first_saf<-5) first_saf=-5;
		if (first_safb<-5) first_safb=-5;
		if (first_safc<-5) first_safc=-5;
		if (second_saf<-5) second_saf=-4.99;
		if (second_safb<-5) second_safb=-4.99;
		if (second_safc<-5) second_safc=-4.99;

		rand_ftrd = rand_ftr ;
		rand_abcd = 1;

		if (globalRNG::rand()%2)
		{
			first_safd=9;
			second_safd=10;
		}

	 

		 



		int cnt_not_zero = 0;


		for(int i=0; i<train_size; i++)
		{
			//tosort[i] = i;

			double vala = newx1[rand_ftr*train_size+i];
			double valb = newx1[rand_ftrb*train_size+i];
			double valc = newx1[rand_ftrc*train_size+i];
			double vald = newx1[rand_ftrd*train_size+i];

		 
			double useit=0;
			double useitb=0;
			double useitc=0;
			double useitd=0;

			 
			if (vala>second_saf) useit = 1; else
			if (vala<first_saf)  useit = 0; else
			useit = (vala-first_saf)/(second_saf-first_saf);

			
			if (valb>second_safb) useitb = 1; else
			if (valb<first_safb)  useitb = 0; else
			useitb = (valb-first_safb)/(second_safb-first_safb);
			if (rand_ab) useitb = 1- useitb;


			if (valc>second_safc) useitc = 1; else
			if (valc<first_safc)  useitc = 0;   else
			useitc = (valc-first_safc)/(second_safc-first_safc);
			if (rand_abc) useitc = 1- useitc;

			if (vald>second_safd) useitd = 1; else
			if (vald<first_safd)  useitd = 0;   else
			useitd = (vald-first_safd)/(second_safd-first_safd);
			if (rand_abcd) useitd = 1- useitd;


			 useit*=useitb; 
			 useit*=useitc; 
			 useit*=useitd; 
			 useit *=rand_size;

		 

			ofc[i] = nfc_lp[i];
			nfc_lp[i] += useit;
			 
				

		}

	 	 
		

			
		int rand_act=0 ;
		int rand_dir=0; 
						 	
		 
				 
	/*	qsort(tosort , train_size , sizeof(int), lp_auccompare );

		int current0=total0;
		int sum=0;
		for(int i1=0; i1<train_size ; i1++)
		{
			if (y1[tosort[i1]]==1) sum+=current0; else current0--;
		}

			double auc = (double)sum/(double)(total0*total1);
 

		sum=0;
		for(int i1=0; i1<(train_size )*0.1; i1++)
		{
			if (y1[tosort[i1]  ]==1) sum+=1;  
		}

		double sens90 = (double)sum/(double)(total1);	*/

		 
		  	int errcntr=0;
			if (1)
			{


				
				double total0[AUCPARTS];
				double total1[AUCPARTS];
				double pauc[AUCPARTS];
				double sens90[AUCPARTS];

				memset(total0, 0, sizeof(total0));
				memset(total1, 0, sizeof(total0));
				memset(pauc, 0, sizeof(total0));
				memset(sens90, 0, sizeof(total0));

 




				int notbetter=0;

				for(int j4=0; j4<AUCPARTS; j4++)
				{
					int j3=  j4 ;
					for(int i=j3*train_size/AUCPARTS; i<(j3+1)*train_size/AUCPARTS; i++)
					{					
						tosort[i] = i;
						if (y1[i]<0.1) total0[j3]++; else total1[j3]++;

					}

					qsort(tosort+j3*train_size/AUCPARTS , train_size/AUCPARTS , sizeof(int), lp_auccompare );

					int current0=(int)total0[j3];
					int sum=0;
					for(int i=j3*train_size/AUCPARTS; i<(j3+1)*train_size/AUCPARTS; i++)
					{
						if (y1[tosort[i]]==1) sum+=current0; else current0--;
					}


					pauc[j3] = (double)sum/(double)(total0[j3]*total1[j3]);



					double part1=0;
					for(int i=j3*train_size/AUCPARTS; i<(j3+0.1)*train_size/AUCPARTS; i++)
					{
						if (y1[tosort[i]]==1) part1+=1.0;
					}


					sens90[j3] = part1/total1[j3];





					if ((j3==0) && (pauc[j3]<=bestpauc[j3])/* || (sens90[j3]<bestsens90[j3])*/)
					{
						errcntr++;
						if (errcntr>0)
						{
							notbetter = 1;
							pauc[j3] = -1e13;
							sens90[j3]=0;
							break;
						}
					} 
				}


				double cnt=0;
				double tot1=0;

				if (rand_size>0)
				{
					for(int i=0; i<train_size; i++)
					{
						if (y1[i]>0) 
						{
							tot1++;
							if (nfc_lp[i]!=ofc[i]) cnt++;
						}
					}
				} else
				{
					for(int i=0; i<train_size; i++)
					{
						if (y1[i]==0) 
						{
							tot1++;
							if (nfc_lp[i]!=ofc[i]) cnt++;
						}
					}
				}
				

				 double rauc = 0;
				 
				 
				 
				 if (rand_size > 0)
				 {
				 rauc = cnt/tot1;
					if (rauc>0.2) rauc = 0.2;
					rauc = rauc*rauc;
					rauc*=0.1; 
				 } else
				 {
					  rauc = cnt/tot1;
				  
				  
						if (rauc>0.2) rauc = 0.2;
						rauc = rauc*rauc;
						rauc*=0.01; 
				 }



				//if (rand_size>0) rauc*=1.001;
				  


			 

				//if (cnt/tot1<0.02) rauc =  -1e13;
				if (rand_size>0) if (cnt<50) rauc = -1e13;
				if (rand_size<0) if (cnt<500) rauc = -1e13;
				 
				for(int j1=0; j1<AUCPARTS; j1++)
				{
					rauc+=pauc[j1]/AUCPARTS;
				}

			 

 


				if ((rauc>0) )
				{				
					//printf("Good one.... now adding to best list\n", bestlineno);
					
					best_full_info_lin_plus.lines[bestlineno].l_rand_act = rand_act;
					best_full_info_lin_plus.lines[bestlineno].l_rand_dir = rand_dir;
					best_full_info_lin_plus.lines[bestlineno].l_rand_ftr = rand_ftr;
					best_full_info_lin_plus.lines[bestlineno].l_rand_ftrb = rand_ftrb;
					best_full_info_lin_plus.lines[bestlineno].l_rand_ftrc = rand_ftrc;
					best_full_info_lin_plus.lines[bestlineno].l_rand_ftrd = rand_ftrd;
					best_full_info_lin_plus.lines[bestlineno].l_rand_ab = rand_ab;
					best_full_info_lin_plus.lines[bestlineno].l_rand_abc = rand_abc;
					best_full_info_lin_plus.lines[bestlineno].l_rand_abcd = rand_abcd;
					best_full_info_lin_plus.lines[bestlineno].l_first_saf = first_saf;
					best_full_info_lin_plus.lines[bestlineno].l_second_saf = second_saf;
					best_full_info_lin_plus.lines[bestlineno].l_first_safb = first_safb;
					best_full_info_lin_plus.lines[bestlineno].l_second_safb = second_safb;
					best_full_info_lin_plus.lines[bestlineno].l_first_safc = first_safc;
					best_full_info_lin_plus.lines[bestlineno].l_second_safc = second_safc;
					best_full_info_lin_plus.lines[bestlineno].l_first_safd = first_safd;
					best_full_info_lin_plus.lines[bestlineno].l_second_safd = second_safd;
					best_full_info_lin_plus.lines[bestlineno].l_rand_size = rand_size;
					bestlineno++;
				}




				if ((rauc>bestrauc )  )
				{			



					for(int j1=0; j1<AUCPARTS; j1++)
					{
				 		 if (bestpauc[j1]<pauc[j1]*0.999) bestpauc[j1]=pauc[j1]*0.999;
					}
					
					 
					 
					 
					bestrauc=rauc ;
					 
					 

					 
				//	printf("%d New record: %.5f %.5f %.5f %.5f %.5f %.5f\n", it, auc, rauc, cnt1, rand_size, auc0, auc1);

					printf("%d New record: %.5f %.4f [", it, rauc, cnt/tot1);
					lastrecordit = it;
					for(int j1=0; j1<AUCPARTS; j1++)
					{
						printf("%d:%.3f:%.3f (%.3f:%.3f)",j1,pauc[j1],bestpauc[j1], sens90[j1], bestsens90[j1]);
					}
					printf("]\n");

					int lineno = full_info_lin_plus->linesno;

					full_info_lin_plus->lines[lineno].l_rand_act = rand_act;
					full_info_lin_plus->lines[lineno].l_rand_dir = rand_dir;
					full_info_lin_plus->lines[lineno].l_rand_ftr = rand_ftr;
					full_info_lin_plus->lines[lineno].l_rand_ftrb = rand_ftrb;
					full_info_lin_plus->lines[lineno].l_rand_ftrc = rand_ftrc;
					full_info_lin_plus->lines[lineno].l_rand_ftrd = rand_ftrd;
					full_info_lin_plus->lines[lineno].l_rand_ab = rand_ab;
					full_info_lin_plus->lines[lineno].l_rand_abc = rand_abc;
					full_info_lin_plus->lines[lineno].l_rand_abcd = rand_abcd;
					full_info_lin_plus->lines[lineno].l_first_saf = first_saf;
					full_info_lin_plus->lines[lineno].l_second_saf = second_saf;
					full_info_lin_plus->lines[lineno].l_first_safb = first_safb;
					full_info_lin_plus->lines[lineno].l_second_safb = second_safb;
					full_info_lin_plus->lines[lineno].l_first_safc = first_safc;
					full_info_lin_plus->lines[lineno].l_second_safc = second_safc;
					full_info_lin_plus->lines[lineno].l_first_safd = first_safd;
					full_info_lin_plus->lines[lineno].l_second_safd = second_safd;
					full_info_lin_plus->lines[lineno].l_rand_size = rand_size;

					if (lineno==LP_MAXLINES-1) break;

					//full_info_lin_plus->linesno = full_info_lin_plus->linesno+1;


					for(int i=0; i<train_size; i++)
					{
						savenfc_lp[i] = nfc_lp[i];
						 nfc_lp[i] = ofc[i];	 
							 
					}
				} else
				{
					for(int i=0; i<train_size; i++)
					{
						 nfc_lp[i] = ofc[i];	 					 
					}
				}




			}		
		 else
		{
		//	printf("%.5f\n", auc);
			for(int i=0; i<train_size; i++)
			{
				 nfc_lp[i] = ofc[i];	 
						 
			}
		}
	}




	for(int i=0; i<train_size; i++)
	{
			savenfc_lp[i] = nfc_lp[i];
	}

	double besterr=1e13;

	/*double addval = -4;
	for (addval=-4; addval<4; addval+=0.05)
	{
		double err=0;
		for(int i1=train_size; i1<(train_size ); i1++)
		{
			err += exp(-(y1[i1]-0.5)*(nfc_lp[i1]+addval));
		}

		printf("Addval: %.4f    Err: %.4f\n", addval, err);
		if (err>besterr)
		{
			addval=addval-0.05;
			for(int i=0; i<train_size;i++)
			{
				nfc_lp[i]+=addval;
				
			}
			break;
		} else
		{
			besterr =err;

		}

	}*/


	//for(int i=0; i<train_size; i++)
	//	nfc_lp[i]=addvalue;

	time(&end);
	fprintf(stderr, "Elpased time: %d\n", (long int)(end)-(long int)(start));

	return;

}

int write_full_linear_plus(full_linear_plus_info_t *linear_info, char *fname) {

	FILE *fp ;
	if ((fp = safe_fopen(fname,"w",false))==NULL) {
		fprintf(stderr,"Cannot open %s for reading\n",fname) ;
		return -1 ;
	}

	write_full_linear_plus(linear_info,fp) ;
	fclose(fp) ;
	return 0;
}

void write_full_linear_plus(full_linear_plus_info_t *linear_info, FILE *fp) {

	fprintf(fp,"NFTRS = %d\n",linear_info->nftrs);
	for(int j=0; j<linear_info->nftrs; j++)
	{
		fprintf(fp,"AVR %d = %f\n",j, linear_info->avrgs[j]);
		fprintf(fp,"STD %d = %f\n",j, linear_info->stdevs[j]) ;
	}
 

	fprintf(fp,"Lines = %d\n",linear_info->linesno);
	for (int i=0; i<linear_info->linesno; i++)
	{
		fprintf(fp,"l_rand_ftr %d = %d\n",i,linear_info->lines[i].l_rand_ftr) ; 
		fprintf(fp,"l_rand_ftrb %d = %d\n",i,linear_info->lines[i].l_rand_ftrb) ;
		fprintf(fp,"l_rand_ftrc %d = %d\n",i,linear_info->lines[i].l_rand_ftrc) ;
		fprintf(fp,"l_rand_ftrd %d = %d\n",i,linear_info->lines[i].l_rand_ftrd) ;
		fprintf(fp,"l_rand_act %d = %d\n",i,linear_info->lines[i].l_rand_act) ; 
		fprintf(fp,"l_rand_dir %d = %d\n",i,linear_info->lines[i].l_rand_dir) ; 
		fprintf(fp,"l_rand_ab %d = %d\n",i,linear_info->lines[i].l_rand_ab) ; 
		fprintf(fp,"l_rand_abc %d = %d\n",i,linear_info->lines[i].l_rand_abc) ;
		fprintf(fp,"l_rand_abcd %d = %d\n",i,linear_info->lines[i].l_rand_abcd) ;
		fprintf(fp,"l_first_saf %d = %f\n",i,linear_info->lines[i].l_first_saf) ; 
		fprintf(fp,"l_second_saf %d = %f\n",i,linear_info->lines[i].l_second_saf) ; 
		fprintf(fp,"l_first_safb %d = %f\n",i,linear_info->lines[i].l_first_safb) ; 
		fprintf(fp,"l_second_safb %d = %f\n",i,linear_info->lines[i].l_second_safb) ; 
		fprintf(fp,"l_first_safc %d = %f\n",i,linear_info->lines[i].l_first_safc) ; 
		fprintf(fp,"l_second_safc %d = %f\n",i,linear_info->lines[i].l_second_safc) ; 
		fprintf(fp,"l_first_safd %d = %f\n",i,linear_info->lines[i].l_first_safd) ; 
		fprintf(fp,"l_second_safd %d = %f\n",i,linear_info->lines[i].l_second_safd) ;
		fprintf(fp,"l_rand_size %d = %f\n",i,linear_info->lines[i].l_rand_size) ; 
	}

	fprintf(fp,"min_auc_ftr = %d\n",min_auc_ftr) ; 
	fprintf(stderr, "min_auc_ftr = %d\n",min_auc_ftr) ; 

	return ;
}

int read_full_linear_plus(full_linear_plus_info_t *linear_info, char *fname) {

	fprintf(stderr, "Start reading ....");
	FILE *fp ;
	if ((fp = fopen(fname,"r"))==NULL) {
		fscanf(stderr,"Cannot open %s for reading\n",fname) ;
		return -1 ;
	}

	fscanf(fp,"NFTRS = %d\n",&linear_info->nftrs);

	

	linear_info->avrgs = (double *) calloc(linear_info->nftrs, sizeof(double));
	linear_info->stdevs = (double *) calloc(linear_info->nftrs, sizeof(double));



	for(int j=0; j<linear_info->nftrs; j++)
	{
		 

		int tmp;
		fscanf(fp,"AVR %d = %lf\n",&tmp, &(linear_info->avrgs[j]));
		fscanf(fp,"STD %d = %lf\n",&tmp, &(linear_info->stdevs[j])) ;

		 
	}
	
/*
struct full_linear_plus_info_line_t
{
	int l_rand_ftr;
	int l_rand_act;
	int l_rand_dir;
	int rand_ab;
	double rand_saf;
	double rand_size;
};
*/
 


	fscanf(fp,"Lines = %d\n",&linear_info->linesno);

	linear_info->lines = (full_linear_plus_info_line_t *) calloc(linear_info->linesno, sizeof(full_linear_plus_info_line_t));


	for (int i=0; i<linear_info->linesno; i++)
	{
		int tmp;
		fscanf(fp,"l_rand_ftr %d = %d\n",&tmp,&linear_info->lines[i].l_rand_ftr) ; 
		fscanf(fp,"l_rand_ftrb %d = %d\n",&tmp,&linear_info->lines[i].l_rand_ftrb) ; 
		fscanf(fp,"l_rand_ftrc %d = %d\n",&tmp,&linear_info->lines[i].l_rand_ftrc) ; 
		fscanf(fp,"l_rand_ftrd %d = %d\n",&tmp,&linear_info->lines[i].l_rand_ftrd) ; 
		fscanf(fp,"l_rand_act %d = %d\n",&tmp,&linear_info->lines[i].l_rand_act) ; 
		fscanf(fp,"l_rand_dir %d = %d\n",&tmp,&linear_info->lines[i].l_rand_dir) ; 
		fscanf(fp,"l_rand_ab %d = %d\n",&tmp,&linear_info->lines[i].l_rand_ab) ; 
		fscanf(fp,"l_rand_abc %d = %d\n",&tmp,&linear_info->lines[i].l_rand_abc) ; 
		fscanf(fp,"l_rand_abcd %d = %d\n",&tmp,&linear_info->lines[i].l_rand_abcd) ;
		fscanf(fp,"l_first_saf %d = %lf\n",&tmp,&linear_info->lines[i].l_first_saf) ; 
		fscanf(fp,"l_second_saf %d = %lf\n",&tmp,&linear_info->lines[i].l_second_saf) ; 
		fscanf(fp,"l_first_safb %d = %lf\n",&tmp,&linear_info->lines[i].l_first_safb) ; 
		fscanf(fp,"l_second_safb %d = %lf\n",&tmp,&linear_info->lines[i].l_second_safb) ; 
		fscanf(fp,"l_first_safc %d = %lf\n",&tmp,&linear_info->lines[i].l_first_safc) ; 
		fscanf(fp,"l_second_safc %d = %lf\n",&tmp,&linear_info->lines[i].l_second_safc) ; 
			fscanf(fp,"l_first_safd %d = %lf\n",&tmp,&linear_info->lines[i].l_first_safd) ; 
		fscanf(fp,"l_second_safd %d = %lf\n",&tmp,&linear_info->lines[i].l_second_safd) ; 
		fscanf(fp,"l_rand_size %d = %lf\n",&tmp,&linear_info->lines[i].l_rand_size) ; 
		fprintf(stderr,"%d %.4f\n", i, linear_info->lines[i].l_rand_size);
	}

	fscanf(fp,"min_auc_ftr = %d\n",&min_auc_ftr) ; 
	fprintf(stderr, "min_auc_ftr = %d\n",min_auc_ftr) ; 
	//fprintf(stderr, "Finished reading .... %d...");

	return 0 ;
}
 
int get_linear_plus_f(double *x, double *y, double *w, int nrows, int ncols, char *fname) {
	// Learn
	full_linear_plus_info_t full_info_lin_plus;

	 get_linear_plus(x,y,nrows,ncols,&full_info_lin_plus) ;
	 write_full_linear_plus(&full_info_lin_plus, fname);

	return 0 ;
}

int linear_plus_predict(double *x, double *preds, int nrows, int ncols, char *fname) {

	fprintf(stderr, "In Linear plus predict\n");
	// Read
	 full_linear_plus_info_t linear_plus_struct ;
	if (read_full_linear_plus(&linear_plus_struct,fname) == -1) {
		fprintf(stderr,"Reading failed\n") ;
		return -1 ;
	}

	// Predict
	if (predict_linear_plus(x,nrows,ncols,&linear_plus_struct,preds)==-1) {
		fprintf(stderr,"Linear plus prediction failed\n") ;
		return -1 ;
	}

	// Clearing ?

	return 0 ;
}

int linear_plus_predict(double *x, double *preds, int nrows, int ncols, unsigned char *model) {

	int size ;

	// DeSerializae 
	full_linear_plus_info_t linear_plus_struct ;
	if ((size = linear_plus_deserialize(model,&linear_plus_struct))==-1) {
		fprintf(stderr,"DeSerialization failed\n") ;
		return -1 ;
	}

	// Predict
	if (predict_linear_plus(x,nrows,ncols,&linear_plus_struct,preds)==-1) {
		fprintf(stderr,"Linear plus prediction failed\n") ;
		return -1 ;
	}

	// Clearing ?

	return size ;
}

int predict_linear_plus(double *x, int nrows, int ncols, full_linear_plus_info_t *linear_info, double *preds)
{
	double *newx1 = (double *) calloc(sizeof(double), nrows*ncols);

	double *usevec  = (double *) calloc(sizeof(double), nrows);


	for(int i=0; i<nrows; i++)
		preds[i]=0;


	for(int j=0; j<ncols; j++)
	for(int i=0; i<nrows; i++)
	{
	
		newx1[j*nrows+i]=x[j*nrows+i];
		newx1[j*nrows+i]-=linear_info->avrgs[j];
		newx1[j*nrows+i]/=linear_info->stdevs[j];
	}


	for(int j1=0; j1< linear_info->linesno; j1++)
	{
		
		int rand_act=linear_info->lines[j1].l_rand_act;  
		int rand_dir=linear_info->lines[j1].l_rand_dir;
		int rand_ftr=linear_info->lines[j1].l_rand_ftr;
		int rand_ftrb=linear_info->lines[j1].l_rand_ftrb;
		int rand_ftrc=linear_info->lines[j1].l_rand_ftrc;
		int rand_ftrd=linear_info->lines[j1].l_rand_ftrd;
		int rand_ab=linear_info->lines[j1].l_rand_ab;
		int rand_abc=linear_info->lines[j1].l_rand_abc;
		int rand_abcd=linear_info->lines[j1].l_rand_abcd;
		double first_saf=linear_info->lines[j1].l_first_saf;
		double second_saf=linear_info->lines[j1].l_second_saf;
		double first_safb=linear_info->lines[j1].l_first_safb;
		double second_safb=linear_info->lines[j1].l_second_safb;
		double first_safc=linear_info->lines[j1].l_first_safc;
		double second_safc=linear_info->lines[j1].l_second_safc;
		double rand_size=linear_info->lines[j1].l_rand_size;
		double first_safd=linear_info->lines[j1].l_first_safd;
		double second_safd=linear_info->lines[j1].l_second_safd;
	

		fprintf(stderr,"%.6f\n", rand_size);


		 
		for(int i=0; i<nrows; i++)
		{
			double vala = newx1[rand_ftr*nrows+i];
			double valb = newx1[rand_ftrb*nrows+i];
			double valc = newx1[rand_ftrc*nrows+i];
			double vald = newx1[rand_ftrd*nrows+i];

		 

			double useit=0;
			double useitb=0;
			double useitc=0;
			double useitd=0;

			 
			if (vala>second_saf) useit = 1; else
			if (vala<first_saf)  useit = 0; else
			useit = (vala-first_saf)/(second_saf-first_saf);

			
			if (valb>second_safb) useitb = 1; else
			if (valb<first_safb)  useitb = 0; else
			useitb = (valb-first_safb)/(second_safb-first_safb);
			if (rand_ab) useitb = 1- useitb;


			if (valc>second_safc) useitc = 1; else
			if (valc<first_safc)  useitc = 0;   else
			useitc = (valc-first_safc)/(second_safc-first_safc);
			if (rand_abc) useitc = 1- useitc;

			if (vald>second_safd) useitd = 1; else
			if (vald<first_safd)  useitd = 0;   else
			useitd = (vald-first_safd)/(second_safd-first_safd);
			if (rand_abcd) useitd = 1- useitd;


			 useit*=useitb; 
			 useit*=useitc; 
			 useit*=useitd; 
			 useit *=rand_size;


				

			preds[i]+=useit;

		}

	}

	return 0;
};

// (De)Serialization

size_t get_linear_plus_size(full_linear_plus_info_t *linear_info) {

	size_t size = 2 * sizeof(int) ; // NFTRS + LinesNo
	size += linear_info->nftrs * 2 * sizeof(double) ; // Avrgs + Stdevs
	size += linear_info->linesno * sizeof(full_linear_plus_info_line_t) ; // Lines 

	return size ;
}

int linear_plus_serialize(full_linear_plus_info_t *linear_info, unsigned char *linear_data) {

	size_t idx = 0 ;
	memcpy(linear_data,&(linear_info->nftrs),sizeof(int)) ; idx += sizeof(int) ;
	memcpy(linear_data+idx,&(linear_info->linesno),sizeof(int)) ; idx += sizeof(int) ;
	memcpy(linear_data+idx,linear_info->avrgs,linear_info->nftrs*sizeof(double)) ; idx += linear_info->nftrs*sizeof(double) ;
	memcpy(linear_data+idx,linear_info->stdevs,linear_info->nftrs*sizeof(double)) ; idx += linear_info->nftrs*sizeof(double) ;
	memcpy(linear_data+idx,linear_info->lines,linear_info->linesno*sizeof(full_linear_plus_info_line_t)) ; idx += linear_info->linesno*sizeof(full_linear_plus_info_line_t) ;

	return (int) idx ;
}

int linear_plus_deserialize(unsigned char *linear_data, full_linear_plus_info_t *linear_info) {

	size_t idx = 0 ;
	memcpy(&(linear_info->nftrs),linear_data,sizeof(int)) ; idx += sizeof(int) ;
	memcpy(&(linear_info->linesno),linear_data+idx,sizeof(int)) ; idx += sizeof(int) ;

	linear_info->avrgs = (double *) malloc(linear_info->nftrs*sizeof(double)) ;
	linear_info->stdevs = (double *) malloc(linear_info->nftrs*sizeof(double)) ;
	linear_info->lines = (full_linear_plus_info_line_t *) malloc(linear_info->linesno*sizeof(full_linear_plus_info_line_t)) ;
	if (linear_info->avrgs==NULL || linear_info->stdevs==NULL || linear_info->lines==NULL) {
		fprintf(stderr,"Allocation Failed\n") ;
		return -1 ;
	}

	memcpy(linear_info->avrgs,linear_data+idx,linear_info->nftrs*sizeof(double)) ; idx += linear_info->nftrs*sizeof(double) ;
	memcpy(linear_info->stdevs,linear_data+idx,linear_info->nftrs*sizeof(double)) ; idx += linear_info->nftrs*sizeof(double) ;
	memcpy(linear_info->lines,linear_data+idx,linear_info->linesno*sizeof(full_linear_plus_info_line_t)) ; idx += linear_info->linesno*sizeof(full_linear_plus_info_line_t) ;

	return (int) idx ;
}

void clear_linear_plus(full_linear_plus_info_t *linear_info) {
	free(linear_info->avrgs) ;
	free(linear_info->stdevs) ;
	free(linear_info->lines) ;
}

