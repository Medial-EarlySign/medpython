// ROCR utilities : Analyze regression/classification results
# include "ROCR.h"

int generate_graph(float *labels, float *predictions, int nsamples, char *xmode, char *ymode, int resolution, float **x, float **y) {

	// Check modes
	int x_mode,y_mode ;

	if (strcmp(xmode,"rec")==0 || strcmp(xmode,"sens")==0) {
		x_mode = 1 ;
	} else if (strcmp(xmode,"fpr")==0) {
		x_mode = 2 ;
	} else if (strcmp(xmode,"rpp")==0) {
		x_mode = 3 ;
	} else if (strcmp(xmode,"scr")==0) {
		x_mode = 4 ;
	} else {
		printf("error: Unknown x mode : %s in ROCR\n",xmode) ;
		return -1 ;
	}

	if (strcmp(ymode,"tpr")==0) {
		y_mode = 1 ;
	} else if (strcmp(ymode,"ppv")==0 || strcmp(ymode,"prec")==0) {
		y_mode = 2 ;
	} else if (strcmp(ymode,"lift")==0) {
		y_mode = 3 ;
	} else {
		printf("error : Unknown y mode : %s in ROCR\n",ymode) ;
		return -1 ;
	}


	// Order Predictions
	float *values ;
	int   *inds ;

	if ((values = (float *) malloc (sizeof (float)*nsamples)) == NULL) {
		printf ("error : cannot allocate %d prediction values",nsamples) ;
		return -1 ;
	}

	if ((inds = (int *) malloc (sizeof (float)*nsamples)) == NULL) {
		printf ("error : cannot allocate %d indices",nsamples) ;
		return -1 ;
	}

	for (int i=0; i<nsamples; i++) {
		values[i] = predictions[i] ;
		inds[i] = i ;
	}

	for(int i=0; i<nsamples; i++) {
		for (int j=0; j<i; j++)
		{
			if (values[i]<values[j])
			{
				float temp_value = values[i];
				values[i] = values[j];
				values[j] = temp_value;

				int temp_ind = inds[i];
				inds[i] = inds[j];
				inds[j] = temp_ind;
			}
		}
	}

	// Alocate
	float *init_x ;
	float *init_y ;

	if ((init_x = (float *) malloc (sizeof (float)*nsamples)) == NULL) {
		printf ("error : cannot allocate x for %d samples",nsamples) ;
		return -1 ;
	}

	if ((init_y = (float *) malloc (sizeof (float)*nsamples)) == NULL) {
		printf ("error : cannot allocate y for %d samples",nsamples) ;
		return -1 ;
	}

	// P and N
	int pos = 0;
	int neg = 0 ;

	for (int i=0; i<nsamples; i++) {
		if (labels[i] > 0.0) {
			pos++ ;
		} else {
			neg++ ;
		}
	}	
	
	// Calculate FP and TP
	float *fp ;
	float *tp ;

	if ((fp = (float *) malloc (sizeof (float)*nsamples)) == NULL) {
		printf ("error : cannot allocate fp for %d samples",nsamples) ;
		return -1 ;
	}

	if ((tp = (float *) malloc (sizeof (float)*nsamples)) == NULL) {
		printf ("error : cannot allocate tp for %d samples",nsamples) ;
		return -1 ;
	}


	for (int i=nsamples-1; i>=0; i--) {		
		if (labels[inds[i]] > 0.0) {
			if (i==nsamples-1) {
				tp[i] = 1 ;
				fp[i] = 0 ;
			} else {
				tp[i] = tp[i+1] + 1 ;
				fp[i] = fp[i+1] ;
			}
		} else {
			if (i==nsamples-1) {
				tp[i] = 0 ;
				fp[i] = 1 ;
			} else {
				tp[i] = tp[i+1] ;
				fp[i] = fp[i+1] + 1 ;
			}
		}
	}

	// Calculate X and Y
	float min_x=0 ;
	float min_y=0 ;
	float max_x=0 ;
	float max_y=0 ;

	for (int i=0; i<nsamples; i++) {
		if (x_mode == 1) { // TPR,Recall,Sensitivity
			init_x[i] = tp[i]/pos ;
		} else if (x_mode == 2) { // FPR
			init_x[i] = fp[i]/neg ;
		} else if (x_mode == 3) { // RPP
			init_x[i] = tp[i]/pos ;
		} else if (x_mode == 4) { // SCR
			init_x[i] = predictions[inds[i]] ;
		}
			

		if (i==0 || init_x[i]>max_x)
			max_x = init_x[i] ;
		if (i==0 || init_x[i]<min_x)
			min_x = init_x[i] ;

		if (y_mode ==2) { // PPV,Precision
			init_y[i] = tp[i]/(tp[i]+fp[i]) ;
		} else if (y_mode == 1) { // TPR
			init_y[i] = tp[i]/pos ;
		} else if (y_mode == 3) { // Lift
			init_y[i] = (tp[i]/pos)/((tp[i]+fp[i])/(pos+neg)) ;
		}

		if (i==0 || init_y[i]>max_y)
			max_y = init_y[i] ;
		if (i==0 || init_y[i]<min_y)
			max_y = init_y[i] ;
	}

	// Take care of required resolution
	if (resolution == 0) {
		if (((*x) = (float *) malloc (sizeof (float)*nsamples)) == NULL) {
			printf ("error : cannot allocate x for %d samples",nsamples) ;
			return -1 ;
		}
		memcpy(*x,init_x,sizeof(float)*nsamples) ;

		if (((*y) = (float *) malloc (sizeof (float)*nsamples)) == NULL) {
			printf ("error : cannot allocate y for %d samples",nsamples) ;
			return -1 ;
		}
		memcpy(*y,init_y,sizeof(float)*nsamples) ;
	} else {
		if (((*x) = (float *) malloc (sizeof (float)*resolution)) == NULL) {
			printf ("error : cannot allocate x for %d bins",resolution) ;
			return -1 ;
		}

		if (((*y) = (float *) malloc (sizeof (float)*resolution)) == NULL) {
			printf ("error : cannot allocate y for %d bins",resolution) ;
			return -1 ;
		}
		memset(*y,0,sizeof(float)*resolution) ;

		int *num ;
		if ((num = (int *) malloc (sizeof (int)*resolution)) == NULL) {
			printf ("error : cannot allocate counter for %d bins",resolution) ;
			return -1 ;
		}
		memset(num,0,sizeof(int)*resolution) ;

		float min_x = 0 ;
		float max_x = 0 ;
		for (int i=0; i<nsamples; i++) {
			if (i==0 || min_x > init_x[i])
				min_x = init_x[i] ;
			if (i==0 || max_x < init_x[i])
				max_x = init_x[i] ;
		}

		float x_bin_size = (max_x-min_x)/resolution ;

		for (int i=0; i<nsamples; i++) {
			int x_bin ;
			if (init_x[i] == max_x) {
				x_bin = resolution - 1 ;
			} else {
				x_bin = (int)((init_x[i]-min_x)/x_bin_size) ;
			}

			num[x_bin]++ ;
			(*y)[x_bin]+=init_y[i] ;
		}

		for (int i=0; i<resolution; i++) {
			(*x)[i] = min_x + i*x_bin_size ;
			if ((*y)[i] > 0) {
				(*y)[i] /= num[i] ;
			} else if (i==0) {
				(*y)[i] = 0 ;
			} else {
				(*y)[i] = (*y)[i-1] ;
			}
		}
		free(num) ;
	}

	free(values) ;
	free(inds) ;
	free(init_x) ;
	free(init_y) ;
	free(tp) ;

	return 0 ;
}
	
int get_optimal_confusion(float *labels, float *predictions, int nsamples, int confusion[2][2]) {

	// Order Predictions
	float *values ;
	int   *inds ;

	if ((values = (float *) malloc (sizeof (float)*nsamples)) == NULL) {
		printf ("error : cannot allocate %d prediction values",nsamples) ;
		return -1 ;
	}

	if ((inds = (int *) malloc (sizeof (float)*nsamples)) == NULL) {
		printf ("error : cannot allocate %d indices",nsamples) ;
		return -1 ;
	}

	for (int i=0; i<nsamples; i++) {
		values[i] = predictions[i] ;
		inds[i] = i ;
	}

	for(int i=0; i<nsamples; i++) {
		for (int j=0; j<i; j++)
		{
			if (values[i]<values[j])
			{
				float temp_value = values[i];
				values[i] = values[j];
				values[j] = temp_value;

				int temp_ind = inds[i];
				inds[i] = inds[j];
				inds[j] = temp_ind;
			}
		}
	}

	// Get starting point
	int optimal_i = 0 ;
	int optimal_err = 0 ;
	for (int i=0; i<nsamples; i++) {
		if (labels[i] < 0)
			optimal_err++ ;
	}
	confusion[0][0] = confusion[0][1] = 0 ;
	confusion[1][0] = optimal_err ;
	confusion[1][1] = nsamples - optimal_err ;

	int current_confusion[2][2] ;
	memcpy(current_confusion,confusion,2*2*sizeof(int)) ;

	// Search better places
	int current_i = 0 ;
	float current_value = values[inds[current_i]] ;
	int current_err = optimal_err ;

	while (current_i < nsamples) {
//		fprintf(stderr,"%d %f %f\n",inds[current_i],values[current_i],labels[inds[current_i]]) ;
		current_i++ ;
		if (current_i == nsamples || values[current_i] != current_value) {
			if (current_err < optimal_err) {
				optimal_err = current_err ;
				memcpy(confusion,current_confusion,2*2*sizeof(int)) ;
			}
		}

		if (current_i != nsamples) {
			if (labels[inds[current_i]] < 0) {
				current_err -- ;
				current_confusion[1][0]-- ;
				current_confusion[0][0]++ ;
			} else {
				current_err ++ ;
				current_confusion[1][1]-- ;
				current_confusion[0][1]++ ;
			}
		}
	}

	return 0 ;
}

int generate_hist(float *labels, float *predictions, int nsamples, char *mode, int resolution, float **x, float **y, int **count) {

	// Check modes
	int mode_id ;
	if (strcmp(mode,"prec")==0) {
		mode_id = 1 ;
	} else {
		printf("error: Unknown mode : %s in ROCR\n",mode) ;
		return -1 ;
	}

	// Order Predictions
	float *values ;
	int   *inds ;

	if ((values = (float *) malloc (sizeof (float)*nsamples)) == NULL) {
		printf ("error : cannot allocate %d prediction values",nsamples) ;
		return -1 ;
	}

	if ((inds = (int *) malloc (sizeof (float)*nsamples)) == NULL) {
		printf ("error : cannot allocate %d indices",nsamples) ;
		return -1 ;
	}

	for (int i=0; i<nsamples; i++) {
		values[i] = predictions[i] ;
		inds[i] = i ;
	}

	for(int i=0; i<nsamples; i++) {
		for (int j=0; j<i; j++)
		{
			if (values[i]<values[j])
			{
				float temp_value = values[i];
				values[i] = values[j];
				values[j] = temp_value;

				int temp_ind = inds[i];
				inds[i] = inds[j];
				inds[j] = temp_ind;
			}
		}
	}

	// Allocate 
	if (((*x) = (float *) malloc (sizeof (float)*(resolution + 1))) == NULL) {
		printf ("error : cannot allocate x for %d bins",resolution) ;
		return -1 ;
	}

	if (((*count) = (int *) malloc (sizeof (int)*resolution)) == NULL) {
		printf ("error : cannot allocate counter for %d bins",resolution) ;
		return -1 ;
	}

	if (((*y) = (float *) malloc (sizeof (float)*resolution)) == NULL) {
		printf ("error : cannot allocate y for %d bins",resolution) ;
		return -1 ;
	}

	// Loop and Collect
	int bin_size = nsamples/resolution ;
	int nlabled = 0 ;
	int n = 0 ;
	int ibin ;

	for (int i=0; i<nsamples ;i++) {
		if (i == 0) {
			ibin = 0 ;
			(*x)[ibin] = predictions[inds[i]] ;
		} else if (i%bin_size == 0 && i+bin_size <= nsamples) {
			(*y)[ibin] = (float)(nlabled+0.0)/n ;
			(*count)[ibin] = n ;
			nlabled = n = 0 ;

			ibin++ ;
			(*x)[ibin] = predictions[inds[i]] ;
		}

		n++ ;
		if (labels[inds[i]] > 0.0)
			nlabled++ ;
	}

	(*y)[ibin] = (float)(nlabled+0.0)/n ;
	(*count)[ibin] = n ;
	(*x)[ibin+1] = predictions[inds[nsamples-1]] ;

	free(values) ;
	free(inds) ;

	return 0 ;
}

int get_auc(float *labels, float *predictions, int nsamples, int resolution, float *auc) {

	// Order Predictions
	float *values ;
	int   *inds ;

	if ((values = (float *) malloc (sizeof (float)*nsamples)) == NULL) {
		fprintf (stderr,"error : cannot allocate %d prediction values\n",nsamples) ;
		return -1 ;
	}

	if ((inds = (int *) malloc (sizeof (float)*nsamples)) == NULL) {
		fprintf (stderr,"error : cannot allocate %d indices\n",nsamples) ;
		return -1 ;
	}

	for (int i=0; i<nsamples; i++) {
		values[i] = predictions[i] ;
		inds[i] = i ;
	}

	for(int i=0; i<nsamples; i++) {
		for (int j=0; j<i; j++)
		{
			if (values[i]<values[j])
			{
				float temp_value = values[i];
				values[i] = values[j];
				values[j] = temp_value;

				int temp_ind = inds[i];
				inds[i] = inds[j];
				inds[j] = temp_ind;
			}
		}
	}

	// Alocate
	float *init_x ;
	float *init_y ;

	if ((init_x = (float *) malloc (sizeof (float)*nsamples)) == NULL) {
		fprintf (stderr,"error : cannot allocate x for %d samples\n",nsamples) ;
		return -1 ;
	}

	if ((init_y = (float *) malloc (sizeof (float)*nsamples)) == NULL) {
		fprintf (stderr,"error : cannot allocate y for %d samples\n",nsamples) ;
		return -1 ;
	}

	// P and N
	int pos = 0;
	int neg = 0 ;

	for (int i=0; i<nsamples; i++) {
		if (labels[i] > 0.0) {
			pos++ ;
		} else {
			neg++ ;
		}
	}	
	
	// Calculate FP and TP
	float *fp ;
	float *tp ;

	if ((fp = (float *) malloc (sizeof (float)*nsamples)) == NULL) {
		fprintf (stderr,"error : cannot allocate fp for %d samples\n",nsamples) ;
		return -1 ;
	}

	if ((tp = (float *) malloc (sizeof (float)*nsamples)) == NULL) {
		fprintf (stderr,"error : cannot allocate tp for %d samples\n",nsamples) ;
		return -1 ;
	}


	for (int i=nsamples-1; i>=0; i--) {		
		if (labels[inds[i]] > 0.0) {
			if (i==nsamples-1) {
				tp[i] = 1 ;
				fp[i] = 0 ;
			} else {
				tp[i] = tp[i+1] + 1 ;
				fp[i] = fp[i+1] ;
			}
		} else {
			if (i==nsamples-1) {
				tp[i] = 0 ;
				fp[i] = 1 ;
			} else {
				tp[i] = tp[i+1] ;
				fp[i] = fp[i+1] + 1 ;
			}
		}
	}

	// Calculate X and Y

	for (int i=0; i<nsamples; i++) {
		init_x[i] = fp[i]/neg ;
		init_y[i] = tp[i]/pos ;
	}

	// Take care of required resolution
	float *x ;
	float *y ;

	if ((x = (float *) malloc (sizeof (float)*resolution)) == NULL) {
		fprintf (stderr,"error : cannot allocate x for %d bins\n",resolution) ;
		return -1 ;
	}

	if ((y = (float *) malloc (sizeof (float)*resolution)) == NULL) {
		fprintf (stderr,"error : cannot allocate y for %d bins\n",resolution) ;
		return -1 ;
	}
	memset(y,0,sizeof(float)*resolution) ;

	int *num ;
	if ((num = (int *) malloc (sizeof (int)*resolution)) == NULL) {
		fprintf (stderr,"error : cannot allocate counter for %d bins\n",resolution) ;
		return -1 ;
	}
	memset(num,0,sizeof(int)*resolution) ;

	float min_x = 0 ;
	float max_x = 0 ;
	for (int i=0; i<nsamples; i++) {
		if (i==0 || min_x > init_x[i])
			min_x = init_x[i] ;
		if (i==0 || max_x < init_x[i])
			max_x = init_x[i] ;
	}

	float x_bin_size = (max_x-min_x)/resolution ;
	for (int i=0; i<nsamples; i++) {
		int x_bin ;
		if (init_x[i] == max_x) {
			x_bin = resolution - 1 ;
		} else {
			x_bin = (int)((init_x[i]-min_x)/x_bin_size) ;
		}

		num[x_bin]++ ;
		y[x_bin]+=init_y[i] ;
	}
		
	for (int i=0; i<resolution; i++) {
		x[i] = min_x + i*x_bin_size ;
		if (y[i] > 0) {
			y[i] /= num[i] ;
		} else if (i==0) {
			y[i] = 0 ;
		} else {
			y[i] = y[i-1] ;
		}
	}

	// Calculate AUC

	(*auc) = 0 ;
	for (int i=0; i<resolution-1; i++)
		(*auc) += x_bin_size * (y[i] + y[i+1])/2 ;


	free(values) ;
	free(inds) ;
	free(init_x) ;
	free(init_y) ;
	free(tp) ;
	free(fp) ;
	free(x) ;
	free(y) ;
	free(num) ;

	return 0 ;
}

int get_auc(float *labels, float *predictions, int tot_nsamples, int resolution, int *marked, float *auc) {

	// Order Predictions
	float *values ;
	int   *inds ;

	if ((values = (float *) malloc (sizeof (float)*tot_nsamples)) == NULL) {
		fprintf (stderr,"error : cannot allocate %d prediction values\n",tot_nsamples) ;
		return -1 ;
	}

	if ((inds = (int *) malloc (sizeof (float)*tot_nsamples)) == NULL) {
		fprintf (stderr,"error : cannot allocate %d indices\n",tot_nsamples) ;
		return -1 ;
	}

	int nsamples = 0 ;
	for (int i=0; i<tot_nsamples; i++) {
		if (marked[i]) {
			values[nsamples] = predictions[i] ;
			inds[nsamples++] = i ;
		}
	}

	for(int i=0; i<nsamples; i++) {
		for (int j=0; j<i; j++)
		{
			if (values[i]<values[j])
			{
				float temp_value = values[i];
				values[i] = values[j];
				values[j] = temp_value;

				int temp_ind = inds[i];
				inds[i] = inds[j];
				inds[j] = temp_ind;
			}
		}
	}

	// Alocate
	float *init_x ;
	float *init_y ;

	if ((init_x = (float *) malloc (sizeof (float)*nsamples)) == NULL) {
		fprintf (stderr,"error : cannot allocate x for %d samples\n",nsamples) ;
		return -1 ;
	}

	if ((init_y = (float *) malloc (sizeof (float)*nsamples)) == NULL) {
		fprintf (stderr,"error : cannot allocate y for %d samples\n",nsamples) ;
		return -1 ;
	}

	// P and N
	int pos = 0;
	int neg = 0 ;

	for (int i=0; i<nsamples; i++) {
		if (labels[inds[i]] > 0.0) {
			pos++ ;
		} else {
			neg++ ;
		}
	}	
	
	// Calculate FP and TP
	float *fp ;
	float *tp ;

	if ((fp = (float *) malloc (sizeof (float)*nsamples)) == NULL) {
		fprintf (stderr,"error : cannot allocate fp for %d samples\n",nsamples) ;
		return -1 ;
	}

	if ((tp = (float *) malloc (sizeof (float)*nsamples)) == NULL) {
		fprintf (stderr,"error : cannot allocate tp for %d samples\n",nsamples) ;
		return -1 ;
	}


	for (int i=nsamples-1; i>=0; i--) {		
		if (labels[inds[i]] > 0.0) {
			if (i==nsamples-1) {
				tp[i] = 1 ;
				fp[i] = 0 ;
			} else {
				tp[i] = tp[i+1] + 1 ;
				fp[i] = fp[i+1] ;
			}
		} else {
			if (i==nsamples-1) {
				tp[i] = 0 ;
				fp[i] = 1 ;
			} else {
				tp[i] = tp[i+1] ;
				fp[i] = fp[i+1] + 1 ;
			}
		}
	}

	// Calculate X and Y

	for (int i=0; i<nsamples; i++) {
		init_x[i] = fp[i]/neg ;
		init_y[i] = tp[i]/pos ;
	}

	// Take care of required resolution
	float *x ;
	float *y ;

	if ((x = (float *) malloc (sizeof (float)*resolution)) == NULL) {
		fprintf (stderr,"error : cannot allocate x for %d bins\n",resolution) ;
		return -1 ;
	}

	if ((y = (float *) malloc (sizeof (float)*resolution)) == NULL) {
		fprintf (stderr,"error : cannot allocate y for %d bins\n",resolution) ;
		return -1 ;
	}
	memset(y,0,sizeof(float)*resolution) ;

	int *num ;
	if ((num = (int *) malloc (sizeof (int)*resolution)) == NULL) {
		fprintf (stderr,"error : cannot allocate counter for %d bins\n",resolution) ;
		return -1 ;
	}
	memset(num,0,sizeof(int)*resolution) ;

	float min_x = 0 ;
	float max_x = 0 ;
	for (int i=0; i<nsamples; i++) {
		if (i==0 || min_x > init_x[i])
			min_x = init_x[i] ;
		if (i==0 || max_x < init_x[i])
			max_x = init_x[i] ;
	}

	float x_bin_size = (max_x-min_x)/resolution ;
	for (int i=0; i<nsamples; i++) {
		int x_bin ;
		if (init_x[i] == max_x) {
			x_bin = resolution - 1 ;
		} else {
			x_bin = (int)((init_x[i]-min_x)/x_bin_size) ;
		}

		num[x_bin]++ ;
		y[x_bin]+=init_y[i] ;
	}
		
	for (int i=0; i<resolution; i++) {
		x[i] = min_x + i*x_bin_size ;
		if (y[i] > 0) {
			y[i] /= num[i] ;
		} else if (i==0) {
			y[i] = 0 ;
		} else {
			y[i] = y[i-1] ;
		}
	}

	// Calculate AUC

	(*auc) = 0 ;
	for (int i=0; i<resolution-1; i++)
		(*auc) += x_bin_size * (y[i] + y[i+1])/2 ;


	free(values) ;
	free(inds) ;
	free(init_x) ;
	free(init_y) ;
	free(tp) ;
	free(fp) ;
	free(x) ;
	free(y) ;
	free(num) ;

	return 0 ;
}
