#define _CRT_SECURE_NO_WARNINGS

#include "classifiers.h"

//R Runners
int R_learn(double *tx, double *y, int nrows, int ncols) {

	// Create Learner Scripts
	FILE *ofp = NULL;
	int ret = 0;

	if ((ofp = safe_fopen(LEARN_BATCH, "w", false)) == NULL) {
		fprintf( stderr,"error (%d) (%s) writing to file %s", errno, strerror(errno), LEARN_BATCH );
		return -1;
	}

	fprintf(ofp,"library(randomForest) ;\n") ;
	fprintf(ofp,"set.seed(%d) ;\n",RF_SEED) ;
	fprintf(ofp,"data <- read.table(\"%s\",header=T) ;\n",LRN_FILE) ;
	fprintf(ofp,"mdl <- randomForest(Label~.,data=data)\n") ;
	fprintf(ofp,"save(mdl,file=\"%s\") ;\n",MDL_FILE) ;
	fclose(ofp) ;



	// Print Data
	if (tprint_matrix(tx,y,nrows,ncols,LRN_FILE) == -1)
		return -1 ;

	// Run script
	char cmd[MAX_STRING_LEN] ;
	sprintf(cmd,"%s %s %s",R_EXEC,LEARN_BATCH,STD_FILE) ;
	if (system(cmd) != 0) {
		fprintf(stderr,"%s failed\n",cmd) ;
		return -1 ;
	}

	return 0 ;
}

void print_rf_learner (FILE *ofp, int eq_size) {

	fprintf(ofp,"library(randomForest) ;\n") ;
	fprintf(ofp,"set.seed(%d) ;\n",RF_SEED) ;
	fprintf(ofp,"data.all <- read.table(\"%s\",header=T) ;\n",LRN_FILE) ;
	fprintf(ofp,"data <- subset(data.all , select = -Label) ;\n") ;
	fprintf(ofp,"labels <- data.all$Label ;\n") ;

	fprintf(ofp,"pos <- which(labels > 0) ;\n") ;
	fprintf(ofp,"neg <- which(labels <= 0) ;\n") ;
	if (eq_size)
		fprintf(ofp,"inds <- c(pos,sample(neg,length(pos))) ;\n") ;
	else
		fprintf(ofp,"inds <- 1:length(labels) ;\n") ;

	fprintf(ofp,"labels.factor <- as.factor(labels) ;\n") ;
	fprintf(ofp,"mdl <- randomForest(data[inds,],labels.factor[inds])\n") ;
	fprintf(ofp,"save(mdl,file=\"%s\") ;\n",MDL_FILE) ;
	fprintf(ofp,"sink(\"%s\") ;\n",TREE_FILE) ;
	fprintf(ofp,"print(mdl$forest$ntree) ;\n") ;
	fprintf(ofp,"for (i in 1:mdl$forest$ntree) {\n") ;
	fprintf(ofp,"	tree <- getTree(mdl,i) ;\n") ;
	fprintf(ofp,"	print(attr(tree,\"dim\")[[1]]) ;\n") ;
	fprintf(ofp,"	print(tree) ;\n") ;
	fprintf(ofp,"}\n") ;
	fprintf(ofp,"sink() ;\n") ;
}

void print_rf_learner (FILE *ofp, char *rf_file, int eq_size) {

	fprintf(ofp,"library(randomForest) ;\n") ;
	fprintf(ofp,"set.seed(%d) ;\n",RF_SEED) ;
	fprintf(ofp,"data.all <- read.table(\"%s\",header=T) ;\n",LRN_FILE) ;
	fprintf(ofp,"data <- subset(data.all , select = -Label) ;\n") ;
	fprintf(ofp,"labels <- data.all$Label ;\n") ;

	fprintf(ofp,"pos <- which(labels > 0) ;\n") ;
	fprintf(ofp,"neg <- which(labels <= 0) ;\n") ;
	if (eq_size)
		fprintf(ofp,"inds <- c(pos,sample(neg,length(pos))) ;\n") ;
	else
		fprintf(ofp,"inds <- 1:length(labels) ;\n") ;

	fprintf(ofp,"labels.factor <- as.factor(labels) ;\n") ;
	fprintf(ofp,"mdl <- randomForest(data[inds,],labels.factor[inds])\n") ;
	fprintf(ofp,"save(mdl,file=\"%s\") ;\n",rf_file) ;
}

void print_lm_learner (FILE *ofp) {

	fprintf(ofp,"data.all <- read.table(\"%s\",header=T) ;\n",LRN_FILE) ;
	fprintf(ofp,"data <- subset(data.all , select = -Label) ;\n") ;
	fprintf(ofp,"labels <- data.all$Label ;\n") ;
	fprintf(ofp,"mdl <- lm(labels ~ .,data);\n") ;
	fprintf(ofp,"save(mdl,file=\"%s\") ;\n",MDL_FILE) ;
}

void print_mars_learner (FILE *ofp) {

	fprintf(ofp,"library(mda) ;\n") ;
	fprintf(ofp,"data.all <- read.table(\"%s\",header=T) ;\n",LRN_FILE) ;
	fprintf(ofp,"data <- subset(data.all , select = -Label) ;\n") ;
	fprintf(ofp,"labels <- data.all$Label ;\n") ;
	fprintf(ofp,"mdl <- mars(data,labels,degree=2);\n") ;
	fprintf(ofp,"save(mdl,file=\"%s\") ;\n",MDL_FILE) ;
}

void print_log_learner (FILE *ofp) {

	fprintf(ofp,"data.all <- read.table(\"%s\",header=T) ;\n",LRN_FILE) ;
	fprintf(ofp,"data <- subset(data.all , select = -Label) ;\n") ;
	fprintf(ofp,"labels <- data.all$Label ;\n") ;
	fprintf(ofp,"mdl <- glm(labels ~ .,data,family=\"binomial\");\n") ;
	fprintf(ofp,"save(mdl,file=\"%s\") ;\n",MDL_FILE) ;
}

int R_get_classifiction_importance(double *tx, double *y, int nrows, int ncols, int *ignore, double *importance) {

	// Create script
	// Create Learner Scripts
	FILE *ofp = NULL;
	int ret = 0;

	if ((ofp = safe_fopen(IMP_BATCH, "w", false)) == NULL) {
		fprintf( stderr,"error (%d) (%s) writing to file %s", errno, strerror(errno), IMP_BATCH );
		return -1;
	}

	fprintf(ofp,"library(randomForest) ;\n") ;
	fprintf(ofp,"set.seed(%d) ;\n",RF_SEED) ;
	fprintf(ofp,"data.all <- read.table(\"%s\",header=T) ;\n",LRN_FILE) ;
	fprintf(ofp,"data <- subset(data.all , select = -Label) ;\n") ;
	fprintf(ofp,"labels <- data.all$Label ;\n") ;

	fprintf(ofp,"pos <- which(labels > 0) ;\n") ;
	fprintf(ofp,"neg <- which(labels < 0) ;\n") ;
	fprintf(ofp,"inds <- c(pos,sample(neg,length(pos))) ;\n") ;

	fprintf(ofp,"labels.factor <- as.factor(labels) ;\n") ;
	fprintf(ofp,"mdl <- randomForest(data[inds,],labels.factor[inds],ntree=3000,importance=T) ;\n") ;
	fprintf(ofp,"imp <- mdl$importance ;\n") ;
	fprintf(ofp,"sink(\"%s\") ;\n",IMP_FILE) ;
	fprintf(ofp,"print(as.matrix(imp[,4])) ;\n") ;
	fprintf(ofp,"sink() ;\n") ;
	fclose (ofp) ;

	// Print Data
	if (tprint_matrix(tx,y,nrows,ncols,ignore,LRN_FILE) == -1)
		return -1 ;

	// Run script
	char cmd[MAX_STRING_LEN] ;
	sprintf(cmd,"%s %s %s",R_EXEC,IMP_BATCH,STD_FILE) ;
	if (system(cmd) != 0) {
		fprintf(stderr,"%s failed\n",cmd) ;
		return -1 ;
	}

	// Read importance
	FILE *ifp = NULL;

	if ((ifp = safe_fopen(IMP_FILE, "r", false)) == NULL) {
		fprintf( stderr,"error (%d) (%s) writing to file %s", errno, strerror(errno), IMP_FILE );
		return -1;
	}

	char dummy[MAX_STRING_LEN] ;
	if (fscanf(ifp,"%s\n",dummy) != 1) {
		fprintf(stderr,"Cannot read header from %s\n",IMP_FILE) ;
		return -1 ;
	}

	double imp ;
	int iftr = 0 ;
	while ((ret = fscanf(ifp,"%s %lf\n",dummy,&imp)) != EOF) {
		if (ret != 2) {
			fprintf(stderr,"Cannot read line %d from %s\n",iftr,IMP_FILE) ;
			return -1 ;
		}
		importance[iftr++] = imp ;
	}
	fclose(ifp) ;

	return 0 ;
}

int R_learn_classification(double *tx, double *y, int nrows, int ncols, char *type, int eq_size) {

	// Create Learner Scripts
	FILE *ofp = NULL;
	int ret = 0;

	if ((ofp = safe_fopen(LEARN_BATCH, "w", false)) == NULL) {
		fprintf( stderr,"error (%d) (%s) writing to file %s", errno, strerror(errno), LEARN_BATCH );
		return -1;
	}

	if (strcmp(type,"rf")==0)
		print_rf_learner(ofp,eq_size) ;
	else if (strcmp(type,"lm")==0)
		print_lm_learner(ofp) ;
	else if (strcmp(type,"log")==0)
		print_log_learner(ofp) ;
	else if (strcmp(type,"mars")==0)
		print_mars_learner(ofp) ;
	else {
		fprintf(stderr,"Unknonw type %s\n",type) ;
		return -1 ;
	}
	fclose(ofp) ;

	// Print Data
	if (tprint_matrix(tx,y,nrows,ncols,LRN_FILE) == -1)
		return -1 ;

	// Run script
	char cmd[MAX_STRING_LEN] ;
	sprintf(cmd,"%s %s %s",R_EXEC,LEARN_BATCH,STD_FILE) ;
	if (system(cmd) != 0) {
		fprintf(stderr,"%s failed\n",cmd) ;
		return -1 ;
	}

	return 0 ;
}

int R_learn_classification(double *tx, double *y, int nrows, int ncols, int *ignore, char *type, int eq_size) {

	// Create Learner Scripts
	FILE *ofp = NULL;
	int ret = 0;

	if ((ofp = safe_fopen(LEARN_BATCH, "w", false)) == NULL) {
		fprintf( stderr,"error (%d) (%s) writing to file %s", errno, strerror(errno), LEARN_BATCH );
		return -1;
	}

	if (strcmp(type,"rf")==0)
		print_rf_learner(ofp,eq_size) ;
	else if (strcmp(type,"lm")==0)
		print_lm_learner(ofp) ;
	else {
		fprintf(stderr,"Unknonw type %s\n",type) ;
		return -1 ;
	}
	fclose(ofp) ;

	// Print Data
	if (tprint_matrix(tx,y,nrows,ncols,ignore,LRN_FILE) == -1)
		return -1 ;

	// Run script

	char cmd[MAX_STRING_LEN] ;
	sprintf(cmd,"%s %s %s",R_EXEC,LEARN_BATCH,STD_FILE) ;
	if (system(cmd) != 0) {
		fprintf(stderr,"%s failed\n",cmd) ;
		return -1 ;
	}

	return 0 ;
}

int R_learn_classification(double *tx, double *y, int nrows, int ncols, char *type, char *rf_file, int eq_size) {

	// Create Learner Scripts
	FILE *ofp = NULL;
	int ret = 0;

	if ((ofp = safe_fopen(LEARN_BATCH, "w", false)) == NULL) {
		fprintf( stderr,"error (%d) (%s) writing to file %s", errno, strerror(errno), LEARN_BATCH );
		return -1;
	}

	if (strcmp(type,"rf")==0)
		print_rf_learner(ofp,rf_file,eq_size) ;
	else if (strcmp(type,"lm")==0)
		print_lm_learner(ofp) ;
	else {
		fprintf(stderr,"Unknonw type %s\n",type) ;
		return -1 ;
	}
	fclose(ofp) ;

	// Print Data
	if (tprint_matrix(tx,y,nrows,ncols,LRN_FILE) == -1)
		return -1 ;

	// Run script

	char cmd[MAX_STRING_LEN] ;
	sprintf(cmd,"%s %s %s",R_EXEC,LEARN_BATCH,STD_FILE) ;
	if (system(cmd) != 0) {
		fprintf(stderr,"%s failed\n",cmd) ;
		return -1 ;
	}

	return 0 ;
}

int R_predict(double *tx, double *preds, int nrows, int ncols) {

	// Create predictor Scripts
	FILE *ofp = NULL;
	int ret = 0;

	if ((ofp = safe_fopen(PRED_BATCH, "w", false)) == NULL) {
		fprintf( stderr,"error (%d) (%s) writing to file %s", errno, strerror(errno), PRED_BATCH );
		return -1;
	}

	fprintf(ofp,"library(randomForest) ;\n") ;
	fprintf(ofp,"set.seed(%d) ;\n",RF_SEED) ;
	fprintf(ofp,"data <- read.table(\"%s\",header=T) ;\n",PRD_FILE) ;
	fprintf(ofp,"load(\"%s\") ;\n",MDL_FILE) ;
	fprintf(ofp,"preds <- predict(mdl,data) ;\n") ;
	fprintf(ofp,"sink(\"%s\") ;\n",PRED_FILE) ;
	fprintf(ofp,"print(as.matrix(preds,ncol=1)) ;\n") ;
	fprintf(ofp,"sink() ;\n") ;
	fclose(ofp) ;

	if (file_exists(STD_FILE) && remove(STD_FILE) != 0) {
		fprintf(stderr,"Cannot remove %s\n",STD_FILE) ;
		return -1 ;
	}

	// Print data
	if (tprint_matrix(tx ,nrows,ncols,PRD_FILE) == -1)
		return -1 ;

	// Run script
	char cmd[MAX_STRING_LEN] ;
	sprintf(cmd,"%s %s %s",R_EXEC,PRED_BATCH,STD_FILE) ;
	if (system(cmd) != 0) {
		fprintf(stderr,"%s failed\n",cmd) ;
		return -1 ;
	}

	// Read predictions
	FILE *fp = safe_fopen(PRED_FILE, "rb", false) ;
	if (fp == NULL) {
		fprintf(stderr,"Cannot open %s for reading",PRED_FILE) ;
		return -1 ;
	}

	char head[MAX_STRING_LEN] ;
	if (fscanf(fp,"%s",head)!=1) {
		fprintf(stderr,"Cannot read header\n") ;
		return -1 ;
	}

	int idx ;
	double pred ;
	for (int i=0; i<nrows; i++) {
		if (fscanf(fp,"%d %lf",&idx,&pred)!=2 || idx!=i+1) {
			fprintf(stderr,"Cannot read line %d\n",i) ;
			return -1 ;
		}
		preds[i] = pred ;
	}

	return 0 ;
}

void print_rf_predictor (FILE *ofp) {	
	fprintf(ofp,"library(randomForest) ;\n") ;
	fprintf(ofp,"set.seed(%d) ;\n",RF_SEED) ;
	fprintf(ofp,"data <- read.table(\"%s\",header=T) ;\n",PRD_FILE) ;
	fprintf(ofp,"load(\"%s\") ;\n",MDL_FILE) ;
	fprintf(ofp,"probs <- predict(mdl,data,type=\"prob\") ;\n") ;
	fprintf(ofp,"sink(\"%s\") ;\n",PRED_FILE) ;
	fprintf(ofp,"print(probs) ;\n") ;
	fprintf(ofp,"sink() ;\n") ;
}

void print_rf_predictor (FILE *ofp, char *mdl_file) {	
	fprintf(ofp,"library(randomForest) ;\n") ;
	fprintf(ofp,"set.seed(%d) ;\n",RF_SEED) ;
	fprintf(ofp,"data <- read.table(\"%s\",header=T) ;\n",PRD_FILE) ;
	fprintf(ofp,"load(\"%s\") ;\n",mdl_file) ;
	fprintf(ofp,"probs <- predict(mdl,data,type=\"prob\") ;\n") ;
	fprintf(ofp,"sink(\"%s\") ;\n",PRED_FILE) ;
	fprintf(ofp,"print(probs) ;\n") ;
	fprintf(ofp,"sink() ;\n") ;
}

void print_lm_predictor (FILE *ofp) {

	fprintf(ofp,"data <- read.table(\"%s\",header=T) ;\n",PRD_FILE) ;
	fprintf(ofp,"load(\"%s\") ;\n",MDL_FILE) ;
	fprintf(ofp,"preds <- predict(mdl,data) ;\n") ;
	fprintf(ofp,"sink(\"%s\") ;\n",PRED_FILE) ;
	fprintf(ofp,"print(as.matrix(preds)) ;\n") ;
	fprintf(ofp,"sink() ;\n") ;
}

void print_lm_predictor (FILE *ofp, char *mdl_file) {

	fprintf(ofp,"data <- read.table(\"%s\",header=T) ;\n",PRD_FILE) ;
	fprintf(ofp,"load(\"%s\") ;\n",mdl_file) ;
	fprintf(ofp,"preds <- predict(mdl,data) ;\n") ;
	fprintf(ofp,"sink(\"%s\") ;\n",PRED_FILE) ;
	fprintf(ofp,"print(as.matrix(preds)) ;\n") ;
	fprintf(ofp,"sink() ;\n") ;
}

void print_mars_predictor (FILE *ofp) {

	fprintf(ofp,"library(mda) ;\n") ;
	fprintf(ofp,"data <- read.table(\"%s\",header=T) ;\n",PRD_FILE) ;
	fprintf(ofp,"load(\"%s\") ;\n",MDL_FILE) ;
	fprintf(ofp,"preds <- predict(mdl,data) ;\n") ;
	fprintf(ofp,"sink(\"%s\") ;\n",PRED_FILE) ;
	fprintf(ofp,"print(as.matrix(preds)) ;\n") ;
	fprintf(ofp,"sink() ;\n") ;
}

void print_mars_predictor (FILE *ofp, char *mdl_file) {

	fprintf(ofp,"library(mda) ;\n") ;
	fprintf(ofp,"data <- read.table(\"%s\",header=T) ;\n",PRD_FILE) ;
	fprintf(ofp,"load(\"%s\") ;\n",mdl_file) ;
	fprintf(ofp,"preds <- predict(mdl,data) ;\n") ;
	fprintf(ofp,"sink(\"%s\") ;\n",PRED_FILE) ;
	fprintf(ofp,"print(as.matrix(preds)) ;\n") ;
	fprintf(ofp,"sink() ;\n") ;
}

int rf_predict(double *x, double *preds, int nrows, int ncols, random_forest *rf) {
	return rf_predict(x,preds,nrows,ncols,1,rf) ;
}

int rf_predict(double *tx, double *preds, int nrows, int ncols, int trans_flag, random_forest *rf) {

	// Transpose matrix ; 
	double *x  ;

	if (trans_flag) {
		if (transpose(tx,&x,nrows,ncols,1)==-1) {
			fprintf(stderr,"Transpose failed\n") ;
			return -1 ;
		}
	} else
		x = tx ;

	// Predict
	int counts[2] ;
	for (int i=0; i<nrows; i++) {
		if (i%10000 == 0)
			fprintf(stderr,"Predicting %d/%d ...",i+1,nrows) ;
		forest_pred(rf,x + i*ncols, counts, 2) ;
		preds[i] = counts[1] / (rf->ntrees + 0.0) ;
	}
	fprintf(stderr,"\n") ;

	if (trans_flag)
		free(x) ;

	return 0;
}

int rf_predict_from_text_forest(double *tx, double *preds, int nrows, int ncols, char *ftree) {

	// Read forest
	random_forest rf ;
	if (read_text_forest(ftree,&rf)==-1) {
		fprintf(stderr,"RF reading failed\n") ;
		return -1 ;
	}

	return rf_predict(tx,preds,nrows,ncols,&rf) ;

}

int rf_predict_from_text_forest(double *tx, double *preds, int nrows, int ncols) {
	return rf_predict_from_text_forest(tx,preds,nrows,ncols,TREE_FILE) ;
}

int rf_predict(double *tx, double *preds, int nrows, int ncols, char *ftree) {

	// Read forest
	random_forest rf ;
	if (read_forest(ftree,&rf)==-1) {
		fprintf(stderr,"RF reading failed\n") ;
		return -1 ;
	}

	return rf_predict(tx,preds,nrows,ncols,&rf) ;
}

int rf_predict(double *tx, double *preds, int nrows, int ncols) {

	return rf_predict(tx,preds,nrows,ncols,TREE_FILE) ;
}

int rf_predict(double *tx, double *preds, int nrows, int ncols, unsigned char *model) {

	int size ;

	// Extract forest
	random_forest rf ;
	if ((size = rf_deserialize(model,&rf)) == -1) {
		fprintf(stderr,"Cannot DeSerialize RF model\n") ;
		return -1 ;
	}

	if (rf_predict(tx,preds,nrows,ncols,&rf) == -1)
		return -1 ;
	else 
		return size ;
}

int get_prediction_matrix(double *tx, int nrows, int ncols, double **preds, int *ntrees) {

	// Transpose matrix ; 
	double *x  ;
	if (transpose(tx,&x,nrows,ncols,1)==-1) {
		fprintf(stderr,"Transpose failed\n") ;
		return -1 ;
	}

	// Read forest
	random_forest rf ;
	if (read_forest(TREE_FILE,&rf)==-1) {
		fprintf(stderr,"RF reading failed\n") ;
		return -1 ;
	}

	// Allocate
	if (((*preds) = (double *) malloc((*ntrees)*nrows*sizeof(double)))==NULL) {
		fprintf(stderr,"Allocation failed\n") ;
		return -1 ;
	}

	// Predict
	for (int i=0; i<nrows; i++) {
		for (int j=0; j<(*ntrees); j++) {
			(*preds)[XIDX(i,j,*ntrees)] = (double) tree_pred(rf.trees[j],x + i*ncols) ;
		}
	}

	free(x) ;
	return 0 ;
}

int R_predict_class(double *tx, double *preds, int nrows, int ncols, char *type) {

	// Create predictor Scripts
	FILE *ofp = NULL;
	int ret = 0;

	if ((ofp = safe_fopen(PRED_BATCH, "w", false)) == NULL) {
		fprintf( stderr,"error (%d) (%s) writing to file %s", errno, strerror(errno), PRED_BATCH );
		return -1;
	}

	if (strcmp(type,"rf")==0)
		print_rf_predictor(ofp) ;
	else if (strcmp(type,"lm")==0 || strcmp(type,"log")==0)
		print_lm_predictor(ofp) ;
	else if (strcmp(type,"mars")==0)
		print_mars_predictor(ofp) ;
	else {
		fprintf(stderr,"Unknonw type %s\n",type) ;
		return -1 ;
	}
	fclose(ofp) ;

	if (file_exists(STD_FILE) && remove(STD_FILE) != 0) {
		fprintf(stderr,"Cannot remove %s\n",STD_FILE) ;
		return -1 ;
	}

	int first_line = 0 ;
	while (first_line < nrows) {
		int last_line = (first_line + MAX_R_ROWS > nrows ) ? (nrows-1) : (first_line + MAX_R_ROWS - 1) ;
		fprintf(stderr,"predicting %d - %d / %d\n",first_line,last_line,nrows) ;

		// Print Data
		if (tprint_matrix(tx,nrows,first_line,last_line,ncols,PRD_FILE) == -1)
			return -1 ;
	
		// Run Script
		char cmd[MAX_STRING_LEN] ;
		sprintf(cmd,"%s %s %s",R_EXEC,PRED_BATCH,STD_FILE) ;
		if (system(cmd) != 0) {
			fprintf(stderr,"%s failed\n",cmd) ;
			return -1 ;
		}

		// Read predictions
		FILE *fp = safe_fopen(PRED_FILE, "rb", false) ;
		if (fp == NULL) {
			fprintf(stderr,"Cannot open %s for reading",PRED_FILE) ;
			return -1 ;
		}

		double dmy1,dmy2 ;
		char dmy[MAX_STRING_LEN] ;
		if ((strcmp(type,"rf")==0 && fscanf(fp,"%lf %lf",&dmy1,&dmy2)!=2) ||
			(strcmp(type,"lm")==0 && fscanf(fp,"%s",dmy)!=1) ||
			(strcmp(type,"mars")==0 && fscanf(fp,"%s",dmy)!=1)) {
			fprintf(stderr,"Cannot read header\n") ;
			return -1 ;
		}

		int idx ;
		double p1,p ;
		for (int i=first_line; i<=last_line; i++) {
			if ((strcmp(type,"rf")==0 && fscanf(fp,"%d %lf %lf",&idx,&p1,&p)!=3) || 
				(strcmp(type,"lm")==0 && fscanf(fp,"%d %lf",&idx,&p)!=2)|| 
				(strcmp(type,"mars")==0 && fscanf(fp,"%s %lf",dmy,&p)!=2) ) {
				fprintf(stderr,"Cannot read line %d \n",i) ;
				return -1 ;
			}
			preds[i] = p ;
		}
		fclose(fp) ;

		first_line = last_line + 1 ;
	}

	return 0 ;
}

int R_predict_class(double *tx, double *preds, int nrows, int ncols, char *type, char *mdl_file) {

	// Create predictor Scripts
	FILE *ofp = NULL;
	int ret = 0;

	if ((ofp = safe_fopen(PRED_BATCH, "w", false)) == NULL) {
		fprintf( stderr,"error (%d) (%s) writing to file %s", errno, strerror(errno), PRED_BATCH );
		return -1;
	}

	if (strcmp(type,"rf")==0)
		print_rf_predictor(ofp,mdl_file) ;
	else if (strcmp(type,"lm")==0 || strcmp(type,"log")==0)
		print_lm_predictor(ofp,mdl_file) ;
	else {
		fprintf(stderr,"Unknonw type %s\n",type) ;
		return -1 ;
	}
	fclose(ofp) ;

#ifdef _WIN32
	unsigned long id = GetCurrentProcessId();
#else
	unsigned long id = getpid();
#endif
	time_t my_time;
	time(&my_time);

	char stderr_file[MAX_STRING_LEN] ;
#ifdef _WIN32
	sprintf(stderr_file,"%s\\%s.%lld.%ld",STD_DIR,STD_FILE,my_time,id) ;
#else
	sprintf(stderr_file,"%s/%s.%d.%ld",STD_DIR,STD_FILE,my_time,id) ;
#endif

	int first_line = 0 ;
	while (first_line < nrows) {
		int last_line = (first_line + MAX_R_ROWS > nrows ) ? (nrows-1) : (first_line + MAX_R_ROWS - 1) ;
		fprintf(stderr,"predicting %d - %d / %d\n",first_line,last_line,nrows) ;

		// Print Data
		if (tprint_matrix(tx,nrows,first_line,last_line,ncols,PRD_FILE) == -1)
			return -1 ;
	
		// Run Script
		char cmd[MAX_STRING_LEN] ;
		sprintf(cmd,"%s %s %s",R_EXEC,PRED_BATCH,stderr_file) ;
		if (system(cmd) != 0) {
			fprintf(stderr,"%s failed\n",cmd) ;
			return -1 ;
		}

		// Read predictions
		FILE *fp = safe_fopen(PRED_FILE, "rb", false) ;
		if (fp == NULL) {
			fprintf(stderr,"Cannot open %s for reading",PRED_FILE) ;
			return -1 ;
		}

		double dmy1,dmy2 ;
		char dmy[MAX_STRING_LEN] ;
		if ((strcmp(type,"rf")==0 && fscanf(fp,"%lf %lf",&dmy1,&dmy2)!=2) ||
			(strcmp(type,"lm")==0 && fscanf(fp,"%s",dmy)!=1) ||
			(strcmp(type,"mars")==0 && fscanf(fp,"%s",dmy)!=1)) {
			fprintf(stderr,"Cannot read header\n") ;
			return -1 ;
		}

		int idx ;
		double p1,p ;
		for (int i=first_line; i<=last_line; i++) {
			if ((strcmp(type,"rf")==0 && fscanf(fp,"%d %lf %lf",&idx,&p1,&p)!=3) || 
				(strcmp(type,"lm")==0 && fscanf(fp,"%d %lf",&idx,&p)!=2)|| 
				(strcmp(type,"mars")==0 && fscanf(fp,"%s %lf",dmy,&p)!=2) ) {
				return -1 ;
			}
			preds[i] = p ;
		}
		fclose(fp) ;

		first_line = last_line + 1 ;
	}

	return 0 ;
}

int R_predict_class_with_NA(double *tx, double *preds, int nrows, int ncols, char *type, double missing) {

	// Create predictor Scripts
	FILE *ofp = NULL;
	int ret = 0;

	if ((ofp = safe_fopen(PRED_BATCH, "w", false)) == NULL) {
		fprintf( stderr,"error (%d) (%s) writing to file %s", errno, strerror(errno), PRED_BATCH );
		return -1;
	}

	if (strcmp(type,"rf")==0)
		print_rf_predictor(ofp) ;
	else if (strcmp(type,"lm")==0)
		print_lm_predictor(ofp) ;
	else {
		fprintf(stderr,"Unknonw type %s\n",type) ;
		return -1 ;
	}
	fclose(ofp) ;

	if (file_exists(STD_FILE) && remove(STD_FILE) != 0) {
		fprintf(stderr,"Cannot remove %s\n",STD_FILE) ;
		return -1 ;
	}

	// Print Data
	if (tprint_matrix(tx,nrows,ncols,PRD_FILE,missing) == -1)
		return -1 ;
	
	// Run Script
	char cmd[MAX_STRING_LEN] ;
	sprintf(cmd,"%s %s %s",R_EXEC,PRED_BATCH,STD_FILE) ;
	if (system(cmd) != 0) {
		fprintf(stderr,"%s failed\n",cmd) ;
		return -1 ;
	}

	// Read predictions
	FILE *fp = safe_fopen(PRED_FILE, "rb", false) ;
	if (fp == NULL) {
		fprintf(stderr,"Cannot open %s for reading",PRED_FILE) ;
		return -1 ;
	}

	double dmy1,dmy2 ;
	char dmy[MAX_STRING_LEN] ;
	if ((strcmp(type,"rf")==0 && fscanf(fp,"%lf %lf",&dmy1,&dmy2)!=2) ||
		(strcmp(type,"lm")==0 && fscanf(fp,"%s",dmy)!=1) ||
		(strcmp(type,"mars")==0 && fscanf(fp,"%s",dmy)!=1)) {
		fprintf(stderr,"Cannot read header\n") ;
		return -1 ;
	}

	int idx ;
	double p1,p ;
	for (int i=0; i<nrows; i++) {
		if ((strcmp(type,"rf")==0 && fscanf(fp,"%d %lf %lf",&idx,&p1,&p)!=3) || 
			(strcmp(type,"lm")==0 && fscanf(fp,"%d %lf",&idx,&p)!=2)|| 
			(strcmp(type,"mars")==0 && fscanf(fp,"%s %lf",dmy,&p)!=2) ) {
				fprintf(stderr,"Cannot read line %d\n",i) ;
				return -1 ;
		}
		preds[i] = p ;
	}
	fclose(fp) ;

	return 0 ;
}

int R_predict_class(double *tx, double *preds, int nrows, int ncols, int *ignore, char *type) {

	// Create predictor Scripts
	FILE *ofp = NULL;
	int ret = 0;

	if ((ofp = safe_fopen(PRED_BATCH, "w", false)) == NULL) {
		fprintf( stderr,"error (%d) (%s) writing to file %s", errno, strerror(errno), PRED_BATCH );
		return -1;
	}

	if (strcmp(type,"rf")==0)
		print_rf_predictor(ofp) ;
	else if (strcmp(type,"lm")==0)
		print_lm_predictor(ofp) ;
	else {
		fprintf(stderr,"Unknonw type %s\n",type) ;
		return -1 ;
	}
	fclose(ofp) ;

	if (file_exists(STD_FILE) && remove(STD_FILE) != 0) {
		fprintf(stderr,"Cannot remove %s\n",STD_FILE) ;
		return -1 ;
	}

	// Print Data
	if (tprint_matrix(tx,nrows,ncols,PRD_FILE) == -1)
		return -1 ;
	
	// Run Script
	char cmd[MAX_STRING_LEN] ;
	sprintf(cmd,"%s %s %s",R_EXEC,PRED_BATCH,STD_FILE) ;
	if (system(cmd) != 0) {
		fprintf(stderr,"%s failed\n",cmd) ;
		return -1 ;
	}

	// Read predictions
	FILE *fp = safe_fopen(PRED_FILE, "rb", false) ;
	if (fp == NULL) {
		fprintf(stderr,"Cannot open %s for reading",PRED_FILE) ;
		return -1 ;
	}

	double dmy1,dmy2 ;
	char dmy[MAX_STRING_LEN] ;
	if ((strcmp(type,"rf")==0 && fscanf(fp,"%lf %lf",&dmy1,&dmy2)!=2) ||
		(strcmp(type,"lm")==0 && fscanf(fp,"%s",dmy)!=1) ||
		(strcmp(type,"mars")==0 && fscanf(fp,"%s",dmy)!=1)) {
		fprintf(stderr,"Cannot read header\n") ;
		return -1 ;
	}

	int idx ;
	double p1,p ;
	for (int i=0; i<nrows; i++) {
		if ((strcmp(type,"rf")==0 && fscanf(fp,"%d %lf %lf",&idx,&p1,&p)!=3) || 
			(strcmp(type,"lm")==0 && fscanf(fp,"%d %lf",&idx,&p)!=2)|| 
			(strcmp(type,"mars")==0 && fscanf(fp,"%s %lf",dmy,&p)!=2) ) {
				fprintf(stderr,"Cannot read line %d\n",i) ;
				return -1 ;
		}
		preds[i] = p ;
	}
	fclose(fp) ;

	return 0 ;
}

int R_predict_class_with_NA(double *tx, double *preds, int nrows, int ncols, int *ignore, char *type, double missing) {

	// Create predictor Scripts
	FILE *ofp = NULL;
	int ret = 0;

	if ((ofp = safe_fopen(PRED_BATCH, "w", false)) == NULL) {
		fprintf( stderr,"error (%d) (%s) writing to file %s", errno, strerror(errno), PRED_BATCH );
		return -1;
	}

	if (strcmp(type,"rf")==0)
		print_rf_predictor(ofp) ;
	else if (strcmp(type,"lm")==0)
		print_lm_predictor(ofp) ;
	else {
		fprintf(stderr,"Unknonw type %s\n",type) ;
		return -1 ;
	}
	fclose(ofp) ;

	if (file_exists(STD_FILE) && remove(STD_FILE) != 0) {
		fprintf(stderr,"Cannot remove %s\n",STD_FILE) ;
		return -1 ;
	}

	// Print Data
	if (tprint_matrix(tx,nrows,ncols,PRD_FILE,missing) == -1)
		return -1 ;
	
	// Run Script
	char cmd[MAX_STRING_LEN] ;
	sprintf(cmd,"%s %s %s",R_EXEC,PRED_BATCH,STD_FILE) ;
	if (system(cmd) != 0) {
		fprintf(stderr,"%s failed\n",cmd) ;
		return -1 ;
	}

	// Read predictions
	FILE *fp = safe_fopen(PRED_FILE, "rb", false) ;
	if (fp == NULL) {
		fprintf(stderr,"Cannot open %s for reading",PRED_FILE) ;
		return -1 ;
	}

	double dmy1,dmy2 ;
	char dmy[MAX_STRING_LEN] ;
	if ((strcmp(type,"rf")==0 && fscanf(fp,"%lf %lf",&dmy1,&dmy2)!=2) ||
		(strcmp(type,"lm")==0 && fscanf(fp,"%s",dmy)!=1) ||
		(strcmp(type,"mars")==0 && fscanf(fp,"%s",dmy)!=1)) {
		fprintf(stderr,"Cannot read header\n") ;
		return -1 ;
	}

	int idx ;
	double p1,p ;
	for (int i=0; i<nrows; i++) {
		if ((strcmp(type,"rf")==0 && fscanf(fp,"%d %lf %lf",&idx,&p1,&p)!=3) || 
			(strcmp(type,"lm")==0 && fscanf(fp,"%d %lf",&idx,&p)!=2)|| 
			(strcmp(type,"mars")==0 && fscanf(fp,"%s %lf",dmy,&p)!=2) ) {
				fprintf(stderr,"Cannot read line %d\n",i) ;
				return -1 ;
		}
		preds[i] = p ;
	}
	fclose(fp) ;

	return 0 ;
}

int R_get_oob_err (double *tx, double *y, int nrows, int ncols, int *ignore, int ntrees, double *err) {

	// Print matrix
	if (tprint_matrix(tx,y,nrows,ncols,ignore,OOB_FILE) == -1)
		return -1 ;

	// Create script
	FILE *ofp = NULL;
	int ret = 0;

	if ((ofp = safe_fopen(OOB_BATCH, "w", false)) == NULL) {
		fprintf( stderr,"error (%d) (%s) at writing to file %s\n", errno, strerror(errno), OOB_BATCH );
		return -1;
	}

	fprintf(ofp,"library(randomForest) ;\n") ;
	fprintf(ofp,"set.seed(%d) ;\n",RF_SEED) ;
	fprintf(ofp,"data.all <- read.table(\"%s\",header=T) ;\n",OOB_FILE) ;
	fprintf(ofp,"data <- subset(data.all , select = -Label) ;\n") ;
	fprintf(ofp,"labels <- data.all$Label ;\n") ;
	fprintf(ofp,"labels.factor <- as.factor(labels) ;\n") ;
	fprintf(ofp,"mdl <- randomForest(data,labels.factor,ntree=%d) ;\n",ntrees) ;
	fprintf(ofp,"sink(\"%s\");\n",ERR_FILE) ;
	fprintf(ofp,"print(mdl$err.rate[%d,1]); \n",ntrees) ;
	fprintf(ofp,"sink() ; \n") ;
	fclose(ofp) ;

	// Run script
	char cmd[MAX_STRING_LEN] ;
	sprintf(cmd,"%s %s %s",R_EXEC,OOB_BATCH,STD_FILE) ;
	if (system(cmd) != 0) {
		fprintf(stderr,"%s failed\n",cmd) ;
		return -1 ;
	}

	// Read result
	FILE *ifp = NULL ;
	if ((ifp = safe_fopen(ERR_FILE, "r", false)) == NULL) {
		fprintf( stderr,"error (%d) (%s) reading file %s\n", errno, strerror(errno), ERR_FILE );
		return -1;
	}

	char dummy[MAX_STRING_LEN] ;
	if (fscanf(ifp,"%s",dummy) != 1) {
		fprintf(stderr,"Problems reading header of %s\n",ERR_FILE) ;
		return -1 ;
	}

	if (fscanf(ifp,"%lf",err) != 1) {
		fprintf(stderr,"Problems reading error from %s\n",ERR_FILE) ;
		return -1 ;
	}
	fclose(ifp) ;

	return 0 ;
}