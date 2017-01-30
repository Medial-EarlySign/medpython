// Management of trees and predictions.

#include "stdafx.h"
#include "gbm.h"
#include "gbm_structs.h"

#define MISSING_VAL -999999.9

// Prediction
GBMRESULT gbm_pred (double *radX, int cRows, int cCols, int cTrees, double rdInitF, gbm_tree *rTrees, int **rCSplits, int *raiVarType, double *radPredF, double missing_val) {
   	
   unsigned long hr = 0;
   int iTree = 0;
   int iObs = 0;

   int *aiSplitVar = NULL;
   double *adSplitCode = NULL;
   int *aiLeftNode = NULL;
   int *aiRightNode = NULL;
   int *aiMissingNode = NULL;
   int iCurrentNode = 0;
   double dX = 0.0;
   int iCatSplitIndicator = 0;
 
   // initialize with the intercept for only the smallest rcTrees
   for(iObs=0; iObs<cRows; iObs++)      
	   radPredF[iObs] = rdInitF;

   iTree = 0; 
 
   while(iTree<cTrees){ 
	   aiSplitVar    = rTrees[iTree].split_var ; 
	   adSplitCode   = rTrees[iTree].split_point ;
	   aiLeftNode    = rTrees[iTree].left_node;
	   aiRightNode   = rTrees[iTree].right_node;
	   aiMissingNode = rTrees[iTree].missing_node;

	   for(iObs=0; iObs<cRows; iObs++) {
		   iCurrentNode = 0; 
		   while(aiSplitVar[iCurrentNode] != -1) {
			   dX = radX[aiSplitVar[iCurrentNode]*cRows + iObs];
           
			   if(dX==missing_val) { // missing ?
				   iCurrentNode = aiMissingNode[iCurrentNode];
			   } else if (raiVarType[aiSplitVar[iCurrentNode]] == 0) { // continuous?
				   if(dX < adSplitCode[iCurrentNode])
					   iCurrentNode = aiLeftNode[iCurrentNode];
				   else
					   iCurrentNode = aiRightNode[iCurrentNode] ;
			   } else { // categorical
				   iCatSplitIndicator = rCSplits[iCurrentNode][(int)dX];
				   if(iCatSplitIndicator==-1)
					   iCurrentNode = aiLeftNode[iCurrentNode]; 
				   else if(iCatSplitIndicator==1)
					   iCurrentNode = aiRightNode[iCurrentNode]; 
				   else // categorical level not present in training
					   iCurrentNode = aiMissingNode[iCurrentNode];
			   }
		   }

		   radPredF[iObs] += adSplitCode[iCurrentNode]; // add the prediction
	   } // iObs

	   iTree++;
   } // iTree
    

   return GBM_OK;

}

GBMRESULT gbm_pred (double *radX, int cRows, int cCols, int cTrees, double rdInitF, gbm_tree *rTrees, int **rCSplits, int *raiVarType, double *radPredF) {
	return gbm_pred(radX,cRows,cCols,cTrees,rdInitF,rTrees,rCSplits,raiVarType,radPredF,MISSING_VAL) ;
}

GBMRESULT gbm_pred (double *radX, int iObs, int cRows, int rcCols, gbm_tree *rTrees, int itree, int **rCSplits, int *raiVarType, double *radPredF, 
					double missing_val)
{
 
   unsigned long hr = 0;

   int *aiSplitVar = NULL;
   double *adSplitCode = NULL;
   int *aiLeftNode = NULL;
   int *aiRightNode = NULL;
   int *aiMissingNode = NULL;
   int iCurrentNode = 0;
   double dX = 0.0;
   int iCatSplitIndicator = 0;
 
   aiSplitVar    = rTrees[itree].split_var ; 
   adSplitCode   = rTrees[itree].split_point ;
   aiLeftNode    = rTrees[itree].left_node;
   aiRightNode   = rTrees[itree].right_node;
   aiMissingNode = rTrees[itree].missing_node;
   
   iCurrentNode = 0;
   while(aiSplitVar[iCurrentNode] != -1) {
       dX = radX[aiSplitVar[iCurrentNode]*cRows + iObs];
      
	   if(dX==missing_val) {  // missing?
          iCurrentNode = aiMissingNode[iCurrentNode];
	   }else if(raiVarType[aiSplitVar[iCurrentNode]] == 0) {  // continuous?
		   if(dX < adSplitCode[iCurrentNode])               
			   iCurrentNode = aiLeftNode[iCurrentNode];
		   else
			   iCurrentNode = aiRightNode[iCurrentNode];
	   } else { // categorical         
		   iCatSplitIndicator = rCSplits[iCurrentNode][(int)dX];     
		   if(iCatSplitIndicator==-1)
			   iCurrentNode = aiLeftNode[iCurrentNode];
		   else if(iCatSplitIndicator==1)
			   iCurrentNode = aiRightNode[iCurrentNode];
		   else // categorical level not present in training          
			   iCurrentNode = aiMissingNode[iCurrentNode] ;
       }
    }

    *radPredF = adSplitCode[iCurrentNode]; // add the prediction       
	return GBM_OK;
}

GBMRESULT gbm_pred (double *radX, int iObs, int cRows, int rcCols, gbm_tree *rTrees, int itree, int **rCSplits, int *raiVarType, double *radPredF) {
	return gbm_pred(radX,iObs,cRows,rcCols,rTrees,itree,rCSplits,raiVarType,radPredF,MISSING_VAL) ;
}


// Keep a tree in memory
gbm_tree *allocate_gbm_trees(int ntrees) {

	return (gbm_tree *) malloc (ntrees*sizeof(gbm_tree)) ;

}

void free_gbm_trees(gbm_tree *trees) {
	free(trees) ;
}

GBMRESULT set_tree(CGBM *pGBM, int nnodes, VEC_VEC_CATEGORIES& vecSplitCodes, gbm_tree *trees, int itree) {

	trees[itree].nnodes = nnodes ;
	trees[itree].split_var = (int *) malloc(nnodes*sizeof(int)) ;
	trees[itree].split_point = (double *) malloc(nnodes*sizeof(double)) ;
	trees[itree].left_node = (int *) malloc(nnodes*sizeof(int)) ;
	trees[itree].right_node = (int *) malloc(nnodes*sizeof(int)) ;
	trees[itree].missing_node = (int *) malloc(nnodes*sizeof(int)) ;
	trees[itree].error_reduction = (double *) malloc(nnodes*sizeof(double)) ;
	trees[itree].weight = (double *) malloc(nnodes*sizeof(double)) ;
	trees[itree].pred = (double *) malloc(nnodes*sizeof(double)) ;

	if (trees[itree].split_var == NULL || trees[itree].split_point == NULL || trees[itree].left_node == NULL || trees[itree].right_node == NULL || 
		trees[itree].missing_node == NULL || trees[itree].error_reduction == NULL || trees[itree].weight == NULL || trees[itree].pred == NULL)
		return GBM_OUTOFMEMORY;

	
    int hr = gbm_transfer_to_R(pGBM,vecSplitCodes,trees[itree].split_var,trees[itree].split_point,trees[itree].left_node,trees[itree].right_node,trees[itree].missing_node,
		trees[itree].error_reduction,trees[itree].weight,trees[itree].pred,0) ;

	return hr ;
}

int read_gbm_tree(FILE *file, gbm_tree &tree) {

	int nnodes ;
	if (fscanf(file,"Nnodes = %d\n",&nnodes) != 1) {
		fprintf(stderr,"Cannot read Nnodes\n") ;
		return -1 ;
	}
	tree.nnodes = nnodes ;

	tree.split_var = (int *) malloc(nnodes*sizeof(int)) ;
	tree.split_point = (double *) malloc(nnodes*sizeof(double)) ;
	tree.left_node = (int *) malloc(nnodes*sizeof(int)) ;
	tree.right_node = (int *) malloc(nnodes*sizeof(int)) ;
	tree.missing_node = (int *) malloc(nnodes*sizeof(int)) ;
	tree.error_reduction = (double *) malloc(nnodes*sizeof(double)) ;
	tree.weight = (double *) malloc(nnodes*sizeof(double)) ;
	tree.pred = (double *) malloc(nnodes*sizeof(double)) ;

	if (tree.split_var == NULL || tree.split_point == NULL || tree.left_node == NULL || tree.right_node == NULL || tree.missing_node == NULL || tree.error_reduction == NULL ||
		tree.weight == NULL || tree.pred == NULL)
		return -1 ;

	int idx,ival1,ival2,ival3 ;
	double rval ;
	for (int i=0; i<tree.nnodes; i++) {
		if (fscanf(file,"Node %d %d %lf %d %d\n",&idx,&ival1,&rval,&ival2,&ival3)!=5 || idx != i) {
			fprintf(stderr,"Problem reading node %d\n",i) ;
			return -1 ;
		}
		tree.split_var[i] = ival1 ;
		tree.split_point[i] = rval ;
		tree.left_node[i]= ival2 ;
		tree.right_node[i]= ival3 ;
	}

	return 0 ;
}

int read_gbm_tree(FILE *file, gbm_tree *trees, int i) {
	return read_gbm_tree(file,trees[i]) ;
}

void print_gbm_tree(FILE *file, gbm_tree &tree) {

	fprintf(file,"Nnodes = %d\n",tree.nnodes) ;
	for (int i=0; i<tree.nnodes; i++)
		fprintf(file,"Node %d %d %f %d %d\n",i,tree.split_var[i],tree.split_point[i],tree.left_node[i],tree.right_node[i]) ;

	return ;
}

void print_gbm_tree(FILE *file, gbm_tree *trees, int i) {
	print_gbm_tree(file,trees[i]) ;
}

void print_gbm_tree(const string& prefix, FILE *file, gbm_tree &tree) {

	fprintf(file,"%s: Nnodes = %d\n",prefix.c_str(),tree.nnodes) ;
	for (int i=0; i<tree.nnodes; i++)
		fprintf(file,"%s: Node %d %d %f %d %d\n",prefix.c_str(),i,tree.split_var[i],tree.split_point[i],tree.left_node[i],tree.right_node[i]) ;

	return ;
}

void print_gbm_tree(const string& prefix, FILE *file, gbm_tree *trees, int i) {
	print_gbm_tree(prefix,file,trees[i]) ;
}