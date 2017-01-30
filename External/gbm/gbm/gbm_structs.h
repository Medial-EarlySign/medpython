struct gbm_tree {
	int nnodes ;
	int *split_var ;
	double *split_point ;
	int *left_node ;
	int *right_node ;
	int *missing_node ;
	double *error_reduction ;
	double *weight ;
	double *pred ;
} ;

struct cat_splits {
	int n ;
	int **splits ;
} ;