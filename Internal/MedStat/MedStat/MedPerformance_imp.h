// Implementation of template functions for class MedPerformance.


// Init from array
template <typename T, typename S> MedClassifierPerformance::MedClassifierPerformance(T *preds, S *labels, int n) {
	load(preds,labels,n) ;
}

// Load from arrays
template <typename T, typename S> void MedClassifierPerformance::load(T *_preds, S *_labels, int n) {
	
	preds.resize(1) ;
	preds[0].resize(n) ;
	for (int i=0; i<n; i++) {
		pair<float,float> current((float) _preds[i],(float) _labels[i]) ;
		preds[0][i] = current ;
	}

	init() ;
	ShuffleSort() ;
	getPerformanceValues() ;

	PerformancePointers.resize(preds.size()) ;
}
	
// Init from object
template <typename T> MedClassifierPerformance::MedClassifierPerformance(T& inObj) {
	load(inObj) ;
}

// Load from object
template <typename T> void MedClassifierPerformance::load(T& inObj) {

	_load(inObj) ;	
	post_load();
}




