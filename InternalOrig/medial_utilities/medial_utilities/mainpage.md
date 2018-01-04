@mainpage Welcome To Medial Libs Documentation
you can also use the search box to find your file/function

 @section internal Internal Libs
 
 @subsection InfraMed InfraMed
 * MedRepository or MedPidRepository - the main object to fetch data from our repositories

 @subsection MedAlgo MedAlgo
 * MedPredictor - basic predictor interface. you may see all derived classes and models in the link

 @subsection MedProcessTools MedProcessTools
 * MedModel - medial pipeline model for running all process. Pipeline objects:
    * RepProcessor - processing repository. to see all options for json file please reffer to ::RepProcessorTypes
    * FeatureGenerator - generating features from repository. to see all options for json file please reffer to ::FeatureGeneratorTypes
    * FeatureProcessor - processing features from already generated/processed features. to see all options for json file please reffer to ::FeatureProcessorTypes
 * MedSamples - an object that stores our samples - with patient id, prediction time, label and more..
 * MedFeatures - an object that stores our matrix for all samples in MedSamples

  @subsection MedStats MedStats
 * MedBootstrapResult - the main object to run bootstrap on our data adn store the results, it contains the MedBootstrap as configuration for the bootstrap
 * bootstrap.h - the internal bootstrap infrastracture for more general use

 @subsection MedTime MedTime
 * MedTime.h - A library to convert between time units

  @subsection MedFeat MedFeat
 * MedOutcome - example object, need more documentation

 @subsection MedUtils MedUtils
 * MedMat - a more basic data structure to store matrix data
 * MedPlot.h -  A Library to plot graphs in HTML using plotly.js 
 
 @subsection Logger Logger
 * Several macros command for logging

 @section external Extrenal Libs

*  LightGBM
    * [Development Guide](md_External_LightGBM_LightGBM_docs_development.html)
	* [Quick Start](md_External_LightGBM_LightGBM_docs_Quick-Start.html)
	* [Main Documentation page](md_External_LightGBM_LightGBM_README.html)
	* [Parameter Tuning](md_External_LightGBM_LightGBM_docs_Parameters-tuning.html)
*  XGBoost
	* [Main Documentation page](md_External_xgboost_doc_index.html)
	* [Parameter Tuning](md_External_xgboost_doc_parameter.html)
	* [xgboost name space](namespacexgboost.html)
*  %VowpalWabbit - Yahoo Research library for machine learning (mainly for texts), also kernel svm's and more..
	* [vowpalwabitt name space](namespaceVW.html)
*  Eigen
	* EigenSolver

 @section comments General Comments
  * [TODO List](todo.html)
  * This mainpage can be found in $MR_ROOT/Libs/InternalOrig/medial_utilities/medial_utilities/mainpage.md