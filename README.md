# Medial EarlySign  Libraries

## Overview of Medial Infrastructure
Medial Infrastructure is designed to turn the Electronic Medical Record (EMR)—a complex, semi-structured time-series dataset—into a machine-learning-ready resource. 
Unlike images or free text, EMR data can be stored in countless formats, and its "labels" (the outcomes or targets you want to predict) aren’t always obvious. 
We address this by standardizing both the storage and the processing of time-series signals.

The Links refere to MR_Wiki
## Goals
1. **MedRepository: a high-performance EMR time-series store**
    * Fast retrieval of any patient’s full record or a specific signal across all patients.
    * [Unified representation](Generic%20(Universal)%20Signal%20Vectors): each signal consists of zero or more time channels plus zero or more value channels, all tied to a patient ID.
      - Static example: "Birth year" → no time channels, one value channel.
      - Single-time example: "Hemoglobin" → one time channel (test date), one value channel (numeric result).
      - Interval example: "Hospitalization" → two time channels (admission and discharge dates).
    * **Hierarchical support for categorical medical ontologies** 
      - Enables seamless integration and translation between different systems when working with a frozen model or algorithm. 
      - Example: A query for ICD-10 codes starting with "J" (respiratory diseases) will also automatically map to corresponding categories in systems like Epic. When dictionary of mapping between ICD and Epic is added, no need to change the model. 
      - Ontology mappings are managed by [MedDictionary](InfraMed%20Library%20page/MedDictionary), which supports many-to-many hierarchical relationships across coding systems.
2. **Modular processing pipeline (sklearn-style)**
    * **[Rep Processors](Rep%20Processors%20Practical%20Guide/)**: Clean or derive “raw” virtual signals, while preventing leakage of future data
      - Example: Outlier cleaner that omits values only when abnormality is detected by future readings (e.g., a hemoglobin value on 2023-Feb-04 flagged only by a 2023-May-21 test remains until after May 21).
      - Example: Virtual BMI signal computed from weight/height, or imputed when only two of three inputs exist
    * **[Feature Generators](MedProcessTools%20Library/FeatureGenerator/)**: Convert cleaned signals into predictive features.
      - Examples:
        * "Last hemoglobin in past 365 days"
        * "Hemoglobin slope over three years"
        * "COPD diagnosis code during any emergency admission in last three years"
    * **[Feature Processors](Feature%20Generator%20Practical%20Guide/)**: Operate on the feature matrix—imputation, selection, PCA, etc. 
    * **[Predictors/Classifiers](MedAlgo%20Library/)**: LightGBM, XGBoost, or custom algorithms.
    * **[Post-processing](PostProcessors%20Practical%20Guide/)**: Score calibration, explainability layers, fairness adjustments, etc.
3. **JSON-driven pipeline configuration** - Define every processor, feature generator, and model step in a single JSON file. [Json Format](MedModel%20json%20format)
   Example json for training a model:
   
<details>
  <summary>Click to expend</summary>
  
   ```json
   {
	"model_json_version": "2",
	"serialize_learning_set": "0",
	"model_actions": [
		"json:full_rep_processors.json", // Import a json from current folder with other componenets - in this case, outlier cleaners, signal panel completers, etc.
    // Features
		{
			"action_type": "feat_generator",
			"fg_type": "age"
		},
		{
			"action_type": "feat_generator",
			"fg_type": "gender"
		},
		{
			"action_type": "feat_generator",
			"fg_type": "unified_smoking",
			"tags": "smoking",
			"smoking_features": "Current_Smoker, Ex_Smoker, Unknown_Smoker, Never_Smoker, Passive_Smoker, Smok_Days_Since_Quitting , Smok_Pack_Years_Max, Smok_Pack_Years_Last,Smoking_Years,Smoking_Intensity"
		},
		// Cancers in Dx
		{
			"action_type": "feat_generator",
			"fg_type": "basic",
			"type": "category_set",
			"window": [
				"win_from=0;win_to=10950"
			],
			"time_unit": "Days",
			"sets": [
				"ICD9_CODE:140-149,ICD9_CODE:150-159,ICD9_CODE:160-165,ICD9_CODE:170,ICD9_CODE:171,ICD9_CODE:172,ICD9_CODE:174,ICD9_CODE:175,ICD9_CODE:176,ICD9_CODE:179-189,ICD9_CODE:200-208,ICD9_CODE:209.0,ICD9_CODE:209.1,ICD9_CODE:209.2,ICD9_CODE:290.3,ICD9_CODE:230-234"
			],
			"signal": "ICD9_Diagnosis",
			"in_set_name": "Cancers"
		},
    // Statistical features - will take: last, average, min, max, etc. for each time window: 0-180, 0-365. 365-730, 0-1095 prior prediction day in days and for each signal: Hemoglobin, WBC...
    // In total will create: 8*4*4 = 128 features
		{
			"action_type": "feat_generator",
			"fg_type": "basic",
			"type": [
				"last",
				"last_delta",
				"avg",
				"max",
				"min",
				"std",
				"slope",
				"range_width"
			],
			"window": [
				"win_from=0;win_to=180",
				"win_from=0;win_to=365",
				"win_from=365;win_to=730",
				"win_from=0;win_to=1095"
			],
			"time_unit": "Days",
			"tags": "labs_and_measurements,need_imputer,need_norm",
			"signal": [
				"Hemoglobin",
				"WBC",
				"Platelets",
				"Albumin"
			]
		},
		{
			"action_type": "feat_generator",
			"fg_type": "basic",
			"type": [
				"last_time"
			],
			"window": [
				"win_from=0;win_to=180",
				"win_from=0;win_to=365",
				"win_from=365;win_to=730",
				"win_from=0;win_to=1095"
			],
			"time_unit": "Days",
			"tags": "labs_and_measurements,need_imputer,need_norm",
			//Take only panels - to remove repititions:
			"signal": [
				"BMI",
				"Creatinine",
				"WBC",
				"Cholesterol",
				"Glucose",
				"Hemoglobin",
				"Albumin"
			]
		},
		{
			"action_type": "feat_generator",
			"fg_type": "category_depend",
			"signal": "DIAGNOSIS",
			"window": [
				"win_from=0;win_to=10950;tags=numeric.win_0_10950",
				"win_from=0;win_to=365;tags=numeric.win_0_365"
			],
			"time_unit_win": "Days",
			"regex_filter": "ICD10_CODE:.*",
			"min_age": "40",
			"max_age": "90",
			"age_bin": "5",
			"min_code_cnt": "200",
			"fdr": "0.01",
			"lift_below": "0.7",
			"lift_above": "1.3",
			"stat_metric": "mcnemar",
			"max_depth": "50",
			"max_parents": "100",
			"use_fixed_lift": "1",
			"sort_by_chi": "1",
			"verbose": "1",
			"take_top": "50"
		},
		// Feature selector to remove features with 99.9% same value, there are other options, like lasso, by model importance, etc.
		{
			"action_type": "fp_set",
			"members": [
				{
					"fp_type": "remove_deg",
					"percentage": "0.999"
				}
			]
		},
		// Imputer - simple choise of choosing median value by stratifying to age, gender and smoking status - will commit for all features with "need_imputer" tag
		{
			"action_type": "fp_set",
			"members": [
				{
					"fp_type": "imputer",
					"strata": "Age,40,100,5:Gender,1,2,1:Current_Smoker,0,1,1:Ex_Smoker,0,1,1",
					"moment_type": "median",
					"tag": "need_imputer",
					"duplicate": "1"
				}
			]
		},
		// Normalizer - will commit for all features with "need_imputer" tag
		{
			"action_type": "fp_set",
			"members": [
				{
					"fp_type": "normalizer",
					"resolution_only": "0",
					"resolution": "5",
					"tag": "need_norm",
					"duplicate": "1"
				}
			]
		}
	],
	"predictor": "xgb",
	"predictor_params": "tree_method=auto;booster=gbtree;objective=binary:logistic;eta=0.050;alpha=0.000;lambda=0.010;gamma=0.010;max_depth=6;colsample_bytree=0.800;colsample_bylevel=1.000;min_child_weight=10;num_round=200;subsample=0.800" }
```

</details>


4. **Comprehensive evaluation toolkit**
    * [Bootstrap-based](/Medial%20Tools/bootstrap_app/) cohort analysis allows batch testing across thousands of user-defined subgroups (e.g., age 50–80, males only, prediction window of 365 days, COPD patients).
    * Automatically extracts AUC, ROC points at each 1% FPR increment, odds ratios, PPV/NPV, and applies incidence-rate adjustments or KPI weights
    * Includes explainability and fairness audits
5. **Unified API wrapper for production deployment**
    * Ready for productization out of the box, no need to reinvent integration or design a new interface each time. See [AlgoMarker](AlgoMarkers/)
    * Packages the entire end-to-end pipeline (raw time-series ingestion through inference) into a single, stable SDK.
    * Core infrastructure implemented in C++ for performance and portability, with a lightweight [Python wrapper](/Python) for seamless integration.
    * Although powered by C++, the team mainly uses and maintains workflows via the Python SDK, ensuring rapid development and minimal friction. Experienced user might use the C++ API more often, since the python interface is more limited. 


## Installations
Those are installation steps required for all tools and builds:
### 1. Install Compiler and Build Tools (Ubuntu)
To install the essential compiler and build tools, run:
```bash
sudo apt install binutils gcc g++ cmake make swig -y
```
### 2. Install OpenMP Support (Ubuntu)
To enable OpenMP (used for parallel processing), install the following package:
```bash
sudo apt install libgomp1 -y
```
### 3. Install Boost Libraries (Ubuntu)
To install the required Boost components, on Ubuntu 24.04 use:
```bash
sudo apt install libboost-system1.83-dev libboost-filesystem1.83-dev libboost-regex1.83-dev libboost-program-options1.83-dev -y
```
> [!NOTE] On Ubuntu 22.04, Boost version 1.74 is available and is also compatible.
You may also choose to [download and compile Boost manually](https://www.boost.org/users/download/) if you prefer a different version. This project has been tested with Boost versions 1.67 through 1.85, and should work with other versions as well.

### There are 4 components:
1. AlgoMarker - please go to Internal/AlgoMarker and run ```./full_build.sh``` to create the algomarker shared library to pack a model for production
   You might need to recompile Boost library with -fPIC and then edit the ```CMakeLists.txt``` "set(BOOST_ROOT "$ENV{HOME}/boost-pic-install")" to point to your Boost compiles home, please put the compiled libraries in "/libs" and the headers in /include
3. Build python wrapper for this library. please execute ```Internal/MedPyExport/generate_binding/make-simple.sh```, before please make sure you have python3-dev headers, library and numpy installed: ```sudo apt install python3-dev -y``` in ubuntu
4. Build tools and executables that uses this library, please go to MR_Tools for more info, how to compile tools with this library
5. AlgoMarker Wrapper in c++ to expose the library as REST service. There is also an option to write it with python FastAPI and use the AlgoMarker library, but will also added an option to expose a server with a faster c++ library.
To do so, please go to MR_Tools and it will described there too
