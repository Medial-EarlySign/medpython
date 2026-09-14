# Medial Python Binding

## Relase Notes - 1.2.1
* Feature: Update XGBoost to 3.4.1
* Fix Support for SWIG 4.5
* Close all CodeQL potential vulnerabilities

## Relase Notes - 1.2.0
* Bugfix: Compile in low memory system
* Feature: Update LightGBM to 4.6.0
* Feature: Support compilation with newer SWIG >= 4.3
* Feature: Removed all external libraries from source code and using FetchContent or using system library if exists. Important to support conda
* Feature: provide pre-compiled binaries for python 3.14 to windows and mac instead python 3.10

## Relase Notes - 1.1.2
* Feature: Retreive string attributes from samples

## Relase Notes - 1.1.1
* Bugfix: get_sig when using PidRepository data from json
* Feature: add support to retrieve model arch as json
* Feature: filter MedSamples by bootstrap query language

## Relase Notes - 1.1.0
* **Modularized Package Structure**
    - **Decoupled** `medpython-etl` into a standalone package.
    - `medpython` now acts as a **meta-package**, including `medpython-etl` as a core dependency for user convenience without introducing direct code coupling.
* Enhanced `med.Model` **Functionality**
    - **Runtime Model Updates**: Introduced apply_model_changes to modify models on the fly (e.g., removing explainability or lowering batch counts to save memory) without requiring retraining.
    - **Post-Processor Integration**: Added add_post_processors_json_string_to_model to attach and train new components such as calibrators, explainers, or fairness adjustments onto existing/trained models.
    - **Component Management**: New train/remove/print commands for rep-processors and post-processors, allowing direct execution of training functions for specific model components.
    - **Predictor Export**: Added functionality to retrieve the core predictor in standard scikit-learn or XGBoost formats (for supported predictors)
    - **Architecture Visualization**: Added the ability to retrieve the MedModel architecture in a human-readable, printable format.
    - **Usage Validation**: Implemented strict state checks to prevent logical errors (e.g., ensuring apply cannot be called before the model is trained).
    - **Signal Debugging**: Introduced `debug_rep_processor_signal` to inspect how rep-processors transform raw signals, outlier, etc.
* **Platform & Build Improvements**
    - **macOS Compatibility**: Resolved all compilation issues for **macOS** environments using Conda, ensuring a smoother installation process for Apple Silicon and Intel-based Macs.

## Relase Notes - 1.0.6
* Fix code to compile in ARM computer
* Removed ipython dependency
* Readme fixes
* Removed the Boost libraries dependency, using only headers. Simpler compile process

## Relase Notes - 1.0.5
* Adjustments to compile with Alpine/Musl library
* Bugfix issue #5
* Fix code to support and compile in osx

## Relase Notes - 1.0.4
* Fixed all issues to compile also in windows

## Relase Notes - 1.0.3
* Updated documentation

## Release Notes - 1.0.2
* First version deployed on PyPi with CI/CD of github
* Added MedConvert to med library to allow loading of data into PidRepository data format without external tools like `Flow`
* Added compilation ability based on shared library of Boost for self compliation ability

## Release Notes - 1.0.1
* First version - depolying current status
