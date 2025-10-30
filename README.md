# Medial EarlySign Libraries

## Overview of Medial Infrastructure
This is an infrastructure developed by Medial EarlySign to streamline the creation of predictive models using EMR data for clinical applications. Existing tools often fall short for clinical use many Python libraries are not optimized for sparse time series analysis, leading to high memory consumption and, in some cases, performance that is 10–100 times slower than necessary.

Medial Infrastructure is designed to turn the Electronic Medical Record (EMR)-a complex, semi-structured time-series dataset, into a machine-learning-ready resource. Unlike images or free text, EMR data can be stored in countless formats, and its "labels" (the outcomes or targets you want to predict) aren’t always obvious. We address this by standardizing both the storage and the processing of time-series signals.
We can think about this infrastructure as "TensorFlow" of medical data machine learning. 

### Main contributers from recent years:
- [Avi Shoshan](https://www.linkedin.com/in/avi-shoshan-a684933b/)
- [Yaron Kinar](https://www.linkedin.com/in/yaron-kinar-il/)
- [Alon Lanyado](https://www.linkedin.com/in/lanyado/)

### Challenges
- **Variety of Questions**: Risk prediction (e.g., cancer, CKD), compliance, diagnostics, treatment recommendations
- **Medical Data Complexity**: Temporal irregularity, high dimensionality (>100k categories), sparse signals, multiple data types
- **Retrospective Data Issues**: Noise, bias, spurious patterns, policy sensitivity

### Goals
- Avoid reinventing common methodologies each project. Sometimes complicated code/logic with debugging
- Maintain shareable, versioned, regulatory‑compliant pipelines
- Facilitate reproducible transfer from research to product
- Provide end-to-end support: data import → analysis → productization

### Platform Requirements
- **Performance**: Ultra-efficient in memory & time (>100x compare to native python pandas in some cases, mainly in preprocessing)
- **Extensibility**: Rich APIs, configurable pipelines, support new data types
- **Minimal Rewriting & Ease Of Usage**: JSON‑driven configs, unified codebase, python API to the C library
- **Comprehensive**: From "raw" data to model deployment
- **Reproducible & Versioned**: Track data, code, models, and parameters


## Documentation

Please refer to [MR_WIKI](https://medial-earlysign.github.io/MR_Wiki/) for full documentation
- [Installation](https://medial-earlysign.github.io/MR_Wiki/Installation/)
