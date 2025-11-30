# Medial EarlySign Python Library


![Pepy Total Downloads](https://img.shields.io/pepy/dt/medpython)
![PyPI - License](https://img.shields.io/pypi/l/medpython)
![GitHub contributors](https://img.shields.io/github/contributors-anon/Medial-EarlySign/MR_LIBS)
![GitHub commit activity](https://img.shields.io/github/commit-activity/t/Medial-EarlySign/MR_LIBS)

[![GitHub Repo](https://img.shields.io/badge/github-repo-blue?logo=github)](https://github.com/Medial-EarlySign/MR_LIBS)
![GitHub Repo stars](https://img.shields.io/github/stars/Medial-EarlySign/MR_WIKI)


Our platform is designed to transform complex, semi-structured Electronic Medical Records (EMR) into **machine-learning-ready** data and reproducible model pipelines. The framework is optimized for the unique challenges of sparse, time-series EMR data, delivering **low memory usage** and **high-speed processing** at scale.

It was conceived as a **TensorFlow** for machine learning on medical data.

All software is now open-sourced under the MIT license. Some of the models developed by Medial EarlySign that are currently in production are available exclusively through our partners.

The framework was battle-tested in production across multiple healthcare sites and was a key component of an **award-winning** submission to the [CMS AI Health Outcomes Challenge](https://www.cms.gov/priorities/innovation/innovation-models/artificial-intelligence-health-outcomes-challenge).

## Why Use This Platform?

*   **High-Performance Processing:** Engineered for large-scale, sparse EMR time-series data where general-purpose libraries like pandas fall short.
*   **Reusable Pipelines:** Save valuable engineering time by providing shareable, tested pipelines and methods.
*   **Built-in Safeguards:** Mitigate common pitfalls like data leakage and time-series-specific overfitting.
*   **Production-Ready:** Designed for easy deployment using Docker or minimal distroless Linux images.
*   **State of The Art Algorithms** Developed state of the art algorithms for EMR use cases, explainablity and fairness and more...

## Core Components

The platform is built on three key pillars:

*   **MedRepository:** A compact, efficient data repository and API for storing and accessing EMR signals. Querying categorical signals like perscriptions and diagnosis in an easy and efficient API. 
*   **MedModel:** An end-to-end machine learning pipeline that takes data from MedRepository or JSON EMR inputs to produce predictions and explainability outputs. It supports both training and inference.
*   **Medial Tools:** A suite of utilities for training, evaluation, and workflow management, including bootstrap analysis, fairness checks, and explainability.

## Setup

You can quickly install the package using **pip**:

```bash
pip install medpython
```

**System Requirements**

* **Supported Systems**: This pre-built version is available for **modern Linux** distributions (specifically `manylinux2014` equivalents, such as CentOS >= 7 or Ubuntu >= 13.04). The software also compiles in **Windows** but you will need to install Boost yourself. I hope shortl to provide windows pre-compiled builds through conda.
* **Python**: Requires **Python 3.10 through 3.14**

**Compilation for Other Systems**
If you're using an **older Linux** or a **different platform/Python version >= 3.8**, you'll need to **compile the package yourself**.

* **Note on Compilation**: Ensure the **Boost libraries** are installed. For a local setup, set the environment variable `BOOST_DISABLE_STATIC=1` to link against shared Boost libraries (The reason is that your system static libraries weren't compiled with `-fPIC` flag, so you can't use them inside python module):
```bash
export BOOST_DISABLE_STATIC=1
```
You can also set Boost installation directory with `BOOST_ROOT` environment variable if it is not part of the system libraries.
Currently the library was compiled/tested in **linux** and **windows** only.

## Usage
```python
import med
from AlgoMarker import AlgoMarker
from ETL_Infra import prepare_final_signals, prepare_dicts, finish_prepare_load, create_train_signal
```

More information on usage:

* [ETL_Infra](https://medial-earlysign.github.io/MR_Wiki/Tutorials/01.ETL%20Tutorial) - To load data into medpython data format
    - [API Documentation](https://medial-earlysign.github.io/MR_Tools)
* [AlgoMarker](https://medial-earlysign.github.io/MR_Wiki/Infrastructure%20Library/Medial%20Tools/Python/Python%20AlgoMarker%20API%20Server.md#simple-usage-example)
    - From usage in deployment/lighter API interface. Usefull for accessing prebuild models. [API Documentation](https://medial-earlysign.github.io/MR_Tools/AlgoMarker.html)
* [med library](https://medial-earlysign.github.io/MR_Wiki/Infrastructure%20Library/Medial%20Tools/Python/Examples.md)
    - Full library from training new models or using with existing models. [Complete Example: From Data to Model](https://medial-earlysign.github.io/MR_Wiki/#complete-example-from-data-to-model)
    - [API Documentation](https://medial-earlysign.github.io/MR_Tools)

## Getting Started

*   **Build a new model:** Follow the step-by-step [Tutorials](https://medial-earlysign.github.io/MR_Wiki/Tutorials/index.html) to build a model from scratch.
*   **Use an existing model:** Browse the collection of [Models](https://medial-earlysign.github.io/MR_Wiki/Models/index.html).