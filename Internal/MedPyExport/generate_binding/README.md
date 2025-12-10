# Medial EarlySign Python Library


![Pepy Total Downloads](https://img.shields.io/pepy/dt/medpython)
![PyPI - License](https://img.shields.io/pypi/l/medpython)
![GitHub contributors](https://img.shields.io/github/contributors-anon/Medial-EarlySign/medpython)
![GitHub commit activity](https://img.shields.io/github/commit-activity/t/Medial-EarlySign/medpython)

[![GitHub Repo](https://img.shields.io/badge/github-repo-blue?logo=github)](https://github.com/Medial-EarlySign/medpython)
![GitHub Repo stars](https://img.shields.io/github/stars/Medial-EarlySign/medpython)

**A note on our journey:** Medial EarlySign was a company that developed a proprietary platform for machine learning on electronic medical records. Following the company's liquidation, the decision was made to release the core software as an **open-source** project to allow the community to benefit from and build upon this technology. Please feel free to [reach me out](https://medial-earlysign.github.io/MR_Wiki/#community-and-contributions) in any case of an issue. I'm voluntarily holding this, so please be patients.

Our platform is designed to transform complex, semi-structured Electronic Medical Records (EMR) into **machine-learning-ready** data and reproducible model pipelines. The framework is optimized for the unique challenges of sparse, time-series EMR data, delivering **low memory usage** and **high-speed processing** at scale.

It was conceived as a **TensorFlow** for machine learning on medical data.

All software is now open-sourced under the MIT license. Some of the [models](http://medial-earlysign.github.io/MR_Wiki/Models) developed by Medial EarlySign that are currently in production are available exclusively through our partners.

The framework was battle-tested in production across multiple healthcare sites and was a key component of an **award-winning** submission to the [CMS AI Health Outcomes Challenge](https://www.cms.gov/priorities/innovation/innovation-models/artificial-intelligence-health-outcomes-challenge).

## Why Use This Platform?

*   **High-Performance Processing:** Engineered for large-scale, sparse EMR time-series data where general-purpose libraries like pandas fall short.
*   **Reusable Pipelines:** Save valuable engineering time by providing shareable, tested pipelines and methods.
*   **Built-in Safeguards:** Mitigate common pitfalls like data leakage and time-series-specific overfitting.
*   **Production-Ready:** Designed for easy deployment using Docker or minimal distroless Linux images.
*   **State of The Art Algorithms** Developed state of the art algorithms for EMR use cases, explainablity and fairness and more.
    - Read My Journey on [Explainability](https://medial-earlysign.github.io/MR_Wiki/Blog/explainability.html)

## Core Components

The platform is built on three key pillars:

*   **MedRepository:** A compact, efficient data repository and API for storing and accessing EMR signals. Querying categorical signals like perscriptions and diagnosis in an easy and efficient API. 
*   **MedModel:** An end-to-end machine learning pipeline that takes data from MedRepository or JSON EMR inputs to produce predictions and explainability outputs. It supports both training and inference.
*   **Medial Tools:** A suite of utilities for training, evaluation, and workflow management, including bootstrap analysis, fairness checks, and explainability. `MedPython` is the Python alternative to replace all those tools.

## From Raw Data to Insight in Four Simple Steps

Our platform streamlines the development and deployment of clinical predictive models, transforming raw patient data into actionable insights. For live predictions (inference), you can use raw JSON data directly, bypassing the need for an optimized data store.

<img src="https://medial-earlysign.github.io/MR_Wiki/images/MES_Arch.png">

This structured approach ensures that data is processed efficiently, models are built systematically, and the results are both accurate and interpretable.

### The Workflow

**1. Start with Raw Patient Data**

Begin with your data in a simple JSON format.

```json
{
  "patient_id": "1",
  "data": {
    "signals": [
      {
        "code": "Hemoglobin",
        "unit": "g/dL",
        "data": [
          { "timestamp": [20240806], "value": ["14.1"] },
          { "timestamp": [20250806], "value": ["14.5"] }
        ]
      },
      {
        "code": "Diagnosis",
        "data": [
          { "timestamp": [20240701], "value": ["ICD10_CODE:J20"] },
          ...
        ]
      },
      ...
    ]
  }
}
```

Load it into our Optimized Store or use it "as is" in deployment.

**2. Define Your Label**

For each patient, and for any chosen prediction date (after which no future information is provided to the model), specify what the outcome label should be for training or testing

**3. Define Your ML Pipeline**

Configure your entire machine learning workflow from preprocessing and feature engineering to the final model using a single [JSON configuration file](https://medial-earlysign.github.io/MR_Wiki/Infrastructure%20Library/MedModel%20json%20format.html). This approach ensures your experiments are reproducible and easy to version.

**4. Get Explainable Predictions**

Train your model using the Python SDK and generate predictions with clear, interpretable explanations. [Example output](https://medial-earlysign.github.io/MR_Wiki/Tutorials/07.Deployment/index.html#how-to-use-the-deployed-algomarker).

This is an illustration of the final output in a visual format (Our infrastructure returns the data to create this):

<img src="https://medial-earlysign.github.io/MR_Wiki/images/Explainability.png">

For more details, refer to the [Tutorials](https://medial-earlysign.github.io/MR_Wiki/Tutorials).

## Setup

You can quickly install the package using **pip**:

```bash
pip install medpython
```

**System Requirements**

* **Supported Systems**: This pre-built version is available for **modern Linux** distributions (specifically `manylinux2014` equivalents, such as CentOS >= 7 or Ubuntu >= 13.04). The software also compiles in **Windows**, **Alpine** and **Mac** but you will need to install Boost yourself. I hope shortly to provide windows/macos pre-compiled builds.
* **Python**: Requires **Python 3.10 through 3.14**

**Compilation for Other Systems**
If you're using an **older Linux** or a **different platform/Python version >= 3.8**, or a version without a pre-build, you will need to **compile the package yourself**.
More information in [Installation](http://medial-earlysign.github.io/MR_Wiki/Installation/)

## Release Notes
Please refer to [Release Notes](https://github.com/Medial-EarlySign/medpython/discussions/categories/announcements)

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

## Community and Contributions

This is an open-source project, and we welcome contributions from the community.

*   **Report issues or ask questions:** Please use our [Github Discussions](https://github.com/Medial-EarlySign/medpython/discussions).
*   **Contribute to the code:** Check out our repositories:
    *   [MedPython](https://github.com/Medial-EarlySign/medpython): The core infrastructure libraries.
    *   [MR_Tools](https://github.com/Medial-EarlySign/MR_Tools): Tools and pipelines built on top of medpython.
    *   [MR_Scripts](https://github.com/Medial-EarlySign/MR_Scripts): Dev Ops tools