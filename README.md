## 🌊 Project Motivation

Lake Mälaren is Sweden’s third-largest lake and provides drinking water for over **2 million people**. Understanding **how water flows between different basins** is crucial, as some basins receive polluted inflows while others serve as drinking water reservoirs.

While a **physics-based hydrodynamic model** exists, it is:

* Computationally expensive to run.
* Difficult to interpret in terms of specific driving factors.

This project utilizes **statistical and machine-learning models** to:

* Predict flow across selected basin interfaces.
* Identify how **river inflow** and **meteorological conditions** influence transport.
* Enable faster, real-time-capable, and more interpretable analysis.

---

## 🎯 Objectives

* **Predict discharge** across selected interfaces (**A, B, and C**) between lake basins.
* **Feature Engineering:** Use river inflow and meteorological data as primary input features.
* **Model Comparison:** Compare interpretable models (Linear Regression, Decision Trees) against complex models (Neural Networks, XGBoost).
* **Reproducibility:** Preserve time-series structure and ensure a fully automated pipeline.

---

## 📂 Repository Structure

```text
.
├── Meteo_csv/              # Meteorological input data (CSV)
├── Rivers_csv/             # River inflow input data (CSV)
├── Rivers_csv_old/         # Older / unused river CSV files
├── Python_csv/             # Intermediate CSVs used by Python models
│
├── model.py                # Main command-line entry point
├── _model.py               # Training, evaluation, and prediction logic
├── utils.py                # Data loading and preprocessing
├── const.py                # Global constants (paths, filenames)
├── requirements.txt        # Python dependencies
│
├── decision_tree_A.py      # Decision tree model (interface A)
├── neural_network_A.py     # Neural network model (interface A)
│
├── linear_regression_A.R   # Linear regression model (interface A)
├── linear_regression_B.R   # Linear regression model (interface B)
├── linear_regression_C.R   # Linear regression model (interface C)
├── xgboost_C.R             # XGBoost model (interface C)
│
├── time_series_model.R     # Time-series baseline models
├── meteo_average.R         # Meteorological data aggregation
├── river_average.R         # River inflow aggregation
├── flow_at_interface.R     # Flow extraction at basin interfaces
├── flow_at_coordinate.R    # Flow extraction at spatial coordinates
├── flux_interface.R        # Flux calculations across interfaces
│
└── README.md

```

---

## 📊 Data Description

### Hydrodynamic Target Data

* Generated from a physics-based hydrodynamic model.
* Weekly averages from **2000–2024**.
* Flow extracted across **9 interfaces** grouped into **A, B, and C**.

### Meteorological Data (`Meteo_csv/`)

* Variables: Wind speed, wind direction, precipitation, temperature.
* Aggregated to **weekly resolution**.
* A `month` feature is extracted to account for seasonality.

### River Inflow Data (`Rivers_csv/`)

* Daily inflow data aggregated to weekly values.
* Used to capture basin-specific hydrological forcing.

---

## 🧠 Modeling Approach

We evaluated several approaches to find the optimal balance between accuracy and transparency:

* **Linear Regression** (Final selected model)
* **Decision Trees** (Interpretable comparison)
* **Neural Networks** (Performance benchmark)
* **XGBoost** (R-based gradient boosting comparison)

**Key Design Choices:**

* **Ordered 80/20 train–test split:** Preserves temporal order (no shuffling).
* **Multi-output regression:** Handles multiple interfaces simultaneously.
* **Lagged features:** Incorporates delayed effects of precipitation.

---

## 🚀 Python Pipeline Usage

### 1. Installation

```bash
pip install -r requirements.txt

```

### 2. Training and Evaluation

Run the main pipeline to train the model on available data and generate evaluation metrics:

```bash
python model.py -d Rivers_csv -m Meteo_csv -s 0.8 -i A,B,C

```

* Trains a linear regression model.
* Evaluates on held-out test data.
* Saves models, predictions, and figures automatically.

### 3. Prediction Only

Use a previously trained model to generate new predictions:

```bash
python model.py -m Meteo_csv -r Rivers_csv -s 0

```

---

## 📈 Outputs

* **Models:** Saved in `models/linear_regression/` as `.joblib` files.
* **Figures:** Evaluation and prediction plots saved in `figures/`.
* **Predictions:** Raw output data saved in `predictions/linear_regression_predictions.csv`.

---

## ✨ Results Summary

* **Linear Regression** achieved the best balance between performance and interpretability.
* **River inflow** was identified as the strongest predictor of inter-basin exchange.
* **Meteorological variables** provided significant secondary explanatory power.
* More complex models did not significantly outperform linear regression, likely due to the dataset size and the underlying physics of the system.

---

## 🎓 Course Context

This repository accompanies the project report for the **Project in Data Science** course at **Uppsala University**.

**Authors:** Alberto Fernández Comesaña,  Henrik Jonasson, Moutushi Sen

**Project:** *Predicting Inter-Basin Exchange in Lake Mälaren*

---

## 📜 License

For academic and educational use only.

---

