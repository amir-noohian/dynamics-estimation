# 🤖 Dynamic Parameter Estimation for WAM Barrett Robot

This repository provides tools for **dynamic parameter estimation** of the **WAM Barrett robot**, supporting:

- Arbitrary **Degrees of Freedom (DOF)**
- Estimation **with or without gravity compensation**
- Evaluation of estimated parameters for **model-based control** or **simulation**

---

## 📁 Repository Structure

```
.
├── processed_data/
│   ├── processed_data_<data_id>.mat         # Filtered joint data (positions, velocities, etc.)
│   ├── pi_estimated_nograv_<pi_id>.mat      # Estimated dynamic parameters (no gravity)
│
├── wam_data/
│   ├── state_<id>.json                      # Raw joint state data
│   └── gravity_<id>.json                    # Raw gravity vectors
│
├── dataprocessing_nogravity.py              # Python script to filter and save joint data
├── calculate_pi_wam_nogravity.m             # MATLAB script to estimate dynamic parameters
├── estimate_or_evaluate_nograv.m            # MATLAB script for calculation/evaluation
└── README.md
```

---

## 🚀 Features

- ⚙️ Filter raw joint state data (position, velocity, torque, gravity)
- 🔍 Build regressor matrices based on the robot's dynamic model
- 📐 Estimate dynamic parameters (via least-squares)
- 🧪 Evaluate estimation accuracy using NRMSE and torque prediction
- 📊 Visualize recorded vs estimated torques

---

## 🔧 Setup

### Requirements

#### Python

Make sure you have the following Python libraries installed:

- `numpy`
- `scipy`
- `matplotlib`

Install them with:

```bash
pip install numpy scipy matplotlib
```

#### MATLAB

No additional toolboxes are required. All scripts use standard matrix operations and plotting.

---

## 🧪 Example Workflow

### 1. ✅ Preprocess WAM Data (Python)

```bash
python dataprocessing_nogravity.py 4
```

This will load `wam_data/state_4.json` and `gravity_4.json`, filter the signals, and save:

```
processed_data/processed_data_4.mat
```

---

### 2. 🧮 Estimate Parameters in MATLAB

```matlab
calculate_pi_wam_nogravity(4)
```

This loads `processed_data_4.mat`, estimates dynamic parameters, and saves:

```
processed_data/pi_estimated_nograv_4.mat
```

---

### 3. 📊 Evaluate or Calculate Parameters (Flexible)

Use the `estimate_or_evaluate_nograv.m` script to either **calculate** or **evaluate** parameter estimates:

#### ➕ Calculate new parameters using data ID 4 and save with PI ID 7:

```matlab
estimate_or_evaluate_nograv(4, 7, 'calc')
```

#### 📊 Evaluate an existing estimate (e.g., pi_estimated_nograv_7.mat) on data 4:

```matlab
estimate_or_evaluate_nograv(4, 7, 'eval')
```

This supports using different datasets for training and testing.

---

## 📌 Notes

- Gravity compensation is explicitly handled and removed during processing.
- Estimation uses SVD-based least-squares for numerical stability.
- The regressor is designed for a 2-DOF WAM robot but can be extended to higher DOFs.
- Outputs include:
  - Estimated dynamic parameters
  - Reconstructed torque signals
  - NRMSE error metrics
  - Plots comparing measured and estimated torques

---

## 👨‍💻 Author

**Amirhossein Noohian**  
PhD Student – Robotics, Dynamics, and Control  
University of Alberta
