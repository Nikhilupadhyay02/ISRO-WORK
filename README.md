"Development and Implementation of Extended Kalman Filter for Satellite Data and Lunar Rover Navigation"
This project presents a complete MATLAB-based implementation of the Extended Kalman Filter (EKF) tailored for satellite signal estimation and lunar rover navigation. The system was developed as part of a 3000+ line research-grade codebase authored during work with ISRO. It integrates simulation, real-world signal tracking, and dynamic visualization of Kalman filter behavior in harsh, noisy space environments.

📌 Project Overview
🔭 Designed for satellite-based signal tracking (e.g., L1 band)

🚜 Implemented for lunar rover state estimation using noisy sensor data

🧮 Uses an Extended Kalman Filter (EKF) to handle nonlinear system dynamics

📈 Real-time filtering and visualization using MATLAB .fig outputs

🧠 Key Features
Full state estimation: position, velocity, and attitude tracking

Extended Kalman Filter logic with nonlinear system modeling

Integration with simulated and actual space telemetry/sensor data

Robust against high noise, dropout, and sensor drift

Modular design for adaptation to planetary rovers or satellite platforms

📂 Included Files
File Name	Description
kalman_filter_1D.m	Baseline 1D Kalman filter implementation for simulation & testing.
Nikhil_kalman_static (1).fig	Static EKF result showing convergence in a constant-velocity system.
Nikhil_kalman_dynamic_plot (1).fig	Time-series plot of dynamic Kalman estimation vs. measurements.
Nikhil(plot)-SAT.L1 (2).fig	EKF tracking applied to satellite L1 signal — showcasing real-time signal filtering.


## 👨‍💻 Author & Code Ownership

This project, including all MATLAB scripts and visualizations, was entirely authored by **Nikhil Upadhyay** as part of a research and development initiative associated with ISRO.

- ✅ Over **3000+ lines of MATLAB code** written and tested
- ✅ Original implementation of **Extended Kalman Filter** for nonlinear dynamic systems
- ✅ Customized for **satellite data estimation** and **lunar rover navigation**
- ✅ All `.fig` output plots were generated from this authored codebase
