# 1-bit non-iterative photoacoustic computed tomography

This repository contains the MATLAB code and data to reproduce the main results presented in the manuscript "Reaching the noise-limited imaging depth and field of view via single-frame 1-bit non-iterative photoacoustic computed tomography".

---

## 1. System Requirements

### Software Dependencies & Operating Systems
* **Operating System:** The code has been tested on Windows 11. It is expected to be compatible with other operating systems (Linux, macOS).
* **MATLAB Version:** MATLAB R2022a or later is recommended.
* **Required MATLAB Toolboxes:**
    * **k-Wave Toolbox:** This toolbox is **required** to run the code for Figure 2. It can be downloaded from [http://www.k-wave.org/download.php](http://www.k-wave.org/download.php).

### Hardware Requirements
* **CPU:** Standard Intel or AMD processor.
* **RAM:** A minimum of 16 GB of RAM is recommended.
* **GPU:** A CUDA-enabled NVIDIA GPU is required to run the proprietary reconstruction function for Figures 3 and 4.

---

## 2. Installation Guide

1.  Clone or download this repository to your local machine.
2.  Install the **k-Wave Toolbox** by following the instructions on their official website.
3.  Open MATLAB and add this project's folder and all its subfolders to the MATLAB path.

* **Typical Install Time:** Setting up the project path should take less than 1 minute. The installation time for the k-Wave toolbox will vary.

---

## 3. Demo: Reproducing Figure 1b

Due to the stochastic nature of the added random noise, running the 1D 1-bit simulation code will not generate an identical result to the figure in the manuscript. To ensure direct reproducibility, the simulation data used to generate Fig. 1b has been provided.

* **Instructions:**
    1.  Navigate to the project's root directory in MATLAB.
    2.  The script `Fig1b_reproduce.m` will automatically load the required data from the `fig 1b simu data` folder.
    3.  Run the `Fig1b_reproduce.m` script.
* **Expected Output:** The script will generate a figure that is identical to Figure 1b in the manuscript.
* **Expected Run Time:** Approximately 10 seconds on a standard desktop computer.

---

## 4. Reproduction Instructions

### Figure 2
* Run the script `One_bit_DAS_fig2.m`.
* **Note:** This requires the k-Wave Toolbox to be installed and correctly added to the MATLAB path.

### Figures 3 & 4
* The script `One_bit_UBP_fig3_4.m` contains the primary code for data processing, including the binarization and bandpass filtering steps described in the manuscript.
* **Important Limitation:** The final UBP reconstruction cannot be fully executed with the provided code. This is because a key sub-function, `cudaTmDomnRecon3D20240618`, is **not included** in this repository.
* **Reason:** This function is proprietary and used in licensed technologies.
* **Availability:** This function is available from the corresponding author upon reasonable request.

### Videos
* To reproduce the results shown in the accompanying videos, please run the following scripts:
    * `One_bit_1D_videos1_3.m`
    * `Direct_avergae_compare_1bit_video4.m`

---

## License

The code provided in this repository is licensed under the MIT License. See the `LICENSE` file for details.

Please note that this license does not apply to the proprietary function `cudaTmDomnRecon3D20240618`, which is subject to a separate agreement.

---
