# ACAI/GACAI Single Patient Processing Pipeline

## Overview

This repository provides tools to calculate **ACAI** (Asymmetry Corrected for Anatomy Index) and **GACAI** (Global ACAI) metabolic indices for single subjects. These indices are designed to detect focal metabolic asymmetries in FDG-PET imaging, corrected for normal anatomical asymmetry.

The methods are described in the paper:  
> **"Standardised anatomy-corrected FDG-PET asymmetry predicts seizure freedom in focal cortical dysplasia"** *(under review)*

---

## Available Pipelines

We provide **three different pipeline setups** to accommodate various use cases:

| Pipeline | Folder | Description |
|----------|--------|-------------|
| **cat12_ACAI** | `cat12_ACAI/` | Main pipeline using CAT12 for segmentation. **Used in the paper.** |
| **ACAI_with_prior_seg** | `ACAI_with_prior_seg/` | For cases where CAT12 segmentation has already been performed (e.g., with custom/strict parameters). |
| **spm12_ACAI** | `spm12_ACAI/` | Original implementation by Prof. Dupont using SPM12 segmentation. |

---

## 1. cat12_ACAI (Recommended)

This is the **primary pipeline** used in the paper. It performs CAT12 segmentation as part of the processing workflow.

### Description
`process_single_patient.m` is a MATLAB automation script that:
- Initializes SPM12/CAT12
- Runs CAT12 segmentation on the T1 image
- Calculates ACAI and GACAI indices

It acts as a wrapper around the core logic script (`LCN12_calc_ACAI_GACAI_preprocessed.m`).

---

## 2. ACAI_with_prior_seg

This pipeline is designed for cases where **CAT12 segmentation has already been performed** separately.

### When to Use
- When you need to use **custom or stricter segmentation parameters**
- When you have already processed subjects through CAT12 for other analyses
- When you want more control over the segmentation step
- When batch processing where segmentation was done in a separate step

### Workflow
1. Run CAT12 segmentation separately with your preferred parameters
2. Use this pipeline to calculate ACAI/GACAI from the existing segmentation outputs

---

## 3. spm12_ACAI

This is the **original implementation** developed by Prof. Dupont, using SPM12's unified segmentation instead of CAT12.

### When to Use
- When CAT12 is not available
- For comparison with the original method

### Differences from cat12_ACAI
- Uses SPM12's native tissue segmentation (New Segment)
- May produce slightly different results due to segmentation algorithm differences

---

## Prerequisites

### MATLAB Requirements
- **MATLAB** R2023b

### Required Toolboxes
- **SPM12**: Statistical Parametric Mapping toolbox
- **CAT12**: Computational Anatomy Toolbox (installed within `spm12/toolbox/cat12`)
  - *Note: Only required for `cat12_ACAI` and `ACAI_with_prior_seg` pipelines*
  - v12.9 was used in our paper

### Core Scripts
- `LCN12_calc_ACAI_GACAI_preprocessed.m` must be available locally

### Python Dependencies (for post-processing)
- **Python 3.12**
- Required libraries:
  - `nibabel`
  - `nilearn`
  - `numpy`
  - `pandas`
  - `scipy`
  - `openpyxl` (for Excel export)

Install Python dependencies via pip:
```bash
pip install nibabel nilearn numpy pandas scipy openpyxl
```

---

## Input Data Requirements

* **T1 Image**: High-resolution structural MRI (NIfTI format)
* **PET Image**: FDG-PET image (NIfTI format)
  * **Important:** The PET image must be **coregistered** to the T1 image before running this pipeline

---

## Configuration

Open `process_single_patient.m` (in the appropriate pipeline folder) and edit the **CONFIG** section at the top of the file.

### 1. File Paths
* `t1_image_file`: Full absolute path to the subject's **original T1** NIfTI file
* `pet_image_file`: Full absolute path to the subject's **coregistered FDG-PET** NIfTI file

### 2. Environment Paths
* `spm_path`: Root directory of your SPM12 installation
  * *Example:* `C:\MATLAB\R2023b\toolbox\spm12`
  * *Note:* This folder must contain the `tpm` subdirectory and `toolbox/cat12`
* `acai_script_path`: Directory containing the helper script `LCN12_calc_ACAI_GACAI_preprocessed.m`

### 3. Algorithm Settings
* `kernel_type`: The smoothing kernel type
  * Options: `'Gaussian'`, `'Cubic'`, `'Sphere'`
  * Normal option: Gaussian
* `kernel_size`: The dimensions of the kernel
  * Gaussian: **9 mm FWHM is the standard, used in the paper**
  * Cubic: Side length
  * Sphere: Diameter

---

## Usage

### For cat12_ACAI (standard pipeline)

1. Navigate to the `cat12_ACAI/` folder
2. Open `process_single_patient.m` in MATLAB
3. Modify the variables in the `0. CONFIGURATION` section
4. Run the script
5. Use the Jupyter notebook to extract clusters based on ACAI and volume thresholds

### For ACAI_with_prior_seg

1. First, run CAT12 segmentation separately with your desired parameters
2. Navigate to the `ACAI_with_prior_seg/` folder
3. Open the processing script and configure paths to your existing segmentation outputs
4. Run the script
5. Use the Jupyter notebook to extract clusters based on ACAI and volume thresholds

### For spm12_ACAI

1. Navigate to the `spm12_ACAI/` folder
2. Open the processing script in MATLAB
3. Configure as needed and run
4. Use the Jupyter notebook to extract clusters based on ACAI and volume thresholds. Set the atlas filepath as None.

---

## Example Configuration

```matlab
% --- 0. CONFIGURATION ---
t1_image_file    = 'C:\Data\Patient001\T1.nii';
pet_image_file   = 'C:\Data\Patient001\PET_coreg.nii';
spm_path         = 'C:\MATLAB\R2023b\toolbox\spm12'; 
acai_script_path = 'C:\Scripts\ACAI_Toolbox\';
kernel_type      = 'Gaussian';
kernel_size      = 9;          
% -----------------------------
```

---

## Output

The pipeline generates ACAI and GACAI maps and values that can be used for:
- Detecting focal metabolic asymmetries
- Predicting seizure freedom in focal cortical dysplasia
- Localizing epileptogenic zones

---

## Citation

If you use this code, please cite:

> "Standardised anatomy-corrected FDG-PET asymmetry predicts seizure freedom in focal cortical dysplasia" (under review)

---

## Contact

For questions or issues, please contact the authors.
jeroen.gijs@uzleuven.be

## License

MIT license