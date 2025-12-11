% =========================================================================
% batch_process_single_patient_ACAI_preprocessed.m
%
% This script processes a single patient for ACAI analysis.
% !!!! IT ASSUMES CAT12 SEGMENTATION HAS ALREADY BEEN RUN.
% !!!! IT DOES NOT OUTPUT CLUSTERS BUT THE ACAI IMAGE IN PATIENT AND MNI
% SPACE
% !!!! IT ASSUMES MNI VOXEL DIMENSIONS ARE 1 x 1 x 1
% It finds existing CAT12 outputs (p1, y, iy) based on the T1 filename
% and uses them with the coregistered PET to calculate the ACAI map.

% =========================================================================
clear; clc;

% --- 0. CONFIGURATION ---
% Full path to the original T1 NIfTI file.
% This is used to FIND the CAT12 outputs (e.g., mri/p1T1.nii)
t1_image_file = '';
% Full path to the FDG-PET NIfTI file (must be coregistered to the T1)
pet_image_file = '';

% Paths to required toolboxes and scripts
% for example C:\MATLAB\R2023b\toolbox\spm12
spm_path = '';
% Path to the directory containing 'LCN12_calc_ACAI_GACAI_preprocessed.m'
% ...\acai_single_patient_simple
acai_script_path = '';

% ACAI Kernel settings
kernel_type = 'Gaussian'; % 'Gaussian', 'Cubic', 'Sphere'
kernel_size = 9;          % FWHM for Gaussian, side for Cubic, diameter for Sphere
% -----------------------------

% --- 1. Setup Environment ---
fprintf('--- 1. Setting up Environment ---\n');
try
    addpath(spm_path);
    addpath(acai_script_path);
    
    spm_jobman('initcfg');
    spm('defaults', 'PET');
    fprintf('  SPM and script paths added.\n');
catch ME
    fprintf('  ERROR: Could not initialize SPM or add paths.\n');
    fprintf('  Please check your `spm_path` and `acai_script_path` variables.\n');
    rethrow(ME);
end

% --- 2. Define CAT12 Output File Paths ---
fprintf('\n--- 2. Defining CAT12 Output Paths (assuming segmentation is complete) ---\n');
% CAT12 places outputs in an 'mri' subdirectory in the T1's folder
[t1_path, t1_name, t1_ext] = spm_fileparts(t1_image_file);
cat12_mri_path = fullfile(t1_path, 'mri');

% Define paths for the files needed by ACAI
gm_image_native = fullfile(cat12_mri_path, ['p1' t1_name t1_ext]);
deformationfield_forward = fullfile(cat12_mri_path, ['y_' t1_name t1_ext]);
deformationfield_inverse = fullfile(cat12_mri_path, ['iy_' t1_name t1_ext]);
target_GMactivityimage_native = ''; 

fprintf('  Base T1: %s\n', t1_image_file);
fprintf('  Expecting CAT12 ''mri'' folder: %s\n', cat12_mri_path);

% --- 3. Verify All Required Files Exist ---
fprintf('\n--- 3. Verifying All Files for ACAI ---\n');
files_ok = true;
if ~exist(pet_image_file, 'file')
    fprintf('  ERROR: Coregistered PET file not found: %s\n', pet_image_file);
    files_ok = false;
end
if ~exist(gm_image_native, 'file')
    fprintf('  ERROR: CAT12 GM (p1) file not found: %s\n', gm_image_native);
    files_ok = false;
end
if ~exist(deformationfield_forward, 'file')
    fprintf('  ERROR: CAT12 Forward deformation (y_) not found: %s\n', deformationfield_forward);
    files_ok = false;
end
if ~exist(deformationfield_inverse, 'file')
    fprintf('  ERROR: CAT12 Inverse deformation (iy_) not found: %s\n', deformationfield_inverse);
    files_ok = false;
end

if ~files_ok
    fprintf('  --> STOPPING processing due to missing files.\n');
    return;
end

fprintf('  Found All Required Files:\n');
fprintf('    - PET: %s\n', pet_image_file);
fprintf('    - GM (p1): %s\n', gm_image_native);
fprintf('    - Def. Fwd (y_): %s\n', deformationfield_forward);
fprintf('    - Def. Inv (iy_): %s\n', deformationfield_inverse);

% --- 4. Call the ACAI calculation function ---
fprintf('\n--- 4. Calling LCN12_calc_ACAI_GACAI_preprocessed ---\n');
try
    LCN12_calc_ACAI_GACAI_preprocessed(pet_image_file, gm_image_native, ...
        deformationfield_forward, deformationfield_inverse, kernel_type, kernel_size, ...
        target_GMactivityimage_native);
        
    fprintf('  --> SUCCESS: ACAI processing complete.\n');
catch ME
    fprintf('  --> ERROR during ACAI calculation: %s\n', ME.message);
    for k=1:length(ME.stack)
        fprintf('      In %s at line %d\n', ME.stack(k).name, ME.stack(k).line);
    end
end

fprintf('\n--------------------------------------------------\n');
fprintf('Batch processing complete for this patient.\n');