function LCN12_calc_ACAI_GACAI_preprocessed(pet_image_native, gm_image_native, deformationfield_forward, deformationfield_inverse, kernel_type, kernel_size, varargin)
% LCN12_calc_ACAI_GACAI_preprocessed(pet_image_native, gm_image_native, deformationfield_forward, deformationfield_inverse, kernel_type, kernel_size, target_GMactivityimage_native)
%
% Calculates Asymmetry Index (ACAI) based on CAT12 deformation fields
% The images are expected in native space and will be warped to MNI for calculations.
%
% Input:
%   pet_image_native:
%      Full path to the PET image in native subject space.
%   gm_image_native:
%      Full path to the grey matter probability map in native subject space
%      (e.g., p1T1.nii from CAT12 before warping). This image will be warped to MNI.
%   deformationfield_forward:
%      Full path to the forward deformation field (native to MNI, e.g., y_T1.nii).
%   deformationfield_inverse:
%      Full path to the inverse deformation field (MNI to native, e.g., iy_T1.nii).
%   kernel_type:
%      Kernel type: 'Gaussian', 'Cubic' or 'Sphere'. Default: 'Gaussian'.
%   kernel_size:
%      Size of the kernel (in voxels).
%      'Gaussian': FWHM; 'Cubic': cube side; 'Sphere': sphere diameter.
%   Optional: target_GMactivityimage_native:
%      Full path to a PET image of only GM contribution (e.g., from AMAP),
%      in native subject space. If provided, GACAI will also be calculated.
%
% The results will be written to file in the same directory as pet_image_native.
% A subfolder 'MNI_space' will be created for intermediate MNI-space files.
%
% Final outputs:
%   ACAI_native_[pet_basename]_[gm_suffix]_[kernelinfo].nii
%   GACAI_native_[pet_basename]_[gm_suffix]_[kernelinfo].nii (if applicable)
%
% Important:
%  - Requires SPM12 installation and LCN12 helper functions (LCN12_read_image, LCN12_write_image).
%  - Assumes all input images are correctly oriented.
%
% Based on LCN12_calc_ACAI_GACAI.m by Lin Zhou, Kathleen Vunckx, Patrick Dupont.
% Modified for pre-processed inputs.
%__________________________________________________________________________
% =========================================================================
% CREDIT & DISCLAIMER
% =========================================================================
% Original ACAI implementation by Patrick Dupont and the LCN team.
% Comments or questions can be sent to: Patrick.Dupont@kuleuven.be
%
% Important:
%  - requires the installation of CAT12 v12.9 and SPM12 - http://www.fil.ion.ucl.ac.uk/spm/
%  - the path of SPM and this routine should be included in the Matlab path
%
% The package contains software (SPM12) developed under the auspices of The
% Wellcome Department of Imaging Neuroscience, a department of the
% Institute of Neurology at University College London. The copyright of
% this software remains with that of SPM12, see
% http://www.fil.ion.ucl.ac.uk/spm/.
%    
% This routine is supplied as is. 
%
% IMPORTANT REMARKS: 
%   - this is research software.
%   - always check the orientation (especially left/right) of all images
%   - we assume that all images are in the same space and are coregistered!
% =========================================================================

%--- SETTINGS -------------------------------------------------------------
threshold_GM = 0.25;
datatype = 64; % float64, see spm_type.m
%---------------------------------------------------------------------------

calculate_GACAI = false;
target_GMactivityimage_native = '';


if ~isempty(varargin)
    if numel(varargin) >= 1
        target_GMactivityimage_native = varargin{1};
        if ~isempty(target_GMactivityimage_native) && exist(target_GMactivityimage_native, 'file')
            calculate_GACAI = true;
            fprintf('Calculating ACAI and GACAI \n');
        else
            fprintf('Calculating ACAI (target_GMactivityimage_native was empty or does not exist)\n');
            if ~isempty(target_GMactivityimage_native) && ~exist(target_GMactivityimage_native, 'file')
                fprintf('Warning: Specified target_GMactivityimage_native not found: %s\n', target_GMactivityimage_native);
            end
            target_GMactivityimage_native = ''; % Ensure it's cleared if not valid
        end
    end
else
    fprintf('Calculating ACAI \n');
end

% Define kernel
%---------------

if strcmpi(kernel_type,'Cubic')
    namekernel = ['_Cubic' num2str(kernel_size)];
    Kernel = ones(kernel_size,kernel_size,kernel_size);


elseif strcmpi(kernel_type,'Gaussian')
    sd = kernel_size/2.355;
    namekernel = ['_Gauss' num2str(kernel_size)];
    r = (round(2*kernel_size)-1)/2; 
    [x,y,z]=meshgrid(-r:r,-r:r,-r:r);
    h = exp(-(x.*x+y.*y +z.*z)/(2*sd*sd));
    sumh = sum(h(:));     
    if sumh ~= 0
        Kernel  = h/sumh;
    else
        error('Sum of Gaussian kernel is zero. Check kernel_size.');
    end


elseif strcmpi(kernel_type,'Sphere')
    namekernel = ['_Sphere' num2str(kernel_size)];
    r = (kernel_size-1)/2;
    [x_k,y_k,z_k] = meshgrid(-r:r,-r:r,-r:r); 
    Kernel = double((x_k.^2 + y_k.^2 + z_k.^2) <= r^2); 
    if sum(Kernel(:)) == 0 && kernel_size > 0 
        warning('Spherical kernel is all zeros. Check kernel_size. Defaulting to single voxel.');
        Kernel = zeros(kernel_size,kernel_size,kernel_size);
        center = ceil(kernel_size/2);
        if center > 0 && center <= kernel_size
            Kernel(center,center,center) = 1;
        elseif kernel_size == 1
             Kernel(1,1,1) = 1; 
        end
    elseif kernel_size <= 0
        error('Kernel size for Sphere must be positive.');
    end
else
    error(['No valid kernel_type. Select: ' '''Cubic''' ', ' '''Gaussian''' ' or ' '''Sphere''']);
end

spm_jobman('initcfg');
spm('defaults', 'PET');

% Get path and base names for output file naming
[pth_pet_native, pet_base_name_native, pet_ext_native] = fileparts(pet_image_native);
[orig_pth_gm, gm_base_name_native, gm_ext_native] = fileparts(gm_image_native); % Use native GM for suffix and warping

% Create a suffix from the GM image name
gm_suffix = regexprep(gm_base_name_native, '^p1', ''); % e.g., p1SUBJECT_T1 -> SUBJECT_T1
if isempty(gm_suffix) % if regexprep removed everything (e.g. name was just 'p1')
    gm_suffix = gm_base_name_native; % use full name as fallback
end
kernelinfo = [kernel_type(1) num2str(kernel_size)];


% STEP A: Warp native images to MNI space
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
clear matlabbatch;

filelist_to_warp = {pet_image_native; gm_image_native}; % Add gm_image_native to be warped
if calculate_GACAI
    filelist_to_warp{end+1,1} = target_GMactivityimage_native;
end

matlabbatch{1}.spm.spatial.normalise.write.subj.def = cellstr(deformationfield_forward);
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = filelist_to_warp;
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-90 -126 -72; 90 90 108]; 
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1]; 
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4; 
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w_mni_'; 

fprintf('Warping images to MNI space...\n');
spm_jobman('serial', matlabbatch);

% Define paths to MNI-space warped images
pet_image_mni = fullfile(pth_pet_native, ['w_mni_' pet_base_name_native pet_ext_native]);
gm_image_mni = fullfile(orig_pth_gm, ['w_mni_' gm_base_name_native gm_ext_native]); % Path to warped GM

if calculate_GACAI
    [~, gm_act_base_native, gm_act_ext_native] = fileparts(target_GMactivityimage_native);
    target_GMactivityimage_mni = fullfile(pth_pet_native, ['w_mni_' gm_act_base_native gm_act_ext_native]);
end


% STEP B: Calculate ACAI (and GACAI) images in MNI space
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
fprintf('Reading MNI-space PET image: %s\n', pet_image_mni);
[PET_mni, Vref_mni] = LCN12_read_image(pet_image_mni); 

fprintf('Reading MNI-space GM image: %s\n', gm_image_mni);
[GM_mni, ~]  = LCN12_read_image(gm_image_mni, Vref_mni); 

Mask_GM_mni  = (GM_mni > threshold_GM);

% Replace NaN by 0
PET_mni(isnan(PET_mni)) = 0;
GM_mni(isnan(GM_mni)) = 0; 

% Initializations
VI_mni           = zeros(Vref_mni.dim(1:3));
VI_LRflipped_mni = zeros(Vref_mni.dim(1:3));
AsymIndex_mni    = zeros(Vref_mni.dim(1:3));

% Calculation of ACAI map
%-------------------------
PETGM_mni = PET_mni .* GM_mni;


fprintf('Starting convolution for ACAI: PET_mni*GM_mni\n');
cPETGM_mni = convn(PETGM_mni, Kernel, 'same');

fprintf('Starting convolution for ACAI: GM_mni\n');
cGM_mni = convn(GM_mni, Kernel, 'same');

Mask_Temp1_mni = (cGM_mni > 0);
VI_mni(Mask_Temp1_mni) = cPETGM_mni(Mask_Temp1_mni) ./ cGM_mni(Mask_Temp1_mni);
VI_mni(isinf(VI_mni) | isnan(VI_mni)) = 0; 


% Flip left/right in MNI space
[ixg,iyg,izg]  = ind2sub(size(cGM_mni), find(Mask_Temp1_mni)); 
XYZ_vox   = [ixg';iyg';izg']; 
XYZ_mm = Vref_mni.mat(1:3,1:3) * XYZ_vox + repmat(Vref_mni.mat(1:3,4),1,size(XYZ_vox,2));

XYZ_mm_LRflipped = XYZ_mm;
XYZ_mm_LRflipped(1,:) = -XYZ_mm(1,:); 


XYZ_vox_LRflipped_float = inv(Vref_mni.mat(1:3,1:3)) * (XYZ_mm_LRflipped - repmat(Vref_mni.mat(1:3,4),1,size(XYZ_vox,2)));
XYZ_vox_LRflipped = round(XYZ_vox_LRflipped_float);

nr_voxels_to_flip = size(XYZ_vox_LRflipped,2);


for i = 1:nr_voxels_to_flip
    orig_vx = XYZ_vox(1,i); 
    orig_vy = XYZ_vox(2,i); 
    orig_vz = XYZ_vox(3,i); 

    flipped_vx = XYZ_vox_LRflipped(1,i);
    flipped_vy = XYZ_vox_LRflipped(2,i);
    flipped_vz = XYZ_vox_LRflipped(3,i);

    
    if flipped_vx >= 1 && flipped_vx <= Vref_mni.dim(1) && ...
       flipped_vy >= 1 && flipped_vy <= Vref_mni.dim(2) && ...
       flipped_vz >= 1 && flipped_vz <= Vref_mni.dim(3)
        VI_LRflipped_mni(flipped_vx, flipped_vy, flipped_vz) = VI_mni(orig_vx, orig_vy, orig_vz);
    end
end

Mask_calc_mni = (VI_mni ~= 0) & (VI_LRflipped_mni ~= 0) & Mask_GM_mni; 
AsymIndex_mni(Mask_calc_mni) = 200 * (VI_mni(Mask_calc_mni) - VI_LRflipped_mni(Mask_calc_mni)) ./ (VI_mni(Mask_calc_mni) + VI_LRflipped_mni(Mask_calc_mni));
AsymIndex_mni(isinf(AsymIndex_mni) | isnan(AsymIndex_mni)) = 0; 

% Define MNI-space output filename for ACAI
outputfilename_acai_mni = fullfile(pth_pet_native, ['ACAI_mni_' pet_base_name_native '_' gm_suffix '_' kernelinfo '.nii']);


fprintf('Writing MNI-space ACAI image: %s\n', outputfilename_acai_mni);
Vout_acai_mni = LCN12_write_image(AsymIndex_mni, outputfilename_acai_mni, 'MNI Asymmetry Index Image', datatype, Vref_mni);


if calculate_GACAI
    % Reset initializations for GACAI
    VI_mni_gacai           = zeros(Vref_mni.dim(1:3));
    VI_LRflipped_mni_gacai = zeros(Vref_mni.dim(1:3));
    AsymIndex_mni_gacai    = zeros(Vref_mni.dim(1:3));

    
    fprintf('Reading MNI-space PET GM activity image: %s\n', target_GMactivityimage_mni);
    [PETGMact_mni, ~] = LCN12_read_image(target_GMactivityimage_mni, Vref_mni);
    PETGMact_mni(isnan(PETGMact_mni)) = 0;

    % Calculation of GACAI map
    PETGM_for_GACAI = PETGMact_mni .* GM_mni; 

    
    fprintf('Starting convolution for GACAI: PETGMact_mni*GM_mni\n');
    cPETGM_gacai = convn(PETGM_for_GACAI, Kernel, 'same');

    
    VI_mni_gacai(Mask_Temp1_mni) = cPETGM_gacai(Mask_Temp1_mni) ./ cGM_mni(Mask_Temp1_mni);
    VI_mni_gacai(isinf(VI_mni_gacai) | isnan(VI_mni_gacai)) = 0;

    % Flip left/right 
    
    for i = 1:nr_voxels_to_flip
        orig_vx = XYZ_vox(1,i); 
        orig_vy = XYZ_vox(2,i);
        orig_vz = XYZ_vox(3,i);

        flipped_vx = XYZ_vox_LRflipped(1,i);
        flipped_vy = XYZ_vox_LRflipped(2,i);
        flipped_vz = XYZ_vox_LRflipped(3,i);
        
        if flipped_vx >= 1 && flipped_vx <= Vref_mni.dim(1) && ...
           flipped_vy >= 1 && flipped_vy <= Vref_mni.dim(2) && ...
           flipped_vz >= 1 && flipped_vz <= Vref_mni.dim(3)
            VI_LRflipped_mni_gacai(flipped_vx, flipped_vy, flipped_vz) = VI_mni_gacai(orig_vx, orig_vy, orig_vz);
        end
    end
    
    Mask_calc_mni_gacai = (VI_mni_gacai ~= 0) & (VI_LRflipped_mni_gacai ~= 0) & Mask_GM_mni;
    AsymIndex_mni_gacai(Mask_calc_mni_gacai) = 200 * (VI_mni_gacai(Mask_calc_mni_gacai) - VI_LRflipped_mni_gacai(Mask_calc_mni_gacai)) ./ (VI_mni_gacai(Mask_calc_mni_gacai) + VI_LRflipped_mni_gacai(Mask_calc_mni_gacai));
    AsymIndex_mni_gacai(isinf(AsymIndex_mni_gacai) | isnan(AsymIndex_mni_gacai)) = 0;
    
    outputfilename_gacai_mni = fullfile(pth_pet_native, ['GACAI_mni_' pet_base_name_native '_' gm_suffix '_' kernelinfo '.nii']);
    
    
    fprintf('Writing MNI-space GACAI image: %s\n', outputfilename_gacai_mni);
    Vout_gacai_mni = LCN12_write_image(AsymIndex_mni_gacai, outputfilename_gacai_mni, 'MNI GM-Only Asymmetry Index Image', datatype, Vref_mni);
end


% STEP C: Inverse warp ACAI/GACAI to native space
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
clear matlabbatch;

filelist_to_inverse_warp = {outputfilename_acai_mni};

if calculate_GACAI
   filelist_to_inverse_warp{end+1,1} = outputfilename_gacai_mni; 
end

matlabbatch{1}.spm.spatial.normalise.write.subj.def = cellstr(deformationfield_inverse);
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = filelist_to_inverse_warp;

Vpet_native_header = spm_vol(pet_image_native); 
bb_native = spm_get_bbox(Vpet_native_header);
vox_native = abs(diag(Vpet_native_header.mat(1:3,1:3)))'; 

matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-Inf -Inf -Inf
                                                          Inf Inf Inf];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 1;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w_native_'; 


fprintf('Inverse warping ACAI/GACAI images to native space...\n');
spm_jobman('serial', matlabbatch);

% Define paths to temporary native-space (inverse-warped) images
[~, acai_mni_base, acai_mni_ext] = fileparts(outputfilename_acai_mni);
temp_native_acai_path = fullfile(pth_pet_native, ['w_native_' acai_mni_base acai_mni_ext]);


if calculate_GACAI
    [~, gacai_mni_base, gacai_mni_ext] = fileparts(outputfilename_gacai_mni);
    temp_native_gacai_path = fullfile(pth_pet_native, ['w_native_' gacai_mni_base gacai_mni_ext]);
end

% STEP D: Clean up directory and rename final files
%+++++++++++++++++++++++++++++++++++++++++++++++++++
mni_space_folder = fullfile(pth_pet_native, 'MNI_space');

if ~exist(mni_space_folder, 'dir')
    mkdir(mni_space_folder);
end

% Move MNI-space files created by this script

fprintf('Moving MNI-space files to: %s\n', mni_space_folder);
movefile(pet_image_mni, fullfile(mni_space_folder, [pet_base_name_native pet_ext_native])); 
movefile(gm_image_mni, fullfile(mni_space_folder, [gm_base_name_native gm_ext_native])); % Move the warped GM to MNI_space
movefile(outputfilename_acai_mni, fullfile(mni_space_folder, ['ACAI_mni_' pet_base_name_native '_' gm_suffix '_' kernelinfo '.nii']));


if calculate_GACAI
    [~, gm_act_base_native_cleanup, gm_act_ext_native_cleanup] = fileparts(target_GMactivityimage_native); 
    movefile(target_GMactivityimage_mni, fullfile(mni_space_folder, [gm_act_base_native_cleanup gm_act_ext_native_cleanup])); 
    movefile(outputfilename_gacai_mni, fullfile(mni_space_folder, ['GACAI_mni_' pet_base_name_native '_' gm_suffix '_' kernelinfo '.nii']));
end

% Rename final native space ACAI/GACAI images
final_native_acai_name = fullfile(pth_pet_native, ['ACAI_native_' pet_base_name_native '_' gm_suffix '_' kernelinfo '.nii']);
movefile(temp_native_acai_path, final_native_acai_name);

fprintf('Final native ACAI image: %s\n', final_native_acai_name);


if calculate_GACAI
    final_native_gacai_name = fullfile(pth_pet_native, ['GACAI_native_' pet_base_name_native '_' gm_suffix '_' kernelinfo '.nii']);
    movefile(temp_native_gacai_path, final_native_gacai_name);
    fprintf('Final native GACAI image: %s\n', final_native_gacai_name);
end


fprintf('Processing complete.\n');

end