function LCN12_calc_ACAI_GACAI(target_image,MRI_image,kernel_type,kernel_size,target_GMactivityimage)
% [ACAIimage,GACAIimage] = LCN12_calc_ACAI_GACAI(target_image,MRI_image,kernel_type,kernel_size,target_GMactivityimage)
% Calculating the AsymmetryIndex(AI) using prior anatomical information
%
% Input:
%   Target image: 
%      full name of the image whose asymmetry index needs to be calculated (for example a FDG image)
%   MRI_image:
%      full name of the 3D high resolution T1 weighted MRI image
%   kernel_type:
%      kernel type: 'Gaussian', 'Cubic' or 'Sphere'
%      Default kernel: 'Gaussian'.
%
%   kernel_size:
%      Size of the kernel (in voxels).
%      'Gaussian': FWHM
%      'Cubic': cube side
%      'Sphere': sphere diameter   
%   Optional: target_GMactivitymap = image of only the GM contribution.
%             This image is created during the AMAP procedure. We assume
%             that this image is in the same space and has the same
%             dimensions and voxel sizes as target_image.
%
% The results will be written to file.
%   AsymmetryIndex images of the target image (ACAIimage and optional GACAIimage).
%   The description of the method to calculate the ACAIimage is described
%   in Zhou et al. NeuroImage 44 (2009), 35-42
%   The GACAIimage is an improved method taking into account only the GM
%   contribution in every voxel.
%
% Important:
%  - requires the installation of SPM12 - http://www.fil.ion.ucl.ac.uk/spm/
%  - the path of SPM and this routine should be included in the Matlab path
%
% The package contains software (SPM12) developed under the auspices of The
% Wellcome Department of Imaging Neuroscience, a department of the
% Institute of Neurology at University College London. The copyright of
% this software remains with that of SPM12, see
% http://www.fil.ion.ucl.ac.uk/spm/.
%    
% This routine is supplied as is. 
% Comments or questions can be send to:
% Patrick.Dupont@kuleuven.be
%
% IMPORTANT REMARKS: 
%   - this is research software.
%   - always check the orientation (especially left/right) of all images
%   - we assume that all images are in the same space and have the same
%   dimensions!
%
%__________________________________________________________________________
%
% authors: Lin Zhou, Kathleen Vunckx, Patrick Dupont
% date:    July 2017
% history: 
%__________________________________________________________________________
% @(#)LCN12_calc_ACAI_GACAI.m       v0.1          last modified: 2017/07/18

%--- SETTINGS -------------------------------------------------------------
threshold_GM = 0.25;
% fill in the spm_path for example C:\MATLAB\R2023b\toolbox\spm12
spm_path = '';
datatype = 64; % float64, see spm_type.m
%---------------------------------------------------------------------------

if nargin < 5
   fprintf('calculating ACAI \n');
elseif nargin == 5
   fprintf('calculating ACAI and GACAI \n');
end

% define kernel
%---------------
if strcmp(kernel_type,'Cubic')
   namekernel = ['_Cubic' num2str(kernel_size)];
   Kernel = ones(kernel_size,kernel_size,kernel_size);
elseif strcmp(kernel_type,'Gaussian')
   sd = kernel_size/2.355;
   namekernel = ['_Gauss' num2str(kernel_size)];
   r = (round(2*kernel_size)-1)/2; % sampling in kernel is twice the kernel_size
   [x,y,z]=meshgrid(-r:r,-r:r,-r:r);
   h = exp(-(x.*x+y.*y +z.*z)/(2*sd*sd));
   sumh = sum(h(:));     
   if sumh ~= 0
      Kernel  = h/sumh;
   end
elseif strcmp(kernel_type,'Sphere')
   namekernel = ['_Sphere' num2str(kernel_size)];
   r = (kernel_size-1)/2;
   [x,y,z] = meshgrid(-r:r,-r:r,-r:r);
   h = x.*x+y.*y+z.*z;
   for i=1:kernel_size
       for j=1:kernel_size 
           for k=1:kernel_size
               if h(i,j,k) > (r^2)
                  h(i,j,k)=0; 
               else 
                  h(i,j,k)=1;
               end 
           end 
       end
   end      
   Kernel = h;
else
   disp(['No valid kernel_type. Select: ' '''Cubic''' ', ' '''Gaussian''' ' or ' '''Sphere'''])
end

spm_jobman('initcfg');

% STEP 1: coregistration of the MRI and the target image
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
matlabbatch{1}.spm.spatial.coreg.estimate.ref = cellstr(target_image);
matlabbatch{1}.spm.spatial.coreg.estimate.source = cellstr(MRI_image);
matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

spm('defaults', 'PET');
spm_jobman('serial', matlabbatch);

clear matlabbatch

% STEP 2: segmentation of the MRI
%++++++++++++++++++++++++++++++++
matlabbatch{1}.spm.spatial.preproc.channel.vols = cellstr(MRI_image);
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = cellstr(fullfile(fullfile(spm_path,'tpm'),'TPM.nii,1'));
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = cellstr(fullfile(fullfile(spm_path,'tpm'),'TPM.nii,2'));
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = cellstr(fullfile(fullfile(spm_path,'tpm'),'TPM.nii,3'));
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = cellstr(fullfile(fullfile(spm_path,'tpm'),'TPM.nii,4'));
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = cellstr(fullfile(fullfile(spm_path,'tpm'),'TPM.nii,5'));
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = cellstr(fullfile(fullfile(spm_path,'tpm'),'TPM.nii,6'));
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];
spm_jobman('serial', matlabbatch);

% get the name of the deformation fields
%---------------------------------------
[pth,filename,~] = fileparts(MRI_image);
filename1 = ['y_' filename '.nii'];
filename2 = ['iy_' filename '.nii'];
deformationfield_forward = fullfile(pth,filename1);
deformationfield_inverse = fullfile(pth,filename2); 

% STEP 3: warping of target_image image, GM and target_GMactivityimage 
%         to MNI based upon previous step
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
clear matlabbatch

filename_GM = fullfile(pth,['c1' filename '.nii']);
filelist = {};
filelist{1,1} = target_image;
filelist{2,1} = filename_GM;
if nargin == 5
   filelist{3,1} = target_GMactivityimage;
end

matlabbatch{1}.spm.spatial.normalise.write.subj.def = cellstr(deformationfield_forward);
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = filelist;
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-90 -126 -71
                                                          90 90 109];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
spm_jobman('serial', matlabbatch);

% STEP 4: calculate ACAI and GACAI images
%++++++++++++++++++++++++++++++++++++++++
[pth,filename,ext] = fileparts(target_image);
filename_wPET = fullfile(pth,['w' filename ext ',1']);
[pth,filename,ext] = fileparts(filename_GM);
filename_wGM = fullfile(pth,['w' filename ext ',1']);
if nargin == 5
   [pth,filename,ext] = fileparts(target_GMactivityimage);
   filename_wGMact = fullfile(pth,['w' filename ext ',1']);
end
 
% read the images
%----------------
[PET,Vref] = LCN12_read_image(filename_wPET);
[GM,~]  = LCN12_read_image(filename_wGM,Vref);

Mask_GM  = (GM > threshold_GM);

% replace NaN by 0 (usually at the border of the image)
%-------------------------------------------------------
PET(isnan(PET)) = 0;

% initializations 
VI           = zeros(Vref.dim(1:3));
VI_LRflipped = zeros(Vref.dim(1:3));
AsymIndex    = zeros(Vref.dim(1:3));

% calculation of ACAI map
%-------------------------
PETGM = PET.*GM;

fprintf('starting convolution PET\n');
cPETGM = convn(PETGM, Kernel, 'same');
fprintf('starting convolution GM\n');
cGM = convn(GM, Kernel, 'same');

Mask_Temp1 = (cGM > 0);
VI(Mask_Temp1) = cPETGM(Mask_Temp1)./cGM(Mask_Temp1);

% flip left/right
[ixg,iyg,izg]  = ind2sub(size(cGM),find(cGM~=0)); % choose voxels within mask
XYZ   = [ixg';iyg';izg']; % 3 x N matrix with coordinates of the voxels selected within the mask
XYZmm = (Vref.mat(1:3,1:3))*XYZ +repmat(Vref.mat(1:3,4),1,size(XYZ,2)); % world coordinates (eg MNI coordinates)

XYZmm_LRflipped = XYZmm;
XYZmm_LRflipped(1,:) = -XYZmm(1,:);
XYZ_LRflipped = round(inv(Vref.mat(1:3,1:3))*(XYZmm_LRflipped - repmat(Vref.mat(1:3,4),1,size(XYZ,2))));
nr_voxels = size(XYZ_LRflipped,2);
for i = 1:nr_voxels
    VI_LRflipped(XYZ_LRflipped(1,i),XYZ_LRflipped(2,i),XYZ_LRflipped(3,i)) = VI(XYZ(1,i),XYZ(2,i),XYZ(3,i));
end

Mask = (VI ~= 0).*(VI_LRflipped~= 0).*Mask_GM; 

AsymIndex(Mask>0) = 200 * (VI(Mask>0) - VI_LRflipped(Mask>0))./(VI(Mask>0) + VI_LRflipped(Mask>0));

%--------------------------------------------
[pth, filename, ~] = fileparts(target_image);
[pth2, filename2, ~] = fileparts(MRI_image);

% Ensure filename2 is long enough
if length(filename2) >= 14
    suffix = filename2(end-13:end);
else
    % Handle the case where filename2 is too short
    suffix = filename2; % or use a default value
end

kernelinfo = [kernel_type(1) num2str(kernel_size)];
outputfilename = fullfile(pth, ['wACAI_' filename '_' suffix '_' kernelinfo '.nii']);


% Write ACAI image
%------------------
Vout = LCN12_write_image(AsymIndex,outputfilename,'AsymmetryIndex Image',datatype,Vref);

if nargin == 5
   % initializations 
   VI           = zeros(Vref.dim(1:3));
   VI_LRflipped = zeros(Vref.dim(1:3));
   AsymIndex    = zeros(Vref.dim(1:3));

   % read the image
   %---------------
   [PETGMact,Vref] = LCN12_read_image(filename_wGMact,Vref);


   % replace NaN by 0 (usually at the border of the image)
   %-------------------------------------------------------
   PETGMact(isnan(PETGMact)) = 0;


   % calculation of GACAI map
   %-------------------------
   PETGM = PETGMact.*GM;

   fprintf('starting convolution PETGMact\n');
   cPETGM = convn(PETGM, Kernel, 'same');

   VI(Mask_Temp1) = cPETGM(Mask_Temp1)./cGM(Mask_Temp1);

   % flip left/right
   for i = 1:nr_voxels
       VI_LRflipped(XYZ_LRflipped(1,i),XYZ_LRflipped(2,i),XYZ_LRflipped(3,i)) = VI(XYZ(1,i),XYZ(2,i),XYZ(3,i));
   end

   Mask = (VI ~= 0).*(VI_LRflipped~= 0).*Mask_GM; 

   AsymIndex(Mask>0) = 200.*(VI(Mask>0)-VI_LRflipped(Mask>0))./(VI(Mask>0) + VI_LRflipped(Mask>0));
    
   % Set the pathname for the output ACAI image
   %--------------------------------------------
   outputfilename2 = fullfile(pth,['wGACAI_' filename '.nii']);

   % Write ACAI image
   %------------------
   Vout2 = LCN12_write_image(AsymIndex,outputfilename2,'AsymmetryIndex Image',datatype,Vref);
   Vout_unmasked = LCN12_write_image(AsymIndexUnmask, outputfilename_unmasked, 'Unmasked AsymmetryIndex Image', datatype, Vref);

end

% STEP 5: inverse warping to bring result back into patient space
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
clear matlabbatch

filelist = {};
filelist{1,1} = outputfilename;
if nargin == 5
   filelist{2,1} = outputfilename2; 
end
matlabbatch{1}.spm.spatial.normalise.write.subj.def = cellstr(deformationfield_inverse);
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = filelist;
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-Inf -Inf -Inf
                                                          Inf Inf Inf];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [NaN NaN NaN];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 1;
spm_jobman('serial', matlabbatch);

% STEP 6: clean up directory
%+++++++++++++++++++++++++++
cd(pth)
mkdir('MNI_space')
% move warped PET
[pth,filename,ext] = fileparts(filename_wPET);
if strfind(ext,'.img') == 1  
   movefile(fullfile(pth,[filename '.img']),fullfile(fullfile(pth,'MNI_space'),[filename '.img']));
   movefile(fullfile(pth,[filename '.hdr']),fullfile(fullfile(pth,'MNI_space'),[filename '.hdr']));
elseif strfind(ext,'.nii') == 1
   movefile(fullfile(pth,[filename '.nii']),fullfile(fullfile(pth,'MNI_space'),[filename '.nii']));
end

% move warped GM map
[pth,filename,ext] = fileparts(filename_wGM);
if strfind(ext,'.img') == 1  
   movefile(fullfile(pth,[filename '.img']),fullfile(fullfile(pth,'MNI_space'),[filename '.img']));
   movefile(fullfile(pth,[filename '.hdr']),fullfile(fullfile(pth,'MNI_space'),[filename '.hdr']));
elseif strfind(ext,'.nii') == 1
   movefile(fullfile(pth,[filename '.nii']),fullfile(fullfile(pth,'MNI_space'),[filename '.nii']));
end

% move warped ACAI map
[pth,filename,ext] = fileparts(outputfilename);
if strfind(ext,'.img') == 1  
   movefile(fullfile(pth,[filename '.img']),fullfile(fullfile(pth,'MNI_space'),[filename '.img']));
   movefile(fullfile(pth,[filename '.hdr']),fullfile(fullfile(pth,'MNI_space'),[filename '.hdr']));
elseif strfind(ext,'.nii') == 1
   movefile(fullfile(pth,[filename '.nii']),fullfile(fullfile(pth,'MNI_space'),[filename '.nii']));
end

% move deformation fields
[pth,filename,ext] = fileparts(deformationfield_forward);
movefile(deformationfield_forward,fullfile(fullfile(pth,'MNI_space'),[filename '.nii']));
[pth,filename,ext] = fileparts(deformationfield_inverse);
movefile(deformationfield_inverse,fullfile(fullfile(pth,'MNI_space'),[filename '.nii']));

% move warped PETGMact
if nargin == 5
   [pth,filename,ext] = fileparts(filename_wGMact);
   if strfind(ext,'.img') == 1  
      movefile(fullfile(pth,[filename '.img']),fullfile(fullfile(pth,'MNI_space'),[filename '.img']));
      movefile(fullfile(pth,[filename '.hdr']),fullfile(fullfile(pth,'MNI_space'),[filename '.hdr']));
   elseif strfind(ext,'.nii') == 1
      movefile(fullfile(pth,[filename '.nii']),fullfile(fullfile(pth,'MNI_space'),[filename '.nii']));
   end
   [pth,filename,ext] = fileparts(outputfilename2);
   if strfind(ext,'.img') == 1  
      movefile(fullfile(pth,[filename '.img']),fullfile(fullfile(pth,'MNI_space'),[filename '.img']));
      movefile(fullfile(pth,[filename '.hdr']),fullfile(fullfile(pth,'MNI_space'),[filename '.hdr']));
   elseif strfind(ext,'.nii') == 1
      movefile(fullfile(pth,[filename '.nii']),fullfile(fullfile(pth,'MNI_space'),[filename '.nii']));
   end
end

% rename ACAI and GACAI
[pth,filename,~] = fileparts(outputfilename);
movefile(fullfile(pth,['w' filename '.nii']),fullfile(pth,[filename(2:end) '.nii']));
if nargin == 5
   [pth,filename,~] = fileparts(outputfilename2);
   movefile(fullfile(pth,['w' filename '.nii']),fullfile(pth,[filename(2:end) '.nii']));
end
if nargin == 4
   
end