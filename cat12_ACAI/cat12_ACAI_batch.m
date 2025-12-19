function matlabbatch = cat12_ACAI_batch(t1_image_file)
% =========================================================================
% cat12_ACAI_batch.m
% using the specific parameters required for the ACAI workflow.
% This version:
%   - auto-detects SPM/CAT12 install paths
%   - enables CAT12 expert mode when available
%   - resolves shooting template for both old/new CAT12 layouts
% All processing options/values MATCH original EXACTLY.
% =========================================================================

% 1. check string
if iscell(t1_image_file)
    t1_image_file = t1_image_file{1};
end

% 2. Append ',1' if it is missing
if ~endsWith(t1_image_file, ',1')
    t1_image_file = [t1_image_file ',1'];
end

% 3. Wrap back into a cell for the SPM batch structure
t1_image_file = {t1_image_file};

% --- Resolve SPM/CAT12 and enable expert mode ---
if isempty(which('spm'))
    error('SPM12 not found on MATLAB path. Add it first, e.g., addpath(''C:\path\to\spm12'');');
end
spm_dir = spm('Dir');
cat_dir = fullfile(spm_dir,'toolbox','cat12');
if ~isfolder(cat_dir)
    error('CAT12 not found under SPM: %s', cat_dir);
end

if ~isempty(which('cat12'))
    try
        cat12('expert');
        fprintf('[INFO] CAT12 expert mode enabled.\n');
    catch ME
        fprintf(2, '[ERROR] No expert mode found (cat12(''expert'') failed: %s). Continuing without expert mode. CAT12 is probably not CAT12.9.\n', ME.message);
    end
else
    fprintf(2, '[ERROR] CAT12 not found on MATLAB path. Ensure spm12/toolbox/cat12 is installed.\n');
    error('CAT12 missing');
end


% --- TPM and shooting template ---
tpm = fullfile(spm_dir,'tpm','TPM.nii');
if ~isfile(tpm), error('TPM not found: %s', tpm); end

tmpl_new = fullfile(cat_dir,'templates_MNI152NLin2009cAsym','Template_0_GS.nii');
tmpl_old = fullfile(cat_dir,'templates_volumes','Template_0_IXI555_MNI152_GS.nii');
if isfile(tmpl_new)
    shooting_tpm = tmpl_new;
elseif isfile(tmpl_old)
    shooting_tpm = tmpl_old;
else
    error('CAT12 shooting template not found:\n  %s\n  %s', tmpl_new, tmpl_old);
end

% --- CAT12 Batch Definition ---
matlabbatch{1}.spm.tools.cat.estwrite.data = t1_image_file;
matlabbatch{1}.spm.tools.cat.estwrite.data_wmh = {''};
matlabbatch{1}.spm.tools.cat.estwrite.nproc = 1;
matlabbatch{1}.spm.tools.cat.estwrite.useprior = '';
matlabbatch{1}.spm.tools.cat.estwrite.opts.tpm = {tpm};
matlabbatch{1}.spm.tools.cat.estwrite.opts.affreg = 'mni';
matlabbatch{1}.spm.tools.cat.estwrite.opts.biasstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.opts.accstr = 0.5;

matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.restypes.optimal = [1 0.3];
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.setCOM = 1;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.APP = 1070;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.affmod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.NCstr = -Inf;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.LASstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.LASmyostr = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.gcutstr = 2;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.cleanupstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.BVCstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.WMHC = 2;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.SLC = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.mrf = 1;

matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.regmethod.shooting.shootingtpm = {shooting_tpm};
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.regmethod.shooting.regstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.vox = 1.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.bb = 12;

% Note: these "surface" parameters do not enable surface outputs by themselves.
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.pbtres = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.pbtmethod = 'pbtsimple';
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.SRP = 22;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.reduce_mesh = 1;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.vdist = 2;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.scale_cortex = 0.7;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.add_parahipp = 0.1;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.close_parahipp = 1;

matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.experimental = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.new_release = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.lazy = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.ignoreErrors = 1;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.verb = 2;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.print = 2;

matlabbatch{1}.spm.tools.cat.estwrite.output.BIDS.BIDSno = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.surface = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.surf_measures = 1;

matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.neuromorphometrics = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.lpba40 = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.cobra = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.hammers = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.thalamus = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.thalamic_nuclei = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.suit = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.ibsr = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.aal3 = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.mori = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.anatomy3 = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.julichbrain = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.Tian_Subcortex_S4_7T = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.Schaefer2018_100Parcels_17Networks_order = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.Schaefer2018_200Parcels_17Networks_order = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.Schaefer2018_400Parcels_17Networks_order = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.Schaefer2018_600Parcels_17Networks_order = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.ownatlas = {''};

matlabbatch{1}.spm.tools.cat.estwrite.output.GM.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.warped = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.mod = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.dartel = 0;

matlabbatch{1}.spm.tools.cat.estwrite.output.WM.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.warped = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.mod = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.dartel = 0;

matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.warped = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.mod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.dartel = 0;

matlabbatch{1}.spm.tools.cat.estwrite.output.ct.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ct.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ct.dartel = 0;

matlabbatch{1}.spm.tools.cat.estwrite.output.pp.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.pp.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.pp.dartel = 0;

matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.mod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.dartel = 0;

matlabbatch{1}.spm.tools.cat.estwrite.output.SL.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.mod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.dartel = 0;

matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.mod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.dartel = 0;

matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.native = 1;

matlabbatch{1}.spm.tools.cat.estwrite.output.label.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.label.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.label.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.labelnative = 1;

matlabbatch{1}.spm.tools.cat.estwrite.output.bias.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.warped = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.dartel = 0;

matlabbatch{1}.spm.tools.cat.estwrite.output.las.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.las.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.las.dartel = 0;

matlabbatch{1}.spm.tools.cat.estwrite.output.jacobianwarped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.warps = [1 1];
matlabbatch{1}.spm.tools.cat.estwrite.output.rmat = 0;

end
