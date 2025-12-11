%asks for input images before starting the function LCN12_calc_ACAI_GACAI



% target_image
target_image = spm_select(1, 'IMAGE', {'Select the FDG PET image'});

% MRI_image
MRI_image = spm_select(1, 'IMAGE', {'Select the MRI image'});


%kernel_type
kernel_type = spm_input('Kernel_type','+1','Gaussian|Cubic|Sphere',[1,2,3],1);
    if kernel_type ==1
        kernel_type= 'Gaussian';
    elseif kernel_type == 2
        kernel_type = 'Cubic';
    elseif kernel_type ==3
            kernel_type = 'Sphere';
    end;


%kernel_size
kernel_size = spm_input('kernel_size (voxels)','0',[],9); %filenames waarin kernel niet gedefinieerd is, hebben kernelsize6. paper Zhou zegt 9

%target_GMactivityimage & run ACAI script patrick!!!!
GMactivityimage = spm_input('Do you want to select a GM activityimage?','+1','YES|NO',[1,2],2);
if GMactivityimage == 1
    target_GMactivityimage = spm_select(1, 'IMAGE', {'Select the GMactivityimage'});
    LCN12_calc_ACAI_GACAI(target_image,MRI_image,kernel_type,kernel_size,target_GMactivityimage)
elseif GMactivityimage == 2
    LCN12_calc_ACAI_GACAI(target_image,MRI_image,kernel_type,kernel_size)
end;





