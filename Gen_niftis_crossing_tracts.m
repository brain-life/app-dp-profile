function [] = Gen_niftis_crossing_tracts(info, tract_name_1, tract_name_2, other_tracts)
%% INPUT: two crossing tracts (AFQ)
%%  Number    Name            Full Name
%   1         Thal_Rad_L      Left Thalamic Radiation
%   2         Thal_Rad_R      Right Thalamic Radiation
%   3         CST_L           Left Corticospinal 
%   4         CST_R           Right Corticospinal 
%   5         Cing_L          Left Cingulum Cingulate
%   6         Cing_R          Right Cingulum Cingulate
%   7         Hipp_L          Left Cingulum Hippocampus
%   8         Hipp_R          Right Cingulum Hippocampus
%   9         Call_Maj        Callosum Forceps Major
%   10        Call_Min        Callosum Forceps Minor
%   11        IFOF_L          Left IFOF  
%   12        IFOF_R          Right IFOF
%   13        ILF_L           Left ILF
%   14        ILF_R           Right ILF
%   15        SLF_L           Left SLF
%   16        SLF_R           Right SLF
%   17        Unc_L           Left Uncinate
%   18        Unc_R           Right Uncinate
%   19        ARC_L           Left Arcuate
%   20        ARC_R           Right Arcuate

%% INPUT: two crossing tracts (Dan tracts)
%%  Number    Name            Full Name
%   1       Thal_Rad_L        Left Thalamic Radiation
%   2       Thal_Rad_R        Right Thalamic Radiation
%   3       CST_L             Left Corticospinal
%   4       CST_R             Right Corticospinal
%   5       Cing_L            Left Cingulum Cingulate
%   6       Cing_R            Right Cingulum Cingulate
%   7       Hipp_L            Left Cingulum Hippocampus
%   8       Hipp_R            Right Cingulum Hippocampus
%   9       Call_Maj          Callosum Forceps Major
%   10      Call_Min          Callosum Forceps Minor
%   11      IFOF_L            Left IFOF
%   12      IFOF_R            Right IFOF
%   13      Unc_L             Left Uncinate
%   14      Unc_R             Right Uncinate
%   15      ARC_L             Left Arcuate
%   16      ARC_R             Right Arcuate
%   17      VOF_L             Left VOF
%   18      VOF_R             Right VOF
%   19      pARC_L            Left pArc
%   20      pARC_R            Right pArc
%   21      TPC_L             Left TPC
%   22      TPC_R             Right TPC
%   23      MdLF-SPL_L        Left MdLF-SPL
%   24      MdLF-SPL_R        Right MdLF-SPL
%   25      MdLF-Ang_L        Left MdLF-Ang
%   26      MdLF-Ang_R        Right MdLF-Ang
%   27      Meyer_L           Left Meyer
%   28      Meyer_R           Right Meyer
%   29      Baum_L            Left Baum
%   30      Baum_R            Right Baum
%   31      SLF1_L            Left SLF1
%   32      SLF1_R            Right SLF1
%   33      SLF2_L            Left SLF2
%   34      SLF2_R            Right SLF2
%   35      SLF3_L            Left SLF3
%   36      SLF3_R            Right SLF3
%   37      ILF_L             Left ILF
%   38      ILF_R             Right ILF



%% Get tract numbers
other_tract_number = [];
if info.segmentation_type == 'AFQ'
    tract_number_1 = Get_tract_number(tract_name_1);
    tract_number_2 = Get_tract_number(tract_name_2);
    for i=1:size(other_tracts,2)
        other_tract_number = [other_tract_number, Get_tract_number(other_tracts{i})];
    end
else
    tract_number_1 = Get_tract_number_Dan(tract_name_1);
    tract_number_2 = Get_tract_number_Dan(tract_name_2);
    for i=1:size(other_tracts,2)
        other_tract_number = [other_tract_number, Get_tract_number(other_tracts{i})];
    end
end

%% Set the Path for the output
dataOutputPath = info.output.niftis;

%% load fe structure
%FileName = deblank(ls(fullfile(dataInputPath,strcat('fe_*.mat'))));
FileName = info.input.optimal;
load(FileName);
%dwiFile = fe.path.dwifile;
dwiFile = info.input.dwi_path;

%% Save .b file (mrtrix format)
%bvecsFile = strcat(dwiFile(1:end-6),'bvecs');
%bvalsFile = strcat(dwiFile(1:end-6),'bvals');

%% TODO - is this necessary?
%% Generate nifti for original data 
ni = niftiRead(dwiFile);
fName = fullfile(dataOutputPath,strcat('original.nii.gz'));
niftiWrite(ni,fName);

%% Generate nifti using the predicted diffusion signal based on the model fit to the original data
% fe is already loaded
name = fullfile(dataOutputPath,strcat('pred_full.nii.gz'));
coords = fe.roi.coords; % Get the coordinates of the nodes in each voxel of the connectome
%dwi = feGet(fe, 'dwi'); % load dwi structure
dwi = dwiLoad(info.input.dwi_path); % load dwi structure
diff_signal = feGet(fe,'pred full');
Generate_nifti(ni,name,coords,dwi,diff_signal);

%% Load classification file
load(info.input.classification_path);

%% Generate niftis for single tracts prediction
tract_set   = [tract_number_1, tract_number_2];
tract_names = {tract_name_1, tract_name_2};
for i=1:length(tract_set)
    Gen_nifti_single_tract(fe,classification,tract_set(i),tract_names{i}, dataOutputPath, ni, info)
end

%% Generate niftis for crossing tracts prediction
Gen_nifti_crossing_tracts(fe,classification, tract_number_1, tract_number_2, other_tract_number, ...
                            strcat(tract_name_1, '_', tract_name_2,'_new'), dataOutputPath, ni, info)

end

function [] = Gen_nifti_single_tract(fe,classification,tract,fgName,dataOutputPath,ni,info)

name = fullfile(dataOutputPath,strcat(fgName,'.nii.gz'));
coords = fe.roi.coords; % Get the coordinates of the nodes in each voxel of the connectome
dwi = dwiLoad(info.input.dwi_path); % load dwi structure
fibers = find(classification.index == tract);
diff_signal = feGet(fe,'pred tract',fibers);
diff_signal(diff_signal==0) = NaN;
Generate_nifti(ni,name,coords,dwi,diff_signal);
end

function [] = Gen_nifti_crossing_tracts(fe,classification,tract1,tract2,other_tract_number,fgName,dataOutputPath,ni,info)

name = fullfile(dataOutputPath,strcat(fgName,'.nii.gz'));
coords = fe.roi.coords; % Get the coordinates of the nodes in each voxel of the connectome
dwi = dwiLoad(info.input.dwi_path); % load dwi structure
fibers1 = find(classification.index == tract1);
fibers2 = find(classification.index == tract2);
fibers = [fibers1; fibers2];
for i=1:length(other_tract_number)
    fibers = [fibers ;find(classification.index == other_tract_number(i))];
end
diff_signal = feGet(fe,'pred tract',fibers);
diff_signal(diff_signal==0) = NaN;
Generate_nifti(ni,name,coords,dwi,diff_signal);
end

function [] = Generate_nifti(ni_in,name,coords,dwi,dwisignal)

ni_out = ni_in;
ni_out.fname = name;

bvals = dwi.bvals;
indexes = find(bvals~=0);
b0indexes = find(bvals==0);

% Copy original S0 values
b0_data = nan(size(b0indexes,1),size(coords,1));
for ivx = 1:size(coords,1)
    b0_data(:,ivx) = ni_out.data(coords(ivx,1),coords(ivx,2),coords(ivx,3),b0indexes);
end
ni_out.data = nan(size(ni_in.data));

% Replace Nans with b0_data
ni_out.data = feReplaceImageValues(ni_out.data,b0_data,coords,b0indexes);

% Replace Nans with dw_vals
ni_out.data = feReplaceImageValues(ni_out.data,dwisignal,coords,indexes);

% save nifti to disk
niftiWrite(ni_out,ni_out.fname);
end
