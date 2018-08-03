function [] = Gen_tract_profiles_single(info, tract_name, std_parameter, measure)
%% measure = 'FA', 'MD','AD' or 'RD'
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


%addpath(genpath(info.repo.mba));

dataPath = info.output.niftis;

%% load fe_structure
%FileName = deblank(ls(fullfile(dataPath,strcat('fe_*.mat'))));
%load(FileName);
load(info.input.optimal);

%% Load classification file
ClassFileName = info.input.classification_path;
load(ClassFileName);
%classification.index = class.index;
%classification.names = class.names;


%% Insert classification into fe structure 
%ind_tracts = find(classification.index);
%classification.index = classification.index(ind_tracts);
fe = feSet(fe,'tracts_info',classification); % include tract indices in fe structure

% Obtain tract numbers
if info.segmentation_type == 'AFQ'
    tract = Get_tract_number(tract_name);
else
    tract = Get_tract_number_Dan(tract_name);
end

% Set parameters
%std_parameter = 3;
nameroot = 'nosub';
Nnodes = 50;

if isfield(fe.life.fit, 'weights')
    ind_nnz = find(fe.life.fit.weights);
else
    ind_nnz = unique(fe.life.M.Phi.subs(:,3)); % find indices of nnz fascicles
end

%% Extract fibers tract1
fgName = fe.life.M.tracts{tract}.name; fgName = strrep(fgName,' ','');

ind_tract = fe.life.M.tracts{tract}.ind;
ind_tract_nnz = intersect(ind_tract,ind_nnz);

fgTract = fe.fg.fibers(ind_tract_nnz);
fgex = MyfgCreate_img('name', fgName, 'colorRgb', [1 0 0], 'fibers', fgTract);
% clean fibers with mba
fgcx = mbaComputeFibersOutliers(fgex, std_parameter, std_parameter, 100, 'mean');

%% Compute profile tract1 using measure (FA,MD,etc) based on  tract ONLY
file = deblank(ls(char(fullfile(dataPath,strcat(measure,'s'),strcat(strcat(measure,'_'),tract_name,'.nii.gz')))));
famp = niftiRead(file);
[Meas_tract, SuperFiber, ~, ~] = Compute_FA_AlongFG(fgcx, famp, [], [], Nnodes);

%% Compute tract profile using measure (FA,MD,etc)FA based on original
fileOrig = deblank(ls(char(fullfile(dataPath,strcat(measure,'s'),strcat(measure,'_original.nii.gz')))));
fampOrig = niftiRead(fileOrig);
[Meas_tract_orig, ~]= Compute_FA_AlongFG(fgcx, fampOrig, [], [], Nnodes);

%% Compute tract profile using measure (FA,MD,etc) based on prediction
filePred = deblank(ls(char(fullfile(dataPath,strcat(measure,'s'),strcat(measure,'_pred_full.nii.gz')))));
fampPred = niftiRead(filePred);
[Meas_tract_pred, ~]= Compute_FA_AlongFG(fgcx, fampPred, [], [], Nnodes);

%% Plot tract profile
Gen_profile_plot_single(Meas_tract,'r',Meas_tract_orig,'k', Meas_tract_pred,'y',tract_name, 10, Nnodes, measure)
saveas(gcf, strcat('./results/figures/',measure,'/',measure,'_profile_',tract_name,'.fig'));
saveas(gcf, strcat('./results/figures/',measure,'/',measure,'_profile_',tract_name,'.pdf'));
saveas(gcf, strcat('./results/figures/',measure,'/',measure,'_profile_',tract_name,'.png'));


profiles_data.tract_name = tract_name;
profiles_data.tract = Meas_tract;
profiles_data.tract_orig = Meas_tract_orig;
profiles_data.tract_pred = Meas_tract_pred;

save(strcat('./results/',measure,'_',tract_name,'.mat'), 'profiles_data')
end


function [] = Gen_profile_plot_single(FA_tract, clr1, FA_tract_orig, clrorig, FA_pred, clrp,tract_name, s, Nnodes, measure)

N = size(FA_tract,1);
figure
hold on
if ~isempty(FA_tract)
    h1 = shadedErrorBar(1:Nnodes,nanmean(FA_tract,1),s*nanstd(FA_tract)/sqrt(N),'lineprops',clr1);
end
if ~isempty(FA_tract_orig)
    h3 = shadedErrorBar(1:Nnodes,nanmean(FA_tract_orig,1),s*nanstd(FA_tract_orig)/sqrt(N),'lineprops',clrorig);
end
if ~isempty(FA_pred)
    h4 = shadedErrorBar(1:Nnodes,nanmean(FA_pred,1),s*nanstd(FA_pred)/sqrt(N),'lineprops',clrp);
end

%
%legend([h1.mainLine, h2.mainLine],tract_name1,strcat(tract_name1,'+',tract_name2))
if ~isempty(FA_tract)&&isempty(FA_tract_orig)
   legend([h1.mainLine],tract_name)
elseif isempty(FA_tract)
   legend([h3.mainLine, h3.mainLine],'Pred full','Orig')
else
   legend([h1.mainLine, h3.mainLine, h4.mainLine],tract_name,'Original','Pred Full')
end


%plot([Node_cross Node_cross],[0,0.8],'-k','DisplayName','crossing')

set(gca, 'tickdir','out', 'ticklen',[0.025 0.025], ...
    'box','off','XTick', [0 round(Nnodes)/2 Nnodes], 'FontSize', 12);
xlim(gca,[1 Nnodes]);
if strcmp(measure,'FA')
    ylim(gca,[0 0.8]);
    yticks([0 0.2 0.4 0.6 0.8]);
else
    sup_lim = 1.1*max(nanmean(FA_tract)); % set lim y axis as +10% of maximum value on profile for tract1
    ylim(gca,[0 sup_lim]);
    
end

title_str = tract_name;
newStr = strrep(title_str,'_','-');
title(newStr, 'FontSize', 14)
xlabel('Nodes Along Tract', 'FontSize', 14);
ylabel(strcat(measure,' Value'), 'FontSize', 14);
hold off;

end

