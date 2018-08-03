function [] = Gen_tract_profiles_pair(info, tract_name1, tract_name2, std_parameter, measure)
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
    tract1 = Get_tract_number(tract_name1);
    tract2 = Get_tract_number(tract_name2);
else
    tract1 = Get_tract_number_Dan(tract_name1);
    tract2 = Get_tract_number_Dan(tract_name2);
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
fgName1 = fe.life.M.tracts{tract1}.name; fgName1 = strrep(fgName1,' ','');

ind_tract1 = fe.life.M.tracts{tract1}.ind;
ind_tract_nnz1 = intersect(ind_tract1,ind_nnz);

fgTract1 = fe.fg.fibers(ind_tract_nnz1);
fgex1 = MyfgCreate_img('name', fgName1, 'colorRgb', [1 0 0], 'fibers', fgTract1);
% clean fibers with mba
fgcx1 = mbaComputeFibersOutliers(fgex1, std_parameter, std_parameter, 100, 'mean');


%% Extract fibers tract2
fgName2 = fe.life.M.tracts{tract2}.name; fgName2 = strrep(fgName2,' ','');

ind_tract2 = fe.life.M.tracts{tract2}.ind;
ind_tract_nnz2 = intersect(ind_tract2,ind_nnz);

fgTract2 = fe.fg.fibers(ind_tract_nnz2);
fgex2 = MyfgCreate_img('name', fgName2, 'colorRgb', [0 0 1], 'fibers', fgTract2);
% clean fibers with mba
fgcx2 = mbaComputeFibersOutliers(fgex2, std_parameter, std_parameter, 100, 'mean');

%% Compute profile tract1 using measure (FA,MD,etc) based on  tract 1 ONLY
file1 = deblank(ls(char(fullfile(dataPath,strcat(measure,'s'),strcat(strcat(measure,'_'),tract_name1,'.nii.gz')))));
famp1 = niftiRead(file1);
[Meas_tract1, SuperFiber1, ~, ~] = Compute_FA_AlongFG(fgcx1, famp1, [], [], Nnodes);

%% Compute profile tract1 using measure (FA,MD,etc) based on tract1 + tract2 + other tracts
FileName = strcat(tract_name1, '_',tract_name2);
file12 = deblank(ls(char(fullfile(dataPath,strcat(measure,'s'),strcat(strcat(measure,'_'),FileName,'_new.nii.gz')))));
famp12 = niftiRead(file12);
[Meas_tract1_12, ~]= Compute_FA_AlongFG(fgcx1, famp12, [], [], Nnodes);

%% Compute tract profile using measure (FA,MD,etc)FA based on original
fileOrig = deblank(ls(char(fullfile(dataPath,strcat(measure,'s'),strcat(measure,'_original.nii.gz')))));
fampOrig = niftiRead(fileOrig);
[Meas_tract1_orig, ~]= Compute_FA_AlongFG(fgcx1, fampOrig, [], [], Nnodes);
[Meas_tract2_orig, ~]= Compute_FA_AlongFG(fgcx2, fampOrig, [], [], Nnodes);

%% Compute profile tract2 using measure (FA,MD,etc) based on  tract 2 ONLY
file2 = deblank(ls(char(fullfile(dataPath,strcat(measure,'s'),strcat(strcat(measure,'_'),tract_name2,'.nii.gz')))));
famp2 = niftiRead(file2);
[Meas_tract2, SuperFiber2, ~, ~] = Compute_FA_AlongFG(fgcx2, famp2, [], [], Nnodes);

%% Compute profile tract2 using measure (FA,MD,etc) based on tract1 + tract2 + other tracts
[Meas_tract2_12, ~]= Compute_FA_AlongFG(fgcx2, famp12, [], [], Nnodes);

%% Compute tract profile using measure (FA,MD,etc) based on prediction
filePred = deblank(ls(char(fullfile(dataPath,strcat(measure,'s'),strcat(measure,'_pred_full.nii.gz')))));
fampPred = niftiRead(filePred);
[Meas_tract1_pred, ~]= Compute_FA_AlongFG(fgcx1, fampPred, [], [], Nnodes);
[Meas_tract2_pred, ~]= Compute_FA_AlongFG(fgcx2, fampPred, [], [], Nnodes);

%% Find tract crossing point
A1 = reshape(SuperFiber1.fibers{1},[3,Nnodes,1]);   % 3 x N     -> 3 x N x 1
A1 = repmat(A1,1,1,Nnodes);                         % 3 x N x 1 -> 3 x N x N

A2 = reshape(SuperFiber2.fibers{1},[3,1,Nnodes]);   % 3 x 1 x N -> 3 x 1 x N
A2 = repmat(A2,1,Nnodes,1);                         % 3 x 1 x N -> 3 x N x N

dist = squeeze(sum((A1 - A2).^2,1));
[Y,indrow] = min(dist);
[val,indcol] = min(Y);

Node_cross1 = indrow(indcol);
Node_cross2 = indcol;

%% Plot tract1 profile
Gen_profile_plot_new(Meas_tract1,'r',Meas_tract1_12,'g',Meas_tract1_orig,'k', Meas_tract1_pred,'y',tract_name1, tract_name2, 10, Nnodes, Node_cross1, measure)
saveas(gcf, strcat('./results/figures/',measure,'/',measure,'_profile_',tract_name1,'.fig'));
saveas(gcf, strcat('./results/figures/',measure,'/',measure,'_profile_',tract_name1,'.pdf'));
saveas(gcf, strcat('./results/figures/',measure,'/',measure,'_profile_',tract_name1,'.png'));

%% Plot tract2 profile
Gen_profile_plot_new(Meas_tract2,'b',Meas_tract2_12,'g', Meas_tract2_orig,'k',  Meas_tract2_pred,'y',tract_name2, tract_name1, 10, Nnodes, Node_cross2, measure)
saveas(gcf, strcat('./results/figures/',measure,'/',measure,'_profile_',tract_name2,'.fig'));
saveas(gcf, strcat('./results/figures/',measure,'/',measure,'_profile_',tract_name2,'.pdf'));
saveas(gcf, strcat('./results/figures/',measure,'/',measure,'_profile_',tract_name2,'.png'));

profiles_data.tract_name1 = tract_name1;
profiles_data.tract_name2 = tract_name2;
profiles_data.tract1 = Meas_tract1;
profiles_data.tract2 = Meas_tract2;
profiles_data.tract1_12 = Meas_tract1_12;
profiles_data.tract2_12 = Meas_tract2_12;
profiles_data.tract1_orig = Meas_tract1_orig;
profiles_data.tract2_orig = Meas_tract2_orig;
profiles_data.tract1_pred = Meas_tract1_pred;
profiles_data.tract2_pred = Meas_tract2_pred;

save(strcat('./results/',measure,strcat('_',tract_name1,'_',tract_name2),'.mat'), 'profiles_data')
end

function [] = Gen_profile_plot(FA_tract1, clr1, FA_tract1_12, clr12, FA_tract1_orig, clrorig, tract_name1, tract_name2, s, Nnodes, Node_cross)

N = size(FA_tract1_12,1);
figure
hold on
if ~isempty(FA_tract1_12)
    h2 = shadedErrorBar(1:Nnodes,nanmean(FA_tract1_12,1),s*nanstd(FA_tract1_12)/sqrt(N),'lineprops',clr12);
end
if ~isempty(FA_tract1)
    h1 = shadedErrorBar(1:Nnodes,nanmean(FA_tract1,1),s*nanstd(FA_tract1)/sqrt(N),'lineprops',clr1);
end
if ~isempty(FA_tract1_orig)
    h3 = shadedErrorBar(1:Nnodes,nanmean(FA_tract1_orig,1),s*nanstd(FA_tract1_orig)/sqrt(N),'lineprops',clrorig);
end

%
%legend([h1.mainLine, h2.mainLine],tract_name1,strcat(tract_name1,'+',tract_name2))
if ~isempty(FA_tract1)&&isempty(FA_tract1_orig)
   legend([h1.mainLine, h2.mainLine],tract_name1,strcat(tract_name1,'+',tract_name2))
elseif isempty(FA_tract1)
   legend([h2.mainLine, h3.mainLine],'Pred full','Orig')
else
   legend([h1.mainLine, h2.mainLine, h3.mainLine],tract_name1,strcat(tract_name1,'+',tract_name2),'Original')
end


%plot([Node_cross Node_cross],[0,0.8],'-k','DisplayName','crossing')

set(gca, 'tickdir','out', 'ticklen',[0.025 0.025], ...
    'box','off','XTick', [0 round(Nnodes)/2 Nnodes], 'YTick', [0 0.2 0.4 0.6 0.8], 'FontSize', 12);
xlim(gca,[1 Nnodes]);
ylim(gca,[0 0.8]);

title_str = tract_name1;
newStr = strrep(title_str,'_','-');
title(newStr, 'FontSize', 14)
xlabel('Nodes Along Tract', 'FontSize', 14);
ylabel(strcat(measure,' Value'), 'FontSize', 14);
hold off;

end

function [] = Gen_profile_plot_new(FA_tract1, clr1, FA_tract1_12, clr12, FA_tract1_orig, clrorig, FA_pred, clrp,tract_name1, tract_name2, s, Nnodes, Node_cross, measure)

N = size(FA_tract1_12,1);
figure
hold on
if ~isempty(FA_tract1_12)
    h2 = shadedErrorBar(1:Nnodes,nanmean(FA_tract1_12,1),s*nanstd(FA_tract1_12)/sqrt(N),'lineprops',clr12);
end
if ~isempty(FA_tract1)
    h1 = shadedErrorBar(1:Nnodes,nanmean(FA_tract1,1),s*nanstd(FA_tract1)/sqrt(N),'lineprops',clr1);
end
if ~isempty(FA_tract1_orig)
    h3 = shadedErrorBar(1:Nnodes,nanmean(FA_tract1_orig,1),s*nanstd(FA_tract1_orig)/sqrt(N),'lineprops',clrorig);
end
if ~isempty(FA_pred)
    h4 = shadedErrorBar(1:Nnodes,nanmean(FA_pred,1),s*nanstd(FA_pred)/sqrt(N),'lineprops',clrp);
end

%
%legend([h1.mainLine, h2.mainLine],tract_name1,strcat(tract_name1,'+',tract_name2))
if ~isempty(FA_tract1)&&isempty(FA_tract1_orig)
   legend([h1.mainLine, h2.mainLine],tract_name1,strcat(tract_name1,'+',tract_name2))
elseif isempty(FA_tract1)
   legend([h2.mainLine, h3.mainLine],'Pred full','Orig')
else
   legend([h1.mainLine, h2.mainLine, h3.mainLine, h4.mainLine],tract_name1,strcat(tract_name1,'+',tract_name2),'Original','Pred Full')
end


%plot([Node_cross Node_cross],[0,0.8],'-k','DisplayName','crossing')

set(gca, 'tickdir','out', 'ticklen',[0.025 0.025], ...
    'box','off','XTick', [0 round(Nnodes)/2 Nnodes], 'FontSize', 12);
xlim(gca,[1 Nnodes]);
if strcmp(measure,'FA')
    ylim(gca,[0 0.8]);
    yticks([0 0.2 0.4 0.6 0.8]);
else
    sup_lim = 1.1*max(nanmean(FA_tract1)); % set lim y axis as +10% of maximum value on profile for tract1
    ylim(gca,[0 sup_lim]);
end

title_str = tract_name1;
newStr = strrep(title_str,'_','-');
title(newStr, 'FontSize', 14)
xlabel('Nodes Along Tract', 'FontSize', 14);
ylabel(strcat(measure,' Value'), 'FontSize', 14);
hold off;

end

