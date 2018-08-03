function [FA_profile, SuperFiber, fgClipped, fgResampled] = ...
    Compute_FA_AlongFG(fg, dt, roi1, roi2, numberOfNodes)
%   Compute a weighted average of a variable (FA) in a track segment
%
%  [fa, SuperFiber,fgClipped, fgResampled] = ...
%    dtiComputeDiffusionPropertiesAlongFG(fg, dt, roi1, roi2, numberOfNodes, [dFallOff])
%
%   From a fiber group (fg), and diffusion data (dt), compute the average FA 
%   along the fiber track segment between the ROIs, sampled at numberOfNodes point.
%
% INPUTS:
%       fg            - fiber group structure.
%       dt            - dt6.mat structure or a nifti image.  If a nifti
%                       image is passed in then only 1 value will be 
%                       output and others will be empty
%       roi1          - first ROI for the fg
%       roi2          - second ROI for the fg
%       numberOfNodes - number of samples taken along each fg
%
% OUTPUTS:
%       fa         - Weighted fractional anisotropy
%       SuperFiber - structure containing the core of the fiber group
%       fgClipped  - fiber group clipped to the two ROIs
%       fgResampled- The fiber group that has been resampled to 
%                    numberOfNodes and each fiber has been reoriented to 
%                    start and end in a consitent location  
%
% Adapted from function dtiComputeDiffusionPropertiesAlongFG() in VISTASOFT
% AUTHOR Cesar Caiafa, Feb 2017. 
%
% Copyright 
%   Cesar Caiafa, Soichi Hayashi and Franco Pestilli
% 
% Indiana University 2018
% brainlife.io
%
% CC-BY 3.0 License CREDIT MUST BE GIVEN FOR ALL REUSE.

% If two rois are passed in clip the fiber group to the portion that spans
% between the ROIs
if ~notDefined('roi1') && ~notDefined('roi2')
    fgClipped = dtiClipFiberGroupToROIs(fg,roi1,roi2);
    % compute average FA along clipped fiber tract
    [FA_profile , SuperFiber, fgResampled] = ...
        dtiFiberGroup_FA_Average(fg, dt, numberOfNodes);
else
    % compute average FA along full fiber tract
    [FA_profile, SuperFiber, fgResampled] = ...
        dtiFiberGroup_FA_Average(fg, dt, numberOfNodes);
    % There is no clipped fiber group
    fgClipped = nan;
end

return
