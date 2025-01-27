function [sampData, sampClamped, remData, remClamped, bnet ] = sampleSachsData(bnet, nObs, nIntv, interventionSeq, expNum, remData, remClamped)
% sampleSachsData.m: 
%   Sample data from either the observed or interventional data
%   provided with the dataset (since we do not have a CPD that describes
%   the Sachs network)
%  OUTPUT: 
%       --seqData: a matrix with rows being nodes and columns being data 
%        samples (either obs or interventional). This structure grows 
%       (more columns are added) with each experiment performed.
%       --seqClamped: a binary matrix of same dimension as seqData with 
%       indicators for whether or not a node was intervened on for an
%       individual data sample.
%
%       These outputs are used in runExperiment.m by
%       mkAllFamilyLogMargLik()


%% transpose matrix to make the row search more intuitive
remClampedT = remClamped';
remDataT = remData';
nSamples = nObs;
%% make the intervention (or observational) vector (e.g. 0 1 0 0 0 ... for intervening on node 2)
intvVec = zeros(1,size(remDataT,2));
if interventionSeq(expNum) ~= 0 % if not an observational experiment
    intvNode = interventionSeq(expNum); % set which node to intervene on
    if ~ismember(bnet.eligibleNodes, intvNode)
        error('The selected node cannot be perturbed.')
    end
    intvVec(intvNode) = 1;
    if (intvNode == 9) % PKC has 1200 intervention samples whereas others have 600 
        nSamples = 1200; % change to 1200 if sampling from real Sachs data
    else
        nSamples = nIntv;
    end
end
%% find rows in clamped matrix corresponding to desired intv vector
ix = find(ismember(remClampedT, intvVec,'rows')); % returns row indices of tclamped that match desired intv vec
keepIx = randsample(ix, nSamples);
sampDataT = remDataT(keepIx,:);
sampClampedT = remClampedT(keepIx,:);
%% If later version of code allows for repeated interventions consider uncommenting the below if statement. 
% Check whether there is nSamples worth of data remaining for that type
% of data to sample from. If not remove intvNode from eligibleNodes for 
% future interventions. Note that intvNode only exists if we are not
% on experiment 1 (obs data).  We use (2*nSamples) because length(ix) doesn't yet have the
% nSamples removed for this experiment.

% if (expNum > 1 && length(ix) < (2*nSamples))
%     bnet.eligibleNodes(bnet.eligibleNodes == intvNode) = [];
% end     
%% remove the sampled rows from the data and clamped matrices to update them 
remDataT(keepIx,:) = [];
remClampedT(keepIx,:) = [];
%% transpose back to nNodes x nSamples matrices
remData = remDataT';
remClamped = remClampedT';
sampData = sampDataT';
sampClamped = sampClampedT';
end

