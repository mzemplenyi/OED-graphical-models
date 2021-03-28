function [ bnet, entropy, dg, hamming, seqData, seqClamped, remData, remClamped] = runExperiment2( bnet, nObservationCases, interventions, nInterventionCases, seqData, seqClamped, remData, remClamped )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global simPath;
global Scen;
global experiment;
global mkPlots;
global mkHamming;
global chooseIntervention;
global allowMoreObs;
global MCMCsamples; 
global MCMCburnin;
global MCMCglobalFrac;
global sampleSachs;

entropy = []; % initialize object to store various types of entropy
hamming = [];
if sampleSachs == 0
    [sampData, sampClamped] = mkData(bnet, nObservationCases, interventions, nInterventionCases );
    % note that in this case "interventions" is a cell array indicating
    % which node to intervene on and what to set it to
end
if (sampleSachs == 1 || sampleSachs == 2)
    [sampData, sampClamped, remData, remClamped, bnet ] = sampleSachsData(bnet, nObservationCases, nInterventionCases, interventions, experiment, remData, remClamped);
    % note that in this case "interventions" is a vector equal to
    % interventionSeq
end
seqData = [seqData sampData];
seqClamped = [seqClamped sampClamped];

%% Use DP algorithm to re-learn the structure
maxFanIn = bnet.nNodes - 1; % do not restrict the fan-in of any node

% aflp = mkAllFamilyLogPrior( nNodes, 'maxFanIn', maxFanIn ); % construct the modular prior on all possible local families
% aflml = mkAllFamilyLogMargLik( data, 'nodeArity', repmat(2,1,nNodes), 'impossibleFamilyMask', aflp~=-Inf); % compute marginal likelihood on all possible local families
% ep = computeAllEdgeProb( aflp, aflml ); % compute the marginal edge probabilities using

aflp = mkAllFamilyLogPrior( bnet.nNodes, 'maxFanIn', maxFanIn); % construct the modular prior on all possible local families; another option: use 'priorType', 'flat'
aflml = mkAllFamilyLogMargLik(seqData, 'nodeArity', repmat(bnet.node_sizes(1),1,bnet.nNodes), 'impossibleFamilyMask', aflp~=-Inf,'interventionType', 'perfect', 'clampedMask',seqClamped); % compute marginal likelihood on all possible local families


ep = computeAllEdgeProb( aflp, aflml ); % compute the marginal edge probabilities using

if chooseIntervention == 9
    PCmat = ep;
    H = entropyDP(bnet.nNodes, ep);
    dg = checkDiagnostics(PCmat, bnet.dag);
    hamming.mean = dg.hd;
    hamming.var = 0;
    entropy.H = H;
    entropy.postEntropy = max(H);
else

%% Use MCMC to correct the bias induced by modular prior
[samples, diagnostics] = sampleDags(@uniformGraphLogPrior, aflml, bnet.dag, 'burnin', MCMCburnin, 'verbose', mkPlots, ...
    'edgeMarginals', ep, 'globalFrac', MCMCglobalFrac, 'thinning', 1, 'nSamples', MCMCsamples);
% [samples, diagnostics] = sampleDags(aflp, aflml, bnet.dag, 'burnin', MCMCburnin, 'verbose', mkPlots, ...
%     'edgeMarginals', ep, 'globalFrac', MCMCglobalFrac, 'thinning', 1, 'nSamples', MCMCsamples);

%% Explore Hamming Distance
if(mkHamming == 1)
    figure('visible','off');
    plot(samples.hammingDist)
    ylabel('Hamming Distance');
    xlabel('Iteration no.');
    saveas(gcf, sprintf('%s/Exp%d_hamming.png',simPath,experiment))
end
%hamming.mean = mean(samples.hammingDist);
%hamming.var = var(samples.hammingDist);

% if allowMoreObs == 1 
%     skelEntropy = skeletonEntropy(bnet, samples);
% else
%     skelEntropy = 0;
% end

% calculate entropy over skeletons no matter how we are choosing
% interventions
%skelEntropy = skeletonEntropy(bnet, samples);

if chooseIntervention == 0 % obs data case
    H = childSetEntropy(bnet, samples); % doesn't matter what the entropy over children or desc is since we choose obs data anyway
    PCmat = samplesToEdgeMarginals(samples);
    % check true positive edge rate
    dg = checkDiagnostics(PCmat, bnet.dag);
elseif chooseIntervention == 1 % fixed seq data case
    H = zeros(1, bnet.nNodes); % doesn't matter what the entropy over children or desc is since we choose randomly
    PCmat = samplesToEdgeMarginals(samples);
    % check true positive edge rate
    dg = checkDiagnostics(PCmat, bnet.dag);    
elseif chooseIntervention == 2 % random case 
    H = zeros(1, bnet.nNodes); % doesn't matter what the entropy over children or desc is since we choose randomly
    PCmat = samplesToEdgeMarginals(samples);
    % check true positive edge rate
    dg = checkDiagnostics(PCmat, bnet.dag);    
elseif chooseIntervention == 3 || chooseIntervention == 4
    H = entropyPW(bnet,samples);
    %PCmat = calcPCmat(bnet, samples);
    PCmat = samplesToEdgeMarginals(samples);
    dg = checkDiagnostics(PCmat, bnet.dag);
elseif chooseIntervention == 5 || chooseIntervention == 6
    H = parentSetEntropy(bnet, samples);
    PCmat = samplesToEdgeMarginals(samples);
    % check true positive edge rate
    dg = checkDiagnostics(PCmat, bnet.dag);
elseif chooseIntervention == 7
    H = descSetEntropy(bnet, samples);
    PCmat = samplesToEdgeMarginals(samples);
    % check true positive edge rate
    dg = checkDiagnostics(PCmat, bnet.dag);
elseif chooseIntervention == 8
    H = childSetEntropy(bnet, samples);
    PCmat = samplesToEdgeMarginals(samples);
    % check true positive edge rate
    dg = checkDiagnostics(PCmat, bnet.dag); 
elseif chooseIntervention == 10
   % H = cpdagEntropy(bnet, samples);
   H = tseEntropy(bnet, samples); % added this on 9/14/19 
   PCmat = samplesToEdgeMarginals(samples);
    % check true positive edge rate
    dg = checkDiagnostics(PCmat, bnet.dag); 
elseif chooseIntervention == 11 % bninfo / fixed seq data case
    H = zeros(1, bnet.nNodes); % doesn't matter what the entropy over children or desc is since we choose according to the given file
    PCmat = samplesToEdgeMarginals(samples);
    % check true positive edge rate
    dg = checkDiagnostics(PCmat, bnet.dag);   
end
hamming.mean = dg.hd;
hamming.var = 0;
postEntropy = posteriorEntropy(samples);
entropy.H = H;
entropy.postEntropy = postEntropy;
%entropy.skelEntropy = skelEntropy;
 
end

