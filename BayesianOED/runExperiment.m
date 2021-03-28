function [ bnet, entropy, dg, hamming, seqData, seqClamped, remData, remClamped] = runExperiment(stg, simPath, expNum, nObservationCases, interventions, nInterventionCases, seqData, seqClamped, remData, remClamped )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
bnet = stg.bnet;
entropy = []; % initialize object to store various types of entropy
hamming = [];
if isequal(stg.dataOrigin, 'simFromCPD')
    [sampData, sampClamped] = mkData(bnet, nObservationCases, interventions, nInterventionCases );
    % note that in this case "interventions" is a cell array indicating
    % which node to intervene on and what to set it to
end
if (isequal(stg.dataOrigin, 'useSachsData') || isequal(stg.dataOrigin, 'sampSimData'))
    [sampData, sampClamped, remData, remClamped, bnet ] = sampleSachsData(bnet, nObservationCases, nInterventionCases, interventions, expNum, remData, remClamped);
    % note that in this case "interventions" is a vector equal to
    % interventionSeq
end
seqData = [seqData sampData];
seqClamped = [seqClamped sampClamped];

%% Use DP algorithm to re-learn the structure
maxFanIn = bnet.nNodes - 1; % do not restrict the fan-in of any node

aflp = mkAllFamilyLogPrior( bnet.nNodes, 'maxFanIn', maxFanIn); % construct the modular prior on all possible local families; another option: use 'priorType', 'flat'
aflml = mkAllFamilyLogMargLik(seqData, 'nodeArity', repmat(bnet.node_sizes(1),1,bnet.nNodes), 'impossibleFamilyMask', aflp~=-Inf,'interventionType', 'perfect', 'clampedMask',seqClamped); % compute marginal likelihood on all possible local families
ep = computeAllEdgeProb( aflp, aflml ); % compute the marginal edge probabilities using

if isequal(stg.method, 'DP') 
    %Li and Leong method using only dynamic programming (no MCMC)
    PCmat = ep;
    H = entropyDP(bnet.nNodes, ep);
    dg = checkDiagnostics(PCmat, bnet.dag, simPath, expNum);
    hamming.mean = dg.hd;
    entropy.H = H;
    entropy.postEntropy = max(H);
else

    %% Use MCMC to correct the bias induced by modular prior
    [samples, ~, ~] = sampleDags(@uniformGraphLogPrior, aflml, bnet.dag, 'burnin', stg.MCMCburnin, ...
        'edgeMarginals', ep, 'globalFrac', stg.MCMCglobalFrac, 'thinning', 1, 'nSamples', stg.MCMCsamples);


    if isequal(stg.method, 'Observational') % obs data case
        H = childSetEntropy(bnet, samples); % doesn't matter what the entropy over children or desc is since we choose obs data anyway
        PCmat = samplesToEdgeMarginals(samples);
        dg = checkDiagnostics(PCmat, bnet.dag, simPath, expNum);
    elseif isequal(stg.method, 'Fixed') % fixed seq data case
        H = zeros(1, bnet.nNodes); % doesn't matter what the entropy over children or desc is since we choose randomly
        PCmat = samplesToEdgeMarginals(samples);
        dg = checkDiagnostics(PCmat, bnet.dag, simPath, expNum);    
    elseif isequal(stg.method, 'Random') % random case 
        H = zeros(1, bnet.nNodes); % doesn't matter what the entropy over children or desc is since we choose randomly
        PCmat = samplesToEdgeMarginals(samples);
        dg = checkDiagnostics(PCmat, bnet.dag, simPath, expNum);    
    elseif isequal(stg.method, 'PairWiseChild')
        H = entropyPW(bnet,samples);
        PCmat = samplesToEdgeMarginals(samples);
        dg = checkDiagnostics(PCmat, bnet.dag, simPath, expNum);
    elseif isequal(stg.method, 'ParentSet')
        H = parentSetEntropy(bnet, samples);
        PCmat = samplesToEdgeMarginals(samples);
        dg = checkDiagnostics(PCmat, bnet.dag, simPath, expNum);
    elseif isequal(stg.method, 'DescSet')
        H = descSetEntropy(bnet, samples);
        PCmat = samplesToEdgeMarginals(samples);
        dg = checkDiagnostics(PCmat, bnet.dag, simPath, expNum);
    elseif isequal(stg.method, 'ChildSet')
        H = childSetEntropy(bnet, samples);
        PCmat = samplesToEdgeMarginals(samples);
        dg = checkDiagnostics(PCmat, bnet.dag, simPath, expNum); 
    elseif isequal(stg.method, 'TSE')
       % H = cpdagEntropy(bnet, samples);
       H = tseEntropy(bnet, samples); % added  9/14/19 
       PCmat = samplesToEdgeMarginals(samples);
        dg = checkDiagnostics(PCmat, bnet.dag, simPath, expNum); 
    elseif isequal(stg.method, 'bninfo') % bninfo / fixed seq data case
        H = zeros(1, bnet.nNodes); % doesn't matter what the entropy over children or desc is since we choose according to the given file
        PCmat = samplesToEdgeMarginals(samples);
        dg = checkDiagnostics(PCmat, bnet.dag, simPath, expNum);  
    elseif isequal(stg.method, 'MEC')
       H = entropyMEC(bnet, samples); % added this on 2/15/21
       PCmat = samplesToEdgeMarginals(samples);
       dg = checkDiagnostics(PCmat, bnet.dag, simPath, expNum);     
    end

    hamming.mean = dg.hd;
    entropy.postEntropy = posteriorEntropy(samples);
    entropy.H = H;
end

