%%
load 'sachsDiscretizedData.mat';
% data is size d*N (nodes * cases), values are {1,2,3}


nNodes = 11;
maxFanIn = 10; 

%[dag, F, H] = sachsTrueDAG();
bnet = mkSachsBnet();
%%
nr = 2; nc = 2;
figure(1);clf
subplot(nr,nc,1)
%imagesc(dag, [0 1]);
imagesc(bnet.dag, [0 1]);
title('ground truth')
colorbar 
%%
% all families log prior (size d * 2^d)
% impossible configurations are -inf
%aflp = mkAllFamilyLogPrior( nNodes, 'maxFanIn', maxFanIn, 'priorType', 'flat' ); 
aflp = mkAllFamilyLogPrior( nNodes, 'maxFanIn', maxFanIn ); 
%%
% all families log  marginal likelihood (size d * 2^d)
aflml = mkAllFamilyLogMargLik( data, 'nodeArity', repmat(3,1,nNodes), ...
			       'impossibleFamilyMask', aflp~=-Inf, ...
			       'priorESS', 1, ...
			       'clampedMask', clamped);
            
%%
% Use DP to compute edge marginals using modular prior
epDP = computeAllEdgeProb( aflp, aflml ); 
%logZ = computeLogZ(aflp, aflml);
%% 
epDP0 = computeAllEdgeProb(0, aflml ); 
%%
ep = computeAllEdgeProb( aflp, aflml ); % Use DP to compute edge marginals using modular prior

%% Use MCMC to correct the bias induced by modular prior
mkPlots = 0;
[samples, diagnostics] = sampleDags(@uniformGraphLogPrior, aflml, bnet.dag, 'initialDag', bnet.dag, 'burnin', 15000, 'verbose', mkPlots, ...
    'edgeMarginals', ep, 'globalFrac', 0.9, 'thinning', 2, 'nSamples', 25000);

mcmcDag = samplesToEdgeMarginals(samples);
%%
figure(1);
subplot(nr,nc,2)
imagesc(mcmcDag, [0 1 ]);
title('edge marginals (DP)')
colorbar;
%% check if ground truth equals estimated dag from mcmc
isequal(bnet.dag, mcmcDag>0.5)
%%
figure(1);
subplot(nr,nc,3)
imagesc(epDP0, [0 1 ]);
title('edge marginals flat prior (DP)')
colorbar;

%%

% Use DP to find exact MAP structure
optimalDAG = computeOptimalDag(aflml); 
%optimalDAG = computeOptimalDag(aflml+aflp);,

figure(1);
subplot(nr,nc,4)
imagesc(optimalDAG, [0 1]);
title('optimal MAP')
colorbar
%%
% ROC curves
[AUC_dp, FPrate_dp, TPrate_dp, thresholds] = cROC(epDP(:), dag(:));
figure(2); clf
plot(FPrate_dp, TPrate_dp, 'r-');
legendstr{1} = sprintf('dp %3.2f', AUC_dp);
legend(legendstr)
%%

if 1
% Visualize graph structures
load('sachsScienceLayout')

figure(3); clf
myDrawGraph(dag, 'labels', labels, 'layout', graphicalLayout)
title('ground truth')

figure(4); clf
myDrawGraph(optimalDAG, 'labels', labels, 'layout', graphicalLayout)
title(sprintf('exact MAP'))
end
