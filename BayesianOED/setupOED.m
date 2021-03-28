%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setupOED.m: script that initiates the Bayesian optimal 
%             experimental design (OED) algorithm described
%             by Zemplenyi and Miller (2021) using the settings 
%             specified at the beginning of this file.
%             Results files to save are specified
%             at the end of this file. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% point to the DAG that will be used to generate data
stg.bnet = mkBnet('ten2'); % name of structure like 'line' or 'tree'
% set simulation settings (stored in 'stg' structure)
%   this object will be passed to other functions
NSIMS = 5;
stg.method  = 'MEC'; % specify partition type or other method, e.g. 'Random' 
stg.scen = 1;
stg.maxExp = 6;
stg.nInitialObs = 1800;
stg.nIntvCases = 600;
stg.seed = 1;

%% MCMC settings
stg.MCMCsamples = 100000; 
stg.MCMCburnin = 150000; 
stg.MCMCglobalFrac = 0.9;

%% setup directory to save results files in 
% first check the path of the working directory
pwd() % should end in '\OED-graphical-models\BayesianOED'

% add the existing 'Results' folder to the working path 
resPath = sprintf('%s/Results/',pwd());

% create a new subdirectory to store this scenario's results
cd(resPath)
dirName = sprintf('%s_Scen%d',stg.method, stg.scen);
    if ~exist(dirName, 'dir') % if directory doesn't exist make the folder
       mkdir(dirName)
       fprintf('Created a new directory named "%s".\n',dirName);
    else 
        % increment scenario number by 1 if directory already exists
        % to prevent overwriting results        
        stg.scen = stg.scen + 1; 
        dirName = sprintf('%s_Scen%d',stg.method, stg.scen);
        fprintf('Directory name already exists for that scenario number.\nCreated a new directory named "%s"\n',dirName); 
        mkdir(dirName)
    end
% now add this subdirectory to the results path
resPath = sprintf('%s%s/',resPath, dirName);

%% make array of nodes for random node interventions
%   sample from the eligible nodes without replacement
if isequal(stg.method, 'Random')
    randNodeSeq = NaN(NSIMS, stg.maxExp-1); 
    s = RandStream('mlfg6331_64');
    for k=1:NSIMS
         randNodeSeq(k,:) = datasample(s, bnet.eligibleNodes, stg.maxExp-1,'Replace',false);
    end
end

rng(stg.seed);
%% Make figure of benchmark / ground truth DAG used to generate data    
figure;
imagesc(bnet.dag, [0 1 ]);
title('Adjacency Matrix for True DAG','Interpreter','latex');
set(gca, 'XTick', 1:bnet.nNodes);
set(gca, 'YTick', 1:bnet.nNodes);
set(gca, 'XTickLabel', bnet.names);
set(gca, 'YTickLabel', bnet.names);
set(gcf, 'Position', [0 0 500 400])
colorbar;
saveas(gcf, sprintf('%s/benchmarkDAG.png',resPath));

%% create settings file
stg.bnet = bnet;
stg.resPath = resPath;
save(sprintf('%s/stg.mat',resPath)); 
%% initialize objects for storing the intervention seq and diagnostic results
intvSeqRes = NaN(NSIMS, stg.maxExp); 
tprRes = NaN(NSIMS, stg.maxExp); 
fprRes = NaN(NSIMS, stg.maxExp); 
tnrRes = NaN(NSIMS, stg.maxExp); 
fnrRes = NaN(NSIMS, stg.maxExp); 
postEntropyRes = NaN(NSIMS, stg.maxExp);
maxHRes = NaN(NSIMS, stg.maxExp);
hammingMeanRes = NaN(NSIMS, stg.maxExp);
hammingVarRes = NaN(NSIMS, stg.maxExp);

%% Loop for running the simulations
for sim = 1:NSIMS
    sprintf('\n Starting simulation %d of %d. \n', sim, NSIMS)
    [interventionSeq, diagnostics, postE, maxH, hammingMean, hammingVar] = startSim(stg, sim);
    intvSeqRes(sim,:) = interventionSeq;
    tprRes(sim,:) = diagnostics.tpr;
    fprRes(sim,:) = diagnostics.fpr;
    tnrRes(sim,:) = diagnostics.tnr;
    fnrRes(sim,:) = diagnostics.fnr;
    postEntropyRes(sim,:) = postE;
    maxHRes(sim,:) = maxH;
    hammingMeanRes(sim,:) = hammingMean;
    hammingVarRes(sim,:) = hammingVar;
    
    csvwrite(sprintf('%s/intvSeq_Results.csv',resPath),intvSeqRes);
    csvwrite(sprintf('%s/tpr_Results.csv',resPath),tprRes);
    csvwrite(sprintf('%s/fpr_Results.csv',resPath),fprRes);
    csvwrite(sprintf('%s/tnr_Results.csv',resPath),tnrRes);
    csvwrite(sprintf('%s/fnr_Results.csv',resPath),fnrRes);
    csvwrite(sprintf('%s/postEntropy_Results.csv',resPath),postEntropyRes);
    csvwrite(sprintf('%s/maxH_Results.csv',resPath),maxHRes);    
    csvwrite(sprintf('%s/hammingMean_Results.csv',resPath),hammingMeanRes);
    csvwrite(sprintf('%s/hammingVar_Results.csv',resPath),hammingVarRes);
    if isequal(stg.method, 'Random')
        csvwrite(sprintf('%s/randNodeSeq.csv',resPath),randNodeSeq);
    end
end % end for loop over SIMS

