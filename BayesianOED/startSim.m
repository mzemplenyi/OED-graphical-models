%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% startSim.m: script that initiates the Bayesian OED algorithm
%             using the settings specified at the beginning of
%             this file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% point to the DAG that will be used to generate data
bnet = mkBnet('ten2'); % name of structure like 'line' or 'tree'
% set simulation settings (stored in 'stg' structure)
%   this object will be passed to other functions
stg.method  = 'MEC'; 
stg.seed = 1;
stg.scen = 1;
stg.maxExp = 6;
stg.nInitialObs = 1800;
stg.nIntvCases = 600;
% MCMC settings
stg.MCMCsamples = 100000; 
stg.MCMCburnin = 150000; 
stg.MCMCglobalFrac = 0.9;
global sim;
global sampleSachs; sampleSachs = 1; % 0 = generate from CPD, 
                                     % 1 = sample from 5400 original Sachs, 
                                     % 2 = sample from simulated data
                                     % specified at top of runToySachs2.m


NSIMS = 10;


global scenPath;


global randNodeSeq; randNodeSeq = NaN(NSIMS, stg.maxExp-1); 
s = RandStream('mlfg6331_64');
if sampleSachs == 1
    sachsNodes = [2 4 7 8 9 ]; % (had a note to include 9 twice since there is twice as much data for 9 as the others, but I address by sampling 1200 in sampleSachsData.m if intv = 9)   
    for k=1:NSIMS
        randNodeSeq(k,:) = datasample(s, sachsNodes, stg.maxExp-1,'Replace',false);
    end
end


% make array of random numbers for random node interventions
%randNodeSeq = randi(bnet.nNodes, NSIMS, stg.maxExp); % no restriction
%on eligible nodes
% randNodeSeq = randi(length(bnet.eligibleNodes), NSIMS, stg.maxExp);
% randNodeSeq = bnet.eligibleNodes(randNodeSeq); % map the random numbers to the nodes eligible for intervention
if (sampleSachs == 0 || sampleSachs == 2)
    for k=1:NSIMS
         randNodeSeq(k,:) = datasample(s, bnet.eligibleNodes, stg.maxExp-1,'Replace',false);
    end
end


rng(stg.seed);

%% set directory to save files in 

if strcmp(stg.method, 'Random')
    scenPath = 'C:/Users/Michele/Documents/Jeff Miller/BDAGL/BDAGLMZ/Sim Scenarios/compStructures/Random/';
    cd(scenPath); 
end
if strcmp(stg.method, 'DP')
    scenPath = 'C:/Users/Michele/Documents/Jeff Miller/BDAGL/BDAGLMZ/Sim Scenarios/compStructures/Leong/';
    cd(scenPath); 
end
if strcmp(stg.method, 'MEC')
    scenPath = 'C:/Users/Michele/Documents/Jeff Miller/BDAGL/BDAGLMZ/Sim Scenarios/compStructures/MEC/';
    cd(scenPath); 
end
% create folder to store results
scenFolder = sprintf('Scen%d',stg.scen);
    if ~exist(scenFolder, 'dir') % if it doesn't exist make the folder
       mkdir(scenFolder)
    else 
        stg.scen = stg.scen + 1; % increment scen number by 1 folder already exists
        scenFolder = sprintf('S%d',stg.scen);
        mkdir(scenFolder)
    end

scenPath = sprintf('%s%s/',scenPath, scenFolder);


    
    figure;
    imagesc(bnet.dag, [0 1 ]);
    title('Adjacency Matrix for True DAG','Interpreter','latex');
    set(gca, 'XTick', 1:bnet.nNodes);
    set(gca, 'YTick', 1:bnet.nNodes);
    set(gca, 'XTickLabel', bnet.names);
    set(gca, 'YTickLabel', bnet.names);
    set(gcf, 'Position', [0 0 500 400])
    colorbar;
    saveas(gcf, sprintf('%s/BenchmarkDAG.png',scenPath));

%% create settings file
save(sprintf('%s/settings.mat',scenPath)); 
%% create objects for storing the intervention seq results
% and match true DAG results of the NSIMS simulations
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
    cd('C:\Users\Michele\Documents\Jeff Miller\BDAGL\BDAGLMZ')
    [interventionSeq, diagnostics, postE, maxH, hammingMean, hammingVar] = runToySachs2(sim, bnet);
    intvSeqRes(sim,1:stg.maxExp) = interventionSeq;
    tprRes(sim,1:stg.maxExp) = diagnostics.tpr;
    fprRes(sim,1:stg.maxExp) = diagnostics.fpr;
    tnrRes(sim,1:stg.maxExp) = diagnostics.tnr;
    fnrRes(sim,1:stg.maxExp) = diagnostics.fnr;
    postEntropyRes(sim,1:stg.maxExp) = postE;
    maxHRes(sim,1:stg.maxExp) = maxH;
    hammingMeanRes(sim,1:stg.maxExp) = hammingMean;
    hammingVarRes(sim,1:stg.maxExp) = hammingVar;
    csvwrite(sprintf('%s/intvSeq_Results.csv',scenPath),intvSeqRes);
    csvwrite(sprintf('%s/tpr_Results.csv',scenPath),tprRes);
    csvwrite(sprintf('%s/fpr_Results.csv',scenPath),fprRes);
    csvwrite(sprintf('%s/tnr_Results.csv',scenPath),tnrRes);
    csvwrite(sprintf('%s/fnr_Results.csv',scenPath),fnrRes);
    csvwrite(sprintf('%s/postEntropy_Results.csv',scenPath),postEntropyRes);
    csvwrite(sprintf('%s/maxH_Results.csv',scenPath),maxHRes);
    csvwrite(sprintf('%s/randNodeSeq.csv',scenPath),randNodeSeq);
    csvwrite(sprintf('%s/hammingMean_Results.csv',scenPath),hammingMeanRes);
    csvwrite(sprintf('%s/hammingVar_Results.csv',scenPath),hammingVarRes);
end % end for loop over SIMS

