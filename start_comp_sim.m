%% loop over methods
METHODS  =  [10]; % [9 8 7 10 5 4 2] %[11 10 8 7 5 4 2]; %[5 7 8]; % [2 4 5 7 8 9 10]
for method = 1:length(METHODS) 
%% set initial variables
cd('C:/Users/Michele/Documents/Jeff Miller/BDAGL/BDAGLMZ/')
seed = 1;
global Scen; Scen = 136;
global sim;
global chooseIntervention;
chooseIntervention = METHODS(method); 
    % 0 = only sample observational data
    % 1 = fixed sequence defined in nextIntervention.m
    % 2 = random node intervention
    % 3 = adaptive selection, using pairwise descendant entropy criteria
    % 4 = adaptive selection, using pairwise direct child entropy criteria
    % 5 = parent set entropy as entropy criterion
    % 6 = parent set entropy adjusted via MAP as entropy criterion
    % 7 = descendant set entropy 
    % 8 = child set entropy
    % 9 = Leong (DP algorithm no MCMC, adaptive child criterion)
    % 10 = cpdag entropy (over partition created by cpdags)
    % 11 = bninfo (fixed sequence imported)
global maxExperiment; maxExperiment = 9;
global sampleSachs; sampleSachs = 2; % 0 = generate from CPD, 
                                     % 1 = sample from 5400 original Sachs, 
                                     % 2 = sample from simulated data
                                     % specified at top of runToySachs2.m
global allowRepeatIntv; allowRepeatIntv = 0;
%bnet = mkSachsBnet(); 
bnet = mkBnet('tree8'); % name of structure like 'line' or 'tree'
global nInitialObs; nInitialObs = 250;
global nObsCases; nObsCases = 0; % number of obs to sample after the 1st experiment
global nIntvCases; nIntvCases = 200;
NSIMS = 40;
randomCPD = 2; p = 0.8;% randomCPD = 0 means load preexising settings with path below,
                       % 1 = make a CPD with the specified p; 
                       % 2 = use the one specified in mkBnet.m or not necessary to specify since using existing data  
settingsPath = 'C:/Users/Michele/Documents/Jeff Miller/BDAGL/BDAGLMZ/Sim Scenarios/compStructures/Leong/Scen122/settings.mat';
global allowMoreObs; allowMoreObs = 0;
global MCMCsamples; MCMCsamples = 75000; %20000 
global MCMCburnin; MCMCburnin = 150000; % 5000
global MCMCglobalFrac; MCMCglobalFrac = 0.9;
global psCutoff; psCutoff = 0; % set this to zero if you don't want a cutoff
global psMAP; psMAP = 0; % set this to zero to use the empirical parent set probabilities w/ no adjustment
global scenPath;
global mkPlots; mkPlots = 0;
global mkHamming; mkHamming = 0;
%rng(Scen);% set seed 
global randNodeSeq; randNodeSeq = NaN(NSIMS, maxExperiment-1); 
s = RandStream('mlfg6331_64');
if sampleSachs == 1
    sachsNodes = [2 4 7 8 9 ]; % include 9 twice since there is twice as much data for 9 as the others    
    for k=1:NSIMS
        randNodeSeq(k,:) = datasample(s, sachsNodes, maxExperiment-1,'Replace',false);
    end
end


% make array of random numbers for random node interventions
%randNodeSeq = randi(bnet.nNodes, NSIMS, maxExperiment); % no restriction
%on eligible nodes
% randNodeSeq = randi(length(bnet.eligibleNodes), NSIMS, maxExperiment);
% randNodeSeq = bnet.eligibleNodes(randNodeSeq); % map the random numbers to the nodes eligible for intervention
if (sampleSachs == 0 || sampleSachs == 2)
    for k=1:NSIMS
         randNodeSeq(k,:) = datasample(s, bnet.eligibleNodes, maxExperiment-1,'Replace',false);
    end
end


rng(seed);

%% set directory to save files in 
if chooseIntervention == 0
    scenPath = 'C:/Users/Michele/Documents/Jeff Miller/BDAGL/BDAGLMZ/Sim Scenarios/compStructures/Obs/';
    cd(scenPath); 
end
if chooseIntervention == 1
    scenPath = 'C:/Users/Michele/Documents/Jeff Miller/BDAGL/BDAGLMZ/Sim Scenarios/compStructures/Fixed/';
    cd(scenPath); 
end
if chooseIntervention == 2
    scenPath = 'C:/Users/Michele/Documents/Jeff Miller/BDAGL/BDAGLMZ/Sim Scenarios/compStructures/Random/';
    cd(scenPath); 
end
if chooseIntervention == 3
    scenPath = 'C:/Users/Michele/Documents/Jeff Miller/BDAGL/BDAGLMZ/Sim Scenarios/compStructures/AdaptiveDesc/';
    cd(scenPath); 
end
if chooseIntervention == 4
    scenPath = 'C:/Users/Michele/Documents/Jeff Miller/BDAGL/BDAGLMZ/Sim Scenarios/compStructures/AdaptiveChild/';
    cd(scenPath); 
end
if chooseIntervention == 5
    scenPath = 'C:/Users/Michele/Documents/Jeff Miller/BDAGL/BDAGLMZ/Sim Scenarios/compStructures/ParentSet/';
    cd(scenPath); 
end
if chooseIntervention == 6
    scenPath = 'C:/Users/Michele/Documents/Jeff Miller/BDAGL/BDAGLMZ/Sim Scenarios/compStructures/ParentSetMAP/';
    cd(scenPath); 
end
if chooseIntervention == 7
    scenPath = 'C:/Users/Michele/Documents/Jeff Miller/BDAGL/BDAGLMZ/Sim Scenarios/compStructures/DescSet/';
    cd(scenPath); 
end
if chooseIntervention == 8
    scenPath = 'C:/Users/Michele/Documents/Jeff Miller/BDAGL/BDAGLMZ/Sim Scenarios/compStructures/ChildSet/';
    cd(scenPath); 
end
if chooseIntervention == 9
    scenPath = 'C:/Users/Michele/Documents/Jeff Miller/BDAGL/BDAGLMZ/Sim Scenarios/compStructures/Leong/';
    cd(scenPath); 
end
if chooseIntervention == 10
    scenPath = 'C:/Users/Michele/Documents/Jeff Miller/BDAGL/BDAGLMZ/Sim Scenarios/compStructures/Cpdag/';
    cd(scenPath); 
end
if chooseIntervention == 11
    scenPath = 'C:/Users/Michele/Documents/Jeff Miller/BDAGL/BDAGLMZ/Sim Scenarios/compStructures/Bninfo/';
    cd(scenPath); 
    %sachsSeq = [4 8 2 9 7];
    %randNodeSeq = repmat([sachsSeq], NSIMS, 1)
    randNodeSeq = csvread('C:/Users/Michele/Documents/Jeff Miller/BDAGL/BDAGLMZ/Sim Data/EightTree_200/intvSeq/all_intvSeq.csv');
end
scenFolder = sprintf('Scen%d',Scen);
    if ~exist(scenFolder, 'dir') % if it doesn't exist make the folder
       mkdir(scenFolder)
    else 
        Scen = Scen + 1;
        scenFolder = sprintf('S%d',Scen)
        mkdir(scenFolder)
    end

scenPath = sprintf('%s%s/',scenPath, scenFolder);
%%  I temporarily used this when I ran all 120 permutations of the 5 node intvs for the Sachs data
%randNodeSeq = csvread('C:/Users/Michele/Documents/Jeff Miller/BDAGL/BDAGLMZ/Sim Scenarios/compStructures/Random/Scen83/sachsPerm.csv');

%% make CPD\
if randomCPD == 0
    load(settingsPath, 'bnet'); % make sure CPT is consistent by importing it
elseif randomCPD == 1
    % make random CPD
    % from Bayes Net Toolbox help doc: To control the degree of randomness (entropy), 
    % you can sample each row of the CPT from a Dirichlet(p,p,...) distribution. 
    % If p << 1, this encourages "deterministic" CPTs (one entry near 1, the rest near 0).
    % If p = 1, each entry is drawn from U[0,1]. If p >> 1, the entries will all be near 1/k, 
    % where k is the arity of this node, i.e., each row will be nearly uniform.
    for i=1:bnet.nNodes
        k = bnet.ns(i);
        ps = parents(bnet.dag, i);
        psz = prod(bnet.ns(ps));
        CPT = dirichlet_sample(p*ones(1,k), psz); % this is in the Foreign directory (looks like it was previously called sample_dirichlet)
        bnet.CPD{i} = tabular_CPD(bnet, i, 'CPT', CPT);      
        %bnet.CPD{i} = tabular_CPD(bnet, i, 'CPT', CPT, 'dirichlet_type', 'BDeu' );
    end
% uncomment below if you want to inspect the CPTs for a given node
%     for i=1:bnet.nNodes
%         disp(i)
%         values = struct2cell(bnet.CPD{i}); values{1}
%     end
elseif randomCPD == 2
   
%% Generate a tabular CPD(multinomial model) on each node
%labels = {'asia','tub','smoke','lung','either','bronc','dysp','xray'};
%asia = 1; tub = 2; smoke = 3; lung = 4; either = 5; bronc = 6;
%dysp = 7; xray = 8;


% cptA = [0.3 0.7];
%  cptT = [0.8 0.2; 0.3 0.7];
%  cptS = [0.4 0.6]; 
%  cptL = [0.9 0.1; 0.1 0.9];
%   cptE = zeros(2, 2, 2);
%  cptE(:,:,1) = [ .15 .3 ; .7 .9];
%  cptE(:,:,2) = 1-cptE(:,:,1);
%  cptB = [0.9 0.1; 0.1 0.9];
%  cptD = zeros(2, 2, 2);
%  cptD(:,:,1) = [ .15 .3 ; .7 .9];
%  cptD(:,:,2) = 1-cptD(:,:,1); 
%  cptX = [0.9 0.1; 0.1 0.9];
%  bnet.CPD{1} = tabular_CPD(bnet, 1, 'CPT', cptA);
%  bnet.CPD{2} = tabular_CPD(bnet, 2, 'CPT', cptT);
%  bnet.CPD{3} = tabular_CPD(bnet, 3, 'CPT', cptS);
%  bnet.CPD{4} = tabular_CPD(bnet, 4, 'CPT', cptL);
%  bnet.CPD{5} = tabular_CPD(bnet, 5, 'CPT', cptE);
%  bnet.CPD{6} = tabular_CPD(bnet, 6, 'CPT', cptB);
%  bnet.CPD{7} = tabular_CPD(bnet, 7, 'CPT', cptD);
%  bnet.CPD{8} = tabular_CPD(bnet, 8, 'CPT', cptX);
%     p = .9;
%     cptA = [.5 .5]; % P(A=1) = 0.4, P(A=2) = 0.6
%     cptB = [1-p p ; p 1-p]; %  I think value of parent is the column,so P(B=1|A=1) = 0.1, P(B=1|A=2)=0.9, P(B=2|A=1) = 0.9 
%     cptC = zeros(2, 2, 2);
%     cptC(:,:,1) = [ .15 .3 ; .7 .9];
%     cptC(:,:,2) = 1-cptC(:,:,1);
%     % cptC = [1-p p ; p 1-p];
%     % cptD = [1-p p ; p 1-p];
%     % cptE = [1-p p ; p 1-p];
%     bnet.CPD{1} = tabular_CPD(bnet, 1, 'CPT', cptA);
%     bnet.CPD{2} = tabular_CPD(bnet, 2, 'CPT', cptB);
%     bnet.CPD{3} = tabular_CPD(bnet, 3, 'CPT', cptC);
%     % bnet.CPD{4} = tabular_CPD(bnet, 4, 'CPT', cptD);
%     % bnet.CPD{5} = tabular_CPD(bnet, 5, 'CPT', cptE);
end    
%%
global trueDag;
%uncomment below if you want to make the network images
  
    figure;
    myDrawGraph(bnet.dag, 'labels', bnet.names);
    % title('Cancer Network');
    saveas(gcf, sprintf('%s/network.png',scenPath));
    
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
%variables = who;
%save(sprintf('%s/settings.txt',scenPath), 'NSIMS', '-ascii', '-tabs'); 
save(sprintf('%s/settings.mat',scenPath)); 
%% create objects for storing the intervention seq results
% and match true DAG results of the NSIMS simulations
intvSeqRes = NaN(NSIMS, maxExperiment); 
tprRes = NaN(NSIMS, maxExperiment); 
fprRes = NaN(NSIMS, maxExperiment); 
tnrRes = NaN(NSIMS, maxExperiment); 
fnrRes = NaN(NSIMS, maxExperiment); 
postEntropyRes = NaN(NSIMS, maxExperiment);
%skelEntropyRes = NaN(NSIMS, maxExperiment);
maxHRes = NaN(NSIMS, maxExperiment);
hammingMeanRes = NaN(NSIMS, maxExperiment);
hammingVarRes = NaN(NSIMS, maxExperiment);

%% Loop for running the simulations
for sim = 1:NSIMS
    sprintf('\n Starting simulation %d of %d. \n', sim, NSIMS)
    cd('C:\Users\Michele\Documents\Jeff Miller\BDAGL\BDAGLMZ')
    [interventionSeq, diagnostics, postE, maxH, hammingMean, hammingVar] = runToySachs2(sim, bnet);
    intvSeqRes(sim,1:maxExperiment) = interventionSeq;
    tprRes(sim,1:maxExperiment) = diagnostics.tpr;
    fprRes(sim,1:maxExperiment) = diagnostics.fpr;
    tnrRes(sim,1:maxExperiment) = diagnostics.tnr;
    fnrRes(sim,1:maxExperiment) = diagnostics.fnr;
    postEntropyRes(sim,1:maxExperiment) = postE;
    %skelEntropyRes(sim,:) = skelE;
    maxHRes(sim,1:maxExperiment) = maxH;
    hammingMeanRes(sim,1:maxExperiment) = hammingMean;
    hammingVarRes(sim,1:maxExperiment) = hammingVar;
    csvwrite(sprintf('%s/intvSeq_Results.csv',scenPath),intvSeqRes);
    csvwrite(sprintf('%s/tpr_Results.csv',scenPath),tprRes);
    csvwrite(sprintf('%s/fpr_Results.csv',scenPath),fprRes);
    csvwrite(sprintf('%s/tnr_Results.csv',scenPath),tnrRes);
    csvwrite(sprintf('%s/fnr_Results.csv',scenPath),fnrRes);
    csvwrite(sprintf('%s/postEntropy_Results.csv',scenPath),postEntropyRes);
    %csvwrite(sprintf('%s/skelEntropy_Results.csv',scenPath),skelEntropyRes);
    csvwrite(sprintf('%s/maxH_Results.csv',scenPath),maxHRes);
    csvwrite(sprintf('%s/randNodeSeq.csv',scenPath),randNodeSeq);
    csvwrite(sprintf('%s/hammingMean_Results.csv',scenPath),hammingMeanRes);
    csvwrite(sprintf('%s/hammingVar_Results.csv',scenPath),hammingVarRes);
end % end for loop over SIMS
end % end for loop over METHODS
% csvwrite(sprintf('%s/intvSeq_Results.csv',scenPath),intvSeqRes);
% csvwrite(sprintf('%s/matchTD_Results.csv',scenPath),matchRes);
% csvwrite(sprintf('%s/postEntropy_Results.csv',scenPath),postEntropyRes);
% csvwrite(sprintf('%s/randNodeSeq.csv',scenPath),randNodeSeq);
% csvwrite(sprintf('%s/hammingMean_Results.csv',scenPath),hammingMeanRes);
% csvwrite(sprintf('%s/hammingVar_Results.csv',scenPath),hammingVarRes);
    
%%
%Rdag1 = csvread("C:/Users/Michele/Documents/Jeff Miller/BDAGL/BDAGLMZ/Sim Scenarios/Sachs/Random/Scen11/Sim1/PCmat_Exp7.csv");
    figure;
    imagesc(bnet.dag, [0 1 ]);
    title('Benchmark DAG','Interpreter','latex');
    %set(gca, 'XTick', 1:bnet.nNodes);
    %set(gca, 'YTick', 1:bnet.nNodes);
    set(gca, 'XTickLabel', labels); xtickangle(45);
    set(gca, 'YTickLabel', labels);
    colorbar;
    saveas(gcf, sprintf('%s/BenchmarkDAG.png',scenPath));
