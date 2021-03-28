function [interventionSeq, diagnostics, postE, maxH, hammingMean, hammingVar] = startSim(stg, sim )
% set stop threshold to any neg number if you only want to use maxExp as
% to stop running experiments. Otherwise set this to a small positive number.
stopThreshold = -1; 
intvValue = 1; % value to set node to for an intervention: 1 is low-level, 2 is high-level

simPath = sprintf('%sSim%d',stg.resPath, sim);
mkdir(simPath)
bnet = stg.bnet;

global experiment; experiment = 1;
global nObsCases; 
global nIntvCases;
global sampleSachs;

%% setup data to be sampled from for this simulation
if isequal(stg.dataOrigin, 'sampSimData')% sample from existing simulated data
    % read in data file corresponding to this simulation
    %cd(stg.dataDir)
    dataFile = sprintf('%sdata_%d.csv', stg.dataDir, sim);
    data = csvread(dataFile);
    data = data';
    % read in corresponding 'clamped' file
    clampedFile = sprintf('%sclamped_%d.csv', stg.dataDir, sim);
    clamped = csvread(clampedFile);
    clamped = clamped';
    remData = data; % remaining data to be sampled from (intially all data)
    remClamped = clamped;   
    
elseif isequal(stg.dataOrigin, 'useSachsData')
    %loads 'data' (11x5400 obs) and 'clamped'(11x5400 obs) 
    load(sprintf('%s/sachsDiscretizedData.mat', stg.dataDir); %% move this to runSachs.m file
    % initially all the data is the "remaining" data
    remData = data;
    remClamped = clamped;
    
elseif isequal(stg.dataOrigin, 'simFromCPD')
    remData = [];
    remClamped = [];
end

%% initialize matrices that will store the sequentially
% accumulated data over experiments and the corresponding clamped records
seqData = NaN(bnet.nNodes, 0);
seqClamped = NaN(bnet.nNodes, 0);

%% initialize intervention sequences and diagnostics
%cd 'C:/Users/Michele/Documents/GitHub/OED-graphical-models/BayesianOED/';
interventionSeq = NaN(1,stg.maxExp); % vector to store which experiments are run
interventionSeq(1) = 0; % first experiment is just observational 
tpr = NaN(1,stg.maxExp);
fpr = NaN(1,stg.maxExp);
tnr = NaN(1,stg.maxExp);
fnr = NaN(1,stg.maxExp);
postE = NaN(1,stg.maxExp);
maxH = NaN(1,stg.maxExp);
HExp = NaN(bnet.nNodes, stg.maxExp); 
hammingMean = NaN(1,stg.maxExp);
hammingVar = NaN(1,stg.maxExp);

%% BEGIN WHILE LOOP FOR RUNNING EXPERIMENTS
while(experiment <= stg.maxExp)
    if experiment == 1
    %% Generate only observational data for first experiment
        nObservationCases = stg.nInitialObs; % # observational data cases
        nInterventionCases = 0; % no interventions
        if isequal(stg.dataOrigin, 'simFromCPD')
            interventions = {};
        end
        if (isequal(stg.dataOrigin, 'useSachsData') || isequal(stg.dataOrigin, 'sampSimData'))
            interventions = interventionSeq;
        end
    else 
        % check if entropy stop criterion is met 
         if entropy.postEntropy < stopThreshold
            sprintf('\n Posterior entropy of %d dropped below stop threshold. Not running experiment %d \n', entropy.postEntropy, experiment)
            break;
         end  
        % make interventions cell array
        fprintf('\n Beginning Experiment %d \n', experiment)
        % select node to intervene on based on method
        interveneNode = nextIntervention(HExp(:,experiment-1), bnet.eligibleNodes, stg.method, sim); 
        interventionSeq(experiment) = interveneNode;
        if interveneNode == 0
            nObservationCases = nObsCases; % # observational data cases
            nInterventionCases = 0; % no interventions
            if sampleSachs == 0 % added this logic on 6/20/19
                interventions = {};
            end
            if sampleSachs == 1
                interventions = interventionSeq;
            end
        else       
            nObservationCases = nObsCases; % assumes we only generate obs data on first experiment
            nInterventionCases = nIntvCases;
            %nInterventionCases = nIntvCases*(length(interventionSeq)-1);
            % build up interventions cell array
            if sampleSachs == 0
                interventions = cell(1,1);  % changed this 2/23/19 from length(interventionSeq)-1        
                singleIntv = cell(1, bnet.nNodes);
                singleIntv{interveneNode} = intvValue;
                interventions{1} = singleIntv;
            end
            if (sampleSachs == 1 || sampleSachs == 2)
                interventions = interventionSeq;
            end
        end % end creating obs or intv data for this experiment 
    end
    
    [bnet, entropy, dg, hamming, seqData, seqClamped, remData, remClamped] = runExperiment2(bnet, nObservationCases, interventions, nInterventionCases, seqData, seqClamped, remData, remClamped );   
   
 
    HExp(:, experiment) =  entropy.H;
    if(experiment < stg.maxExp) % need this since on last experiment bnet.eligibleNodes might be emptyset and cause error (esp.when sampleSachsData = 1)
        maxH(experiment) =  max(entropy.H(bnet.eligibleNodes)); 
    end
    postE(experiment) = entropy.postEntropy; 
    tpr(experiment) = dg.tpr;  
    fpr(experiment) = dg.fpr;
    tnr(experiment) = dg.tnr;
    fnr(experiment) = dg.fnr;
    hammingMean(experiment) = hamming.mean;
    hammingVar(experiment) = hamming.var;    
    % MAKE SURE THAT THE INTERVENED ON NODE IS REMOVED FROM SET OF ELIGIBLE
    % NODES
    if (experiment > 1 && sampleSachs == 0)
        bnet.eligibleNodes(bnet.eligibleNodes == interveneNode) = [];
    end
    experiment = experiment + 1;
end % end while loop for running experiments
csvwrite(sprintf('%s/interventionSeq.txt',simPath),interventionSeq);
csvwrite(sprintf('%s/HExp.csv',simPath),HExp);

%% make the interventionSeq vector the same length as stg.maxExps 
if length(interventionSeq) < stg.maxExp
    interventionSeq(end+1:stg.maxExp) = NaN;
    tpr(end+1:stg.maxExp) = NaN;
    fpr(end+1:stg.maxExp) = NaN;
    tnr(end+1:stg.maxExp) = NaN;
    fnr(end+1:stg.maxExp) = NaN;
    postE(end+1:stg.maxExp) = NaN;
    maxH(end+1:stg.maxExp) = NaN;
end    
diagnostics.tpr = tpr;
diagnostics.fpr = fpr;
diagnostics.tnr = tnr;
diagnostics.fnr = fnr;

end % end function 
