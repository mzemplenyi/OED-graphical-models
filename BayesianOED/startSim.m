function [interventionSeq, diagnostics, postE, maxH, hammingMean, hammingVar] = startSim(stg, sim )

simPath = sprintf('%sSim%d',stg.resPath, sim);
mkdir(simPath)
bnet = stg.bnet;


global experiment; experiment = 1;
global nObsCases; global nInitialObs;
global nIntvCases;
global sampleSachs;


if sampleSachs == 2  % sample from simulated Sachs data with interventions on all 11 nodes
    %cd('C:/Users/Michele/Documents/Jeff Miller/BDAGL/BDAGLMZ/Sim Data/EightTree_200/');
    cd('C:/Users/Michele/Documents/Jeff Miller/BDAGL/BDAGLMZ/Sim Data/EightLine/');
    %cd('C:/Users/Michele/Documents/Jeff Miller/BDAGL/BDAGLMZ/Asia Sim Data/300obs300_yes/');
    %cd('C:/Users/Michele/Documents/Jeff Miller/BDAGL/BDAGLMZ/Sachs Sim Data/1800obs600/');
    datafile = sprintf('data_%d.csv', sim);
    data = csvread(datafile);
   %datafile       = dir(sprintf('C:/Users/Michele/Documents/Jeff Miller/BDAGL/BDAGLMZ/Sachs Sim Data/500obs100/data_%d.csv',sim));
    %data = csvread(datafile.name);
    %data = csvread('C:/Users/Michele/Documents/Jeff Miller/BDAGL/BDAGLMZ/Sachs Sim Data/sachsData_noheader20190812.csv');
    data = data';
    %clampedfile       = dir(sprintf('C:/Users/Michele/Documents/Jeff Miller/BDAGL/BDAGLMZ/Sachs Sim Data/500obs100/clamped_%d.csv',sim));
    clampedfile = sprintf('clamped_%d.csv', sim);
    clamped = csvread(clampedfile);
    %clamped = csvread('C:/Users/Michele/Documents/Jeff Miller/BDAGL/BDAGLMZ/Sachs Sim Data/sachsClamped_noheader20190812.csv');
    clamped = clamped';
    remData = data;
    remClamped = clamped;
end   

stopThreshold = -0.1;
uncertainty = 1;
%fixedIntvValue = 1; % 0 = adaptive intvValue, 1 = fixed intvValue
intvValue = 1; % 1 is low-level, 2 is high-level



% initialize matrices that will store the sequential data and clamped data
seqData = NaN(bnet.nNodes, 0);
seqClamped = NaN(bnet.nNodes, 0);

%% 
cd 'C:/Users/Michele/Documents/Jeff Miller/BDAGL/BDAGLMZ/';
interventionSeq = NaN(1,stg.maxExp);interventionSeq(1) = 0; % vector to store which experiments are run, first experiment is just observational 
tpr = NaN(1,stg.maxExp);
fpr = NaN(1,stg.maxExp);
tnr = NaN(1,stg.maxExp);
fnr = NaN(1,stg.maxExp);
postE = NaN(1,stg.maxExp);
%skelE = NaN(1,stg.maxExp);
maxH = NaN(1,stg.maxExp);
HExp = NaN(bnet.nNodes, stg.maxExp); 
hammingMean = NaN(1,stg.maxExp);
hammingVar = NaN(1,stg.maxExp);

while(experiment <= stg.maxExp)
    if experiment == 1
    %% Generate only observational data
        nObservationCases = nInitialObs; % # observational data cases
        nInterventionCases = 0; % no interventions
        if sampleSachs == 0
            interventions = {};
        end
        if (sampleSachs == 1 || sampleSachs == 2)
            interventions = interventionSeq;
        end
    else %% make interventions cell array
        uncertainty = entropy.postEntropy; % posterior entropy on graphs       
        if uncertainty < stopThreshold
            sprintf('\n Uncertainty of %d dropped below stop threshold. Not running experiment %d \n', uncertainty, experiment)
            break;
        end    
        fprintf('\n Beginning Experiment %d \n', experiment)
        interveneNode = nextIntervention(HExp(:,experiment-1), bnet.eligibleNodes); 
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
end
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
