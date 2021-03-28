function [ interveneNode] = nextIntervention(H, eligibleNodes)

%UNTITLED5 Summary of this function goes here
%   The additional argument "eligibleNodes" is a vector that specifies
%   which of the nodes in the network you can actually perturb (or we have
%   interventional data to sample from).

global chooseIntervention;
global experiment;
global randNodeSeq;
global sim;
if chooseIntervention == 0
    interveneNode = 0;
    fprintf('Sampling observational data \n')
elseif chooseIntervention == 1
    interventionSeq = [0 2 7 4 8 9]; %[0 8 9 2 4 7]; %[0 8 9 2 3 7 5 11 10]; % Bninfo selected batch %[0 8 9 4 2 7];
    interveneNode = interventionSeq(experiment);
    fprintf('\n Fixed sequence intervention on node %d \n', interveneNode)    
elseif chooseIntervention == 2
    % get random node from array initialized at start of sim
    interveneNode = randNodeSeq(sim, experiment-1); 
    %interveneNode = randi([1 nNodes],1)
    fprintf('\n Randomly chosen intervention on node %d \n', interveneNode)
elseif chooseIntervention == 11 %Bninfo sequence
    % get random node from array initialized at start of sim
    interveneNode = randNodeSeq(sim, experiment-1); 
    %interveneNode = randi([1 nNodes],1)
    fprintf('\n Bninfo chosen intervention on node %d \n', interveneNode)
else % for any method that involves an entropy measure for each node
    [~, idx] = max(H(eligibleNodes));
    interveneNode = eligibleNodes(idx);
    fprintf('\n Adaptively chosen intervention on node %d, using method %d. \n', interveneNode, chooseIntervention)
end



