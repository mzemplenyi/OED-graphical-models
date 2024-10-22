function [ interveneNode] = nextIntervention(H, eligibleNodes, method, sim)

%UNTITLED5 Summary of this function goes here
%   The additional argument "eligibleNodes" is a vector that specifies
%   which of the nodes in the network you can actually perturb (or we have
%   interventional data to sample from).

global experiment;
global randNodeSeq;
if isequal(method,'Observational')
    interveneNode = 0;
    fprintf('Sampling observational data \n')
elseif isequal(method,'Fixed')
    interventionSeq = [0 2 7 4 8 9]; 
    interveneNode = interventionSeq(experiment);
    fprintf('\n Fixed sequence intervention on node %d \n', interveneNode)    
elseif isequal(method,'Random')
    % get random node from array initialized at start of sim
    interveneNode = randNodeSeq(sim, experiment-1); 
    fprintf('\n Randomly chosen intervention on node %d \n', interveneNode)
elseif isequal(method,'bninfo') %Bninfo sequence
    % get node to intervene on using array initialized at start of sim
    interveneNode = randNodeSeq(sim, experiment-1); 
    fprintf('\n NOTE: For the "bninfo" method, specify the intervention sequence using the randNodeSeq array \n')
    fprintf('\n Bninfo chosen intervention on node %d \n', interveneNode)
else % for any method that involves an entropy measure for each node
    [~, idx] = max(H(eligibleNodes));
    interveneNode = eligibleNodes(idx);
    fprintf('\n Adaptively chosen intervention on node %d, using method: %s. \n', interveneNode, method)
end



