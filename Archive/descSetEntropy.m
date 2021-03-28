function [ H] = descSetEntropy(bnet, samples )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global scenPath; global sim; global experiment;
nNodes = bnet.nNodes;
% global psCutoff;
% global psMAP;
N = samples.nSamples;
dsEntropy = NaN(1,nNodes); % initialize vector to store desc set entropy for each node
for node = 1:nNodes
    % initialize hashtable to store possible desc vectors
    descSets = [];
    descSets.HT = java.util.Hashtable(2^(nNodes-1)); 
    % get keys from the DAGs sampled in the MCMC
    visited = 0;
    keys = samples.HT.keys; 
    while keys.hasMoreElements()
        sk = keys.nextElement();
        dagSampled = char2dag(sk, bnet.nNodes); % adj matrix        
        ssVal = samples.HT.get(sk);
       % if child == 2; figure; myDrawGraph(dagSampled,'labels', bnet.names); title(sprintf('No. times visited: %d', dagVal(1))); end;
        descMat = zeros(nNodes); % need to have a matrix in order to use dag2char
        % find descendants of node (taken from entropy.m)
        paths = expm(dagSampled); 
        ds = paths(node,:)>0; % store desc set for the node of interest (note this set includes the node itself)   
        descMat(:,node) = ds; % arbitrarily putting desc vector in "node" column of matrix -- note don't always put it in the same column (I ran into issues  getting distinct hash keys when I did that)
        dsKey = dag2char(descMat); % psKey is the representation of the psMat adj matrix as a character
		dsValue = descSets.HT.get(dsKey);
        
        if isempty(dsValue)
			count = ssVal(1); % number of times that dag was visited  
		else
			count = dsValue(1) + ssVal(1); %if desc vector already sampled previously, update the count of times it has been visited
        end
        visited = visited + ssVal(1); % by end of while loop this should equal the number of requested MCMC samples
		descSets.HT.put( dsKey, count ); % store parent vector and times visited in the parent vec hashtable	
        
    end % end of looping through hashtable of all sampled graphs
%     allKeys = arrayfun(@char, parentSets.HT.keySet.toArray, 'UniformOutput', false); % print keys
%     cellfun(@(x) parentSets.HT.get(x), allKeys, 'UniformOutput', false) % print key values
%    parentSets.HT
    
    % calculate entropy over desc sets, pass in hashtable for descSets
    % dsEntropy(child) = posteriorEntropy(parentSets); 
    timesVisited = [0];
    nUnique = 0;
    dsKeys = descSets.HT.keys; % keys from hashtable with all parent subsets
    while dsKeys.hasMoreElements()
        nUnique = nUnique + 1;
        dsk = dsKeys.nextElement();
        %dagSampled = char2dag(dagKey, samples.nNodes);
        ssVal = descSets.HT.get(dsk);
        timesVisited(nUnique) =  ssVal(1);
       % nSubsets = nSubsets + dagVal(1);
    end
      
    % empirical probability for each visited desc set out of all graphs
    % sampled
    p = timesVisited / N;
    %csvwrite(sprintf('%s/Sim%d/dsProb_Exp%d.csv',scenPath, sim, experiment),p.');
    
%     if psMAP == 1
%         a = 526; % parameter we set, alpha = a / K, a = 0 reduces to empirical probabilities
%         K = 2^(nNodes-1);
%         % calculate phat_s for all visited parent sets
%         phat_v = (timesVisited + (a/K)) / (N + a); % (this is a vector)
%         % phat_s for unvisited parent sets
%         phat_uv = (a/K) / (N + a); % (this is a scalar)
%         % number of unvisited parent sets ( |S^c| )
%         nUV = K - nUnique;
%         phat_uv_all = repelem(phat_uv, nUV); % (now a vector)
%         % build vector with both phat_v and phat_uv
%         p = [phat_v phat_uv_all];
%         csvwrite(sprintf('%s/Sim%d/psMAP_Exp%d.csv',scenPath, sim, experiment),p.');
%     end % end parentSetMAP
%     % implement cutoff
%     if psCutoff ~= 0
%         p = p(p >= psCutoff); % get just the probabilities that exceed cutoff
%         p = p / sum(p); % renormalize to account for the prob mass that we removed
%     end
    logp = log(p);
    dsEntropy(node)= -1*dot(p,logp);
end % end for loop over nodes
H = dsEntropy;
end % end function

