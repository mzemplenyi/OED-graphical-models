 function [ H] = entropyMEC(bnet, samples )
% Implments the "MEC" method described by Jeff Miller
% For each sampled DAG from the posterior and each candidate intervention
% modify the sampled DAG by removing all edges from parent nodes to the 
% candidate intervention node, e, and then find the Markov equivalence
% class (MEC) that corresponds to that modified DAG. This MEC is 
% represented by a CPDAG where:
% 	If the edge is compelled then 1 on the edge.
%	If the edge is reversible then 1 on the edge and 1 in the reverse edge.

nNodes = bnet.nNodes;
N = samples.nSamples;
mecEntropy = NaN(1,nNodes); % initialize vector to store entropy over the partition of MECs relative to each node
for node = 1:nNodes
    % initialize hashtable to store possible cpdags
    cpdags = [];
    cpdags.HT = java.util.Hashtable(2^(nNodes-1)); 
    % get keys from the DAGs sampled in the MCMC
    visited = 0;
    keys = samples.HT.keys; 
    while keys.hasMoreElements()
        sk = keys.nextElement();
        dagSampled = char2dag(sk, bnet.nNodes); % adj matrix
        dagSampled(:,node) = 0; % set candidate intervention node to have no parents (removes incoming edges)
               
        cpdagSampled = dag_to_cpdag(dagSampled); % dagSampled has been modified to remove edges from pa(e) to e
        sVal = samples.HT.get(sk);
       % if child == 2; figure; myDrawGraph(dagSampled,'labels', bnet.names); title(sprintf('No. times visited: %d', dagVal(1))); end;
       %psMat = zeros(nNodes); % need to have a matrix in order to use dag2char
       % ps = dagSampled(:,child); % get parent set for the child of interest   
      %  psMat(:,child) = ps; % arbitrarily putting parent vector in "child" column of matrix -- note don't always put it in the same column (I ran into issues  getting distinct hash keys when I did that)
        cpdagKey = dag2char(cpdagSampled); % psKey is the representation of the psMat adj matrix as a character
		cpdagValue = cpdags.HT.get(cpdagKey);
        
        if isempty(cpdagValue) % this cpdag has never been visited before 
			count = sVal(1); % number of times the corresponding dag was visited  
		else
			count = cpdagValue(1) + sVal(1); % update the count of times this cpdag has been visited
        end
        visited = visited + sVal(1); % by end of while loop this should equal the number of requested MCMC samples
		cpdags.HT.put( cpdagKey, count ); % store parent vector and times visited in the MEC hashtable	
        
    end % end of looping through hashtable of all sampled graphs
%     allKeys = arrayfun(@char, parentSets.HT.keySet.toArray, 'UniformOutput', false); % print keys
%     cellfun(@(x) parentSets.HT.get(x), allKeys, 'UniformOutput', false) % print key values
%    parentSets.HT
    
    % calculate entropy TS equivalence classes (the cpdags hashtable)

    timesVisited = [0]; % vector that will store total number of times each unique cpdag was visited 
    nPart = 0; % counter as we move through keys of the cpdag HT
    partKeys = cpdags.HT.keys; % all keys from hashtable with the parts of the partition
    while partKeys.hasMoreElements()
        nPart = nPart + 1;
        pKey = partKeys.nextElement();
        %dagSampled = char2dag(dagKey, samples.nNodes);
        partVisited = cpdags.HT.get(pKey);
        timesVisited(nPart) =  partVisited(1);
       % nSubsets = nSubsets + dagVal(1);
    end
      
    % empirical probability for each visited parent set out of all graphs
    % sampled
    p = timesVisited / N;
   % csvwrite(sprintf('%s/Sim%d/cpdagProb_Exp%d.csv',scenPath, sim, experiment),p.');
    

    logp = log(p);
    mecEntropy(node)= -1*dot(p,logp);
end % end for loop over nodes
H = mecEntropy;
end % end function

