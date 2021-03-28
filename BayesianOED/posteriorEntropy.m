function [postEntropy] = posteriorEntropy(samples)
% posteriorEntropy: calculates the entropy over the posterior distribution on graphs


count = [0];
nSamples = 0;
uniqueGC = 0;
keys = samples.HT.keys;
while keys.hasMoreElements()
    uniqueGC = uniqueGC + 1;
    dagKey = keys.nextElement();
    %dagSampled = char2dag(dagKey, samples.nNodes);
    dagVal = samples.HT.get(dagKey);
    count(uniqueGC) =  dagVal(1);
    nSamples = nSamples + dagVal(1);
end

% empirical probability for each visited graph 
postG = count/nSamples;

logPostG = log(postG);
postEntropy = -1*dot(postG,logPostG);
% postEntropy = -1*dot(postG,logPostG) / log(uniqueGC);
end