function [ H] = entropyDP(nNodes, ep)
% Michele Zemplenyi
% entropyDP calculate "non-symmetrical entropy" as described in Li and Leong
% using edge probabilities calculated by DP instead of from MCMC samples


%Hvec = NaN(1,nNodes*nNodes);
edgeProb = reshape(ep, 1, nNodes*nNodes);
calcEntropy = @(x) -x*log2(x) - (1-x)*log2(1-x); % matlab function handle: https://www.mathworks.com/help/matlab/matlab_prog/creating-a-function-handle.html
Hvec = arrayfun(calcEntropy, edgeProb);
Hmat = reshape(Hvec, nNodes, nNodes);
H = sum(Hmat, 2, 'omitnan'); % sum of rows

% %% test
% nNodes = 2;
% ep = [ 0 0.25; 0.75 0];
% Hvec = NaN(1,nNodes*nNodes);
% edgeProb = reshape(ep, 1, nNodes*nNodes);
% Hvec = (-edgeProb .* log(edgeProb)
% Hvec = -(edgeProb .* log2(edgeProb)) + (1-edgeProb).* log2(1-edgeProb);
% -pChild*log2(pChild) - (1-pChild)*log2(1-pChild);
% Hmat = reshape(Hvec, nNodes, nNodes);
% -0.4 * log2(0.4) - 0.6 * log2(0.6)
% -0.6 * log2(0.6) - 0.4 * log2(0.4)

end