function [diagnostics] = checkDiagnostics(PCmat, trueDag, simPath, expNum)
%global trueDag
% checkTPR: calculates the true positive edge rate for the estimated average DAG
% adjacency matrix matches the true DAG
%   input: PCmat is a square nNodes x nNodes matrix with entries equal to the 
%           empirical proportion of sampled graphs in which the column node
%           was a child of the row node
%   output: true positive rate, the proportion of correctly detected edges
%           among the edges in the ground truth network

     % print non-thresholded adj mat for run with just observational data
%     if experiment == 1
%         csvwrite(sprintf('%s/Sim%d/PCmat_Exp%d.csv',scenPath, sim, experiment),PCmat);
%     end    
    % create thresholded matrix for the PCmat
    dim = size(PCmat, 1);    
    thr = 0.5; % threshold for edge probability to become a 1 or 0
    PCvec = reshape(PCmat, 1, []);
    estimate = NaN(dim*dim,1);
    estimate(PCvec >= thr) = 1;
    estimate(PCvec < thr) = 0; 
    
    % 
    gt = reshape(trueDag, 1, []); % ground truth network as vector
    % True Positive Rate of Edge Detection
    pos = find(gt == 1);    
    tp = (sum(estimate(pos) == 1));   
    tpr = tp / length(pos);
    
    % False Positive Rate of Edge Detection
    neg = find(gt == 0);
    fp = (sum(estimate(neg) == 1));
    fpr = fp / length(neg);
    
    % True Negative Rate of Edge Detection
    tn = (sum(estimate(neg) == 0));
    tnr = tn / length(neg);
    
    % False Negative Rate of Edge Detection
    fn = (sum(estimate(pos) == 0));
    fnr = fn / length(pos);
    
    % calculate hamming distance 
    estimate= reshape(estimate, dim, dim);
    hammingDist = nnz(xor(estimate,trueDag)); 
    
    %% store all diagnostics
    diagnostics.tpr = tpr;
    diagnostics.fpr = fpr;
    diagnostics.tnr = tnr;
    diagnostics.fnr = fnr;
    diagnostics.hammingDist = hammingDist;
    
    
%     % True Negative Rate of Edge Detection
%     tnEdges = find(gt == 0);    
%     tnSum = (sum(estimate(tnEdges) == 0));   
%     tnr = tnSum / length(tnEdges);
%     
%     % False Negative Rate of Edge Detection       
%     fnSum = (sum(estimate(tpEdges) == 0));   
%     fnr = fnSum / length(tnEdges);
    
     
    % save a copy of the thresholded PCmat 
    csvwrite(sprintf('%s/PCmat_Exp%d.csv',simPath, expNum),PCmat);
%     if match == 0 % store the PCmat if algorithm arrives at wrong DAG
%         csvwrite(sprintf('%s/Sim%d/PCmat_Exp%d.csv',scenPath, sim, experiment),PCmat);
%     end
end % end function

