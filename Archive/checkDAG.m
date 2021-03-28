function [ match] = checkDAG(PCmat, trueDag)
%global trueDag
global scenPath
global sim
global experiment
% checkDAG: checks whether a thresholded version of the posterior sampled
% adjacency matrix matches the true DAG
%   input: PCmat is a square nNodes x nNodes matrix with entries equal to the 
%           empirical proportion of sampled graphs in which the column node
%           was a child of the row node
%   output: binary indicator of whether the thresholded adj matrix
%           matches the true data generating DAG

     % print non-thresholded adj mat for run with just observational data
    if experiment == 1
        csvwrite(sprintf('%s/Sim%d/PCmat_Exp%d.csv',scenPath, sim, experiment),PCmat);
    end 
    dim = size(PCmat, 1);    
    thr = 0.5; % threshold for edge probability to become a 1 or 0
    PCvec = reshape(PCmat, 1, []);
    estimate = NaN(dim*dim,1);
    estimate(PCvec >= thr) = 1;
    estimate(PCvec < thr) = 0; 
    PCmat = reshape(estimate, dim, dim);  % put back in matrix form for csv
    
    %%% old code     
%     dim = size(PCmat, 1);
%     thr = 0.5; % threshold
%     for r = 1:dim
%         for c = 1:dim
%             if PCmat(r,c) >= thr
%                 PCmat(r,c) = 1;
%             else
%                 PCmat(r,c) = 0;
%             end 
%         end
%     end 
    match = isequal(trueDag, PCmat);
    % save a copy of the thresholded PCmat 
    csvwrite(sprintf('%s/Sim%d/PCmat_Exp%d.csv',scenPath, sim, experiment),PCmat);
%     if match == 0 % store the PCmat if algorithm arrives at wrong DAG
%         csvwrite(sprintf('%s/Sim%d/PCmat_Exp%d.csv',scenPath, sim, experiment),PCmat);
%     end
end % end function

