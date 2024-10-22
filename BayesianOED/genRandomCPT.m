function [bnet] = genRandomCPT(bnet)
% genRandomCPT: generates a conditional probability table (CPT) given
%               a causal network structure (bnet) and a desired 
%               degree of randomness.

% from the Bayes Net Toolbox help documentation: 
%   To control the degree of randomness (entropy), 
%   you can sample each row of the CPT from a Dirichlet(p,p,...) distribution. 
%   If p << 1, this encourages "deterministic" CPTs (one entry near 1, the rest near 0).
%   If p = 1, each entry is drawn from U[0,1]. 
%   If p >> 1, the entries will all be near 1/k, 
%   where k is the arity of this node, i.e., each row will be nearly uniform.
    p = 0.8;
    for i=1:bnet.nNodes
        k = bnet.ns(i);
        ps = parents(bnet.dag, i);
        psz = prod(bnet.ns(ps));
        CPT = dirichlet_sample(p*ones(1,k), psz); % this is in the Foreign directory (looks like it was previously called sample_dirichlet)
        %bnet.CPD{i} = tabular_CPD(bnet, i, 'CPT', CPT);      
        bnet.CPD{i} = tabular_CPD(bnet, i, 'CPT', CPT, 'dirichlet_type', 'BDeu' );
    end
% uncomment below if you want to inspect the CPTs for a given node
%     for i=1:bnet.nNodes
%         disp(i)
%         values = struct2cell(bnet.CPD{i}); values{1}
%     end
end

