
addpath('Sim Scenarios');
addpath('Shared');
addpath('StructureMcmc');
addpath('DP');
addpath('OrderMcmc');
addpath('ExactEnumeration');
addpath('DagHashTable');
addpath('ADTree');
addpath('DBN');
addpath('demos2');
addpath('OptimalMAP');
addpath('Foreign/causalExplorer/');
addpath('Foreign/causalExplorer/PCodes/');
addpath('Foreign');
addpath('Mex');

mkMex('Mex');
mkMex('DagHashTable');
mkMex('DP');
mkMex('OrderMcmc');

if exist('gammaln')~=3
	cd('Foreign');
	mex -compatibleArrayDims gammaln.c minka_mexutil.c minka_util.c
	cd('..');
end 

cd('ADTree');
if  exist('mkADTree')~=3
	mex -O mkADTree.c util.c
end
if  exist('mkContab')~=3
	mex -compatibleArrayDims mkContab.c
end
if  exist('testADTree')~=3
	mex -O testADTree.c
end
if  exist('maxNumADTrees')~=3
	mex -O maxNumADTrees.c
end
cd('..');

