function [S,SIntermediate]=OrderedRandomPartitionGenerator(nodes,layers,types,dependencyMatrix,nullDistribution,varargin)
% Generate partition for a multilayer network
% with specified interlayer dependencies.
%
% Input:
%
%   nodes: Either a scalar specifying the number of physical nodes or
%       a |V_M|x(1+d) matrix specifying each state node in the format
%       [node, aspect_1,...,aspect_d] (the latter format allows for missing
%       nodes in some layers)
%
%   layers: Vector giving the number of elements for each aspect in the
%       format [l_1,...,l_d]. Note that for a multilayer network with a
%       single aspect this is just a scalar giving the number of layers.
%
%   types: 'char' vector specifying the update type for each aspect, 'o' for
%       ordered and 'r' for random.
%
%   dependencyMatrix: A matrix of copying probabilities, corresponding to
%       the flattened interlayer dependency tensor. This matrix is either
%       of size lxl for the layer-coupled case or of size (n*l)x(n*l) in
%       the general case. If state nodes are given explicitly in 'nodes',
%       'dependencyMatrix' should be of size |V_M|x|V_M|, giving only the
%       probabilities for state nodes that are actually present in the
%       network.  The ways the matrix encodes the tensor are as follows:
%
%       layer-coupled case:
%           Mapping for flattening the indices:
%
%           a=aspect_1+l_1*(aspect_2-1)+...+l_1*l_2*...*l_(d-1)*(aspect_d-1)
%
%           where 'dependencyMatrix(a,b)' is the probability that a node in
%           layer b copies its community assignment from the same node in
%           layer a. The rows of 'dependencyMatrix' should sum to a value<1,
%           where 1-sum(transtionMatrix(:,b)) is the probability of
%           resampling the community assignment from the specified null
%           distribution for a node in layer b.
%
%       general case:
%           Mapping for flattening the indeces when state nodes are not
%           explicitly specified:
%
%           i=node+n*(aspect_1-1)+n*l_1*(aspect_2-1)+..._n*l_1*...*l_(d-1)*(aspect_d-1)
%
%           When state nodes are given explicitly, the transition matrix
%           should have the same number of rows as the matrix of state
%           nodes, where the ith state node corresponds to the ith row of
%           the transition matrix.
%
%           'dependencyMatrix(i,j)' is the probability that
%           state node j copies its community assignment from state node i.
%           The rows of 'dependencyMatrix' should sum to a value < 1,
%           where 1-sum(transtionMatrix(:,j)) is the probability of
%           resampling the community assignment from the specified null
%           distribution for node j.
%
%   nullDistribution: A function that takes state nodes (i.e., row vectors
%       of the form [node, aspect_1,...,aspect_d] ) as input and returns a
%       random community assignment.
%
%   initialPartition: Optionally specify starting partition (array of
%       dimension nxl_1x...xl_d).
%
% Output:
%
%   S: Multilayer partition after the last update step (returned as an
%       array of dimension nxl_1x...xl_d. In the case where state nodes are
%       specified directly, missing state nodes have community label 0.
%
%   SIntermediate: Cell array of multilayer partitions at intermediate
%       steps of the Markov chain (see 'IntermediateSteps'). First entry is
%       always the initial partition and last entry is the final partition
%       (the same as 'S').
%
% Options:
%
%   UpdateSteps: [default: 100] number of Gibbs updates before returning
%       partition (i.e., the number of times each state node's community
%       assignment is updated). This defaults to a single update if the
%       provided transition matrix is fully ordered (i.e. has zero
%       lower-triangular part.
%
%   IntermediateSteps: [default: []] number of Gibbs updates
%       between sampled partitions from the same chain returned in
%       'SIntermediate'. The initial partition and final partition are
%       always included in 'SIntermediate'. Useful for looking at
%       autocorrelation features of the Markov chain to verify it is
%       consistent with having converged.
%
%
%
% Version: 1.0.1
% Date: Tue  4 Jul 2017 16:38:06 BST
% Author: Lucas Jeub
% Email: ljeub@iu.edu
%
%
% References:
%
%       [1] Generative benchmark models for mesoscale structure in multilayer
%       networks, M. Bazzi, L. G. S. Jeub, A. Arenas, S. D. Howison, M. A.
%       Porter. arXiv1:608.06196.
%
% Citation:
%
%       If you use this code, please cite as
%       Lucas G. S. Jeub and Marya Bazzi
%       "A generative model for mesoscale structure in multilayer networks
%       implemented in MATLAB," https://github.com/MultilayerBenchmark/MultilayerBenchmark (2016).


% Parse input
parseArgs=inputParser();
addParameter(parseArgs,'IntermediateSteps',[]);
addParameter(parseArgs,'InitialPartition',[]);
addParameter(parseArgs,'UpdateSteps',100);
parse(parseArgs,varargin{:});
options=parseArgs.Results;
options.isset=@(s) ~isempty(options.(s));

% determine network shape and set up state nodes if not given explicitly
l=prod(layers);
if numel(nodes)==1
    n=nodes;
    nodes=ind2subarray([n,layers],1:n*l);
else
    n=max(nodes,[],1);
end
nvm=size(nodes,1);
oind=find(types=='o');
ordered=layers(oind);
ol=prod(ordered);
nodeUpdateBlocks=cell(ol,1);
for i=1:nvm
    oi=subarray2ind(ordered,nodes(i,oind+1));
    nodeUpdateBlocks{oi}=[nodeUpdateBlocks{oi};nodes(i,:)];
end

% setup map from state nodes to rows of transition matrix
transitionMap=zeros(n*prod(layers));
if isequal(size(dependencyMatrix),[l,l])
        transitionMap(subarray2ind([n,layers],nodes))=subarray2ind(layers,nodes(:,2:end));
elseif isequal(size(dependencyMatrix),[nvm,nvm])
    transitionMap(subarray2ind([n,layers],nodes))=(1:nvm)';
else
    error('CommunityStructureGenerator:dependencyMatrix:size','Specified transition matrix is of inconsistent size');
end


% Sample partitions
S=zeros([n,layers]);
if options.isset('IntermediateSteps')
    SIntermediate=cell(ceil(options.UpdateSteps/options.IntermediateSteps)+1,1);
else
    SIntermediate=cell(2,1);
end

if options.isset('InitialPartition')
    S=options.InitialPartition;
else
    for i=1:size(nodes,1)
        S(subarray2ind([n,layers],nodes(i,:)))=nullDistribution(nodes(i,:));
    end
end
SIntermediate{1}=S;
usteps=options.UpdateSteps;

if options.isset('IntermediateSteps')
    totalsteps=0;
    isteps=options.IntermediateSteps;
    it=1;
    while totalsteps+isteps<=usteps
        S=GibbsPartitionSampler(S,nodeUpdateBlocks,transitionMap,dependencyMatrix,nullDistribution,isteps);
        totalsteps=totalsteps+isteps;
        it=it+1;
        SIntermediate{it}=S;
    end

    if totalsteps<usteps
        S=GibbsPartitionSampler(S,nodeUpdateBlocks,transitionMap,dependencyMatrix,nullDistribution,usteps-totalsteps);
        SIntermediate{it+1}=S;
    end
else
    S=GibbsPartitionSampler(S,nodeUpdateBlocks,transitionMap,dependencyMatrix,nullDistribution,usteps);
    SIntermediate{2}=S;
end

end

% Gibbs sampling
function S=GibbsPartitionSampler(S,nodeUpdateBlocks,transitionMap,dependencyMatrix,nullDistribution,steps)
% Run Gibbs sampling
%
% Inputs:
%
%   S: current multilayer partition
%
%   nodes: matrix of statenodes (each state node is a row) with format
%       [node,aspect_1,...,aspect_d]
%
%   transitionMap: vector of indeces mapping statenodes to the
%       corresponding row/column of the dependencyMatrix
%
%   dependencyMatrix: coupling edges (sum(dependencyMatrix,1)<=1)
%
%   nullDistribution: function (statenode->random community assignment)
%
%   steps: number of update steps
%
% Output:
%
%   S: updated partiton after 'steps' passes over all state nodes

not_resample=full(sum(dependencyMatrix,1));
LW=cumsum(dependencyMatrix,1);

size_spec=size(S);
if length(dependencyMatrix)==prod(size_spec(2:end))
    % layer-coupled
    for lind=1:length(nodeUpdateBlocks)
        nodes=nodeUpdateBlocks{lind};
        nodeindmap=subarray2ind(size_spec,nodes);
        for i=1:steps
            for j=randperm(size(nodes,1))
                % nodeind=subarray2ind(size_spec,nodes(j,:));
                nodeind=nodeindmap(j);
                decide=rand();
                if decide<not_resample(transitionMap(nodeind))
                    S(nodeind)=S(nodes(j,1),ind2sub(size_spec(2:end),find(LW(:,transitionMap(nodeind))>decide,1)));
                else
                    S(nodeind)=nullDistribution(nodes(j,:));
                end
            end
        end
    end
else
    % general case
    for lind=1:length(nodeUpdateBlocks)
        nodes=nodeUpdateBlocks{lind};
        nodeindmap=subarray2ind(size_spec,nodes);
        for i=1:steps
            for j=randperm(size(nodes,1))
                % nodeind=subarray2ind(size_spec,nodes(j,:));
                nodeind=nodeindmap(j);
                decide=rand();
                if decide<not_resample(transitionMap(nodeind))
                    S(nodeind)=S(subarray2ind(size_spec,nodes(find(LW(:,transitionMap(nodeind))>decide,1),:)));
                else
                    S(nodeind)=nullDistribution(nodes(j,:));
                end
            end
        end
    end
end
end
