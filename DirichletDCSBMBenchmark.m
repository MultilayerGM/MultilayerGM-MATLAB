function [A,S]=DirichletDCSBMBenchmark(nodes,layers,types,dependencyMatrix,varargin)
% Generate multi-aspect multilayer networks with a DCSBM planted partition model
%
% Input:
%
%   nodes: Number of physical nodes
%
%   layers: Number of layers for each aspect
%
%   types: 'char' vector specifying the update type for each aspect, 'o' for
%       ordered and 'r' for random.
%
%   dependencyMatrix: matrix of copying probabilities
%
% Output:
%
%   A: cell array of adjacency matrices for each layer
%
%   S: planted multilayer partition
%
% Options:
%
%   UpdateSteps: [100] number of Gibbs updates before returning partition
%
%   theta: [1] concentration parameter for Dirichlet null distribution
%
%   communities: [10] number of communities
%
%   q: [1] probability for a community to be active in a given layer
%
%   exponent: [-3] powerlaw exponent for expected degree distribution
%
%   kmin: [3] minimum expected degree
%
%   kmax: [50] maximum expected degree
%
%   mu: [0.1] fraction of random edges
%
%   maxreject: [100] maximum number of rejections before bailing out and
%       issuing a warning (the resulting network has less than the desired
%       number of edges)
%
% see also: PartitionGenerator, DirichletNullDistribution,
% DCSBMNetworkGenerator, MultiplexMatrix, TemporalMatrix,
% BlockMultiplexMatrix
%
% Version: 1.0.1
% Date: Tue  4 Jul 2017 16:38:06 BST
% Author: Lucas Jeub
% Email: ljeub@iu.edu
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


% set defaults for all options
parsePartitionOptions=inputParser();
parsePartitionOptions.KeepUnmatched=true;
addParameter(parsePartitionOptions,'UpdateSteps',[]);
addParameter(parsePartitionOptions,'InitialPartition',[])
parse(parsePartitionOptions,varargin{:});
PartitionOptions=NonDefaultOptions(parsePartitionOptions);

parseNullOptions=inputParser();
parseNullOptions.KeepUnmatched=true;
addParameter(parseNullOptions,'theta',[]);
addParameter(parseNullOptions,'communities',[]);
addParameter(parseNullOptions,'q',[]);
parse(parseNullOptions,parsePartitionOptions.Unmatched);
NullOptions=NonDefaultOptions(parseNullOptions);

parseNetworkOptions=inputParser();
addParameter(parseNetworkOptions,'exponent',[]);
addParameter(parseNetworkOptions,'kmin',[]);
addParameter(parseNetworkOptions,'kmax',[]);
addParameter(parseNetworkOptions,'mu',[]);
addParameter(parseNetworkOptions,'maxreject',[]);
parse(parseNetworkOptions,parseNullOptions.Unmatched);
NetworkOptions=NonDefaultOptions(parseNetworkOptions);

% generate partition using Dirichlet null distribution
S=PartitionGenerator(nodes,layers,types,dependencyMatrix,...
    DirichletNullDistribution(layers,NullOptions),PartitionOptions);

% generate intralayer edges using DCSBM benchmark model
A=DCSBMNetworkGenerator(S,NetworkOptions);

end

function options=NonDefaultOptions(input)
    fields=fieldnames(input.Results);
    options=struct();
    for i=1:numel(fields)
        if ~ismember(fields{i},input.UsingDefaults)
            options.(fields{i})=input.Results.(fields{i});
        end
    end
end
