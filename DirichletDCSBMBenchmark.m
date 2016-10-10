function [A,S]=DirichletDCSBMBenchmark(nodes,layers,transitionMatrix,varargin)
% Generate uniform multiplex benchmark networks
%
% Input:
% 
%   nodes: Number of physical nodes 
%
%   layers: Number of layers for each aspect
%
%   transitionMatrix: matrix of copying probabilities
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

% Version: 
% Date: 
% Author: 
% Email: 

% set defaults for all options (note that we force 'IntermediateSteps' to 
% have its default value since we are not using intermediate results)
options=OptionStruct('UpdateSteps',100,'theta',1,'communities',10,'q',1,...
    'exponent',-3,'kmin',3,'kmax',50,'mu',0.1,'maxreject',100);

% set up options for sub-functions
PartitionOptions=OptionStruct('UpdateSteps','IntermediateSteps');
NullOptions=OptionStruct('theta','communities','q');
NetworkOptions=OptionStruct('exponent','kmin','kmax','mu','maxreject');

options.set(varargin) % error if invalid options specified

PartitionOptions.setvalid(options); % extract options
NullOptions.setvalid(options);
NetworkOptions.setvalid(options);

% generate partition using Dirichlet null distribution
S=PartitionGenerator(nodes,layers,transitionMatrix,...
    DirichletNullDistribution(layers,NullOptions),PartitionOptions);

% generate intralayer edges using DCSBM benchmark model
A=DCSBMNetworkGenerator(S,NetworkOptions);

end


