function P=TemporalCPDependencyMatrix(n_layers,p, pc, nc)
% Generate a uniform temporal dependency matrix
%
% Input:
%
%   n_layers: number of layers
%
%   p: probability for a state node to copy its community assignment from
%       corresponding state node in previous layer when not at change point
%
%   pc: probability for a state node to copy its community assignment from
%       corresponding state node in previous layer when at change point
%
%   nc: number of changepoints
%
% Output:
%
%   P: dependency matrix for use with PartitionGenerator
%
% Note that p,pc<=1
%
% Version: 2.0.0
% Date: Thu 11 Jul 2019 15:24:15 CEST
% Author: Lucas Jeub
% Email: lucasjeub@gmail.com
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

if p>1||p<0
    error('MultilayerBenchmark:TemporalDependencies:p',...
        'Copying probability p out of range')
end

P=zeros(n_layers);
cp_skip=ceil(n_layers/(nc+1));

for i=1:n_layers-1
    P(i,i+1)=p;
end

for i=cp_skip:cp_skip:n_layers-1
    P(i,i+1)=pc;
end
end
