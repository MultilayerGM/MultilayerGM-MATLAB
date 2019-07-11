function P=MultiplexDependencyMatrix(n_layers,p)
% Generate a uniform multiplex dependency matrix
%
% Input:
%
%   n_layers: number of layers
%
%   p: probability for a state node to copy its community assignment from
%       corresponding state nodes in other layers
%
% Output:
%
%   P: dependency matrix for use with PartitionGenerator
%
% Note that p<=1
%
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
%

P(1:n_layers,1:n_layers)=p/(n_layers-1);
P=P-diag(diag(P));

end
