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
% Version: v1.0-beta1
% Date: Thu 24 Nov 2016 15:13:51 EST
% Author: Lucas Jeub
% Email: jeub@maths.ox.ac.uk
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
%       implemented in MATLAB," [insert website] (2016).
%

P(1:n_layers,1:n_layers)=p/(n_layers-1);
P=P-diag(diag(P));

end
