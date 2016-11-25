function P=TemporalDependencyMatrix(n_layers,p)
% Generate a uniform temporal dependency matrix
%
% Input: 
%   
%   n_layers: number of layers
%
%   p: probability for a state node to copy its community assignment from
%       corresponding state node in previous layer
%
% Output:
%
%   P: dependency matrix for use with PartitionGenerator
%
% Note that p<=1
%
% Version: v1.0
% Date: Fri 25 Nov 2016 16:24:28 EST
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
%       implemented in MATLAB," [insert website] (2016).

if p>1||p<0
    error('MultilayerBenchmark:TemporalDependencies:p',...
        'Copying probability p out of range')
end

P=zeros(n_layers);

for i=1:n_layers-1
    P(i,i+1)=p;
end

end
