function weights=DirichletSampler(theta,q,nc)
% Sample weights for community null-distriubtion from a Dirichlet distribution
%
% Input: 
%
%   theta: shape parameter for Dirichlet distribution
%
%   q: activation probability (probability for a given community
%              weight to be non-zero)
%
%   nc: number of communities
%
% Output: 
%
%   weights: sampled weigths (i.e. probabilities of a categorical
%       distribution)
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

weights=zeros(nc,1);
active=find(rand(nc,1)<=q);

weights(active)=gamrnd(theta,1,length(active),1);
weights=weights./sum(weights);

end

