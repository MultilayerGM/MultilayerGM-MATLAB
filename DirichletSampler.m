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
% Version: 2.0.0
% Date: Thu 11 Jul 2019 15:24:15 CEST
% Author: Lucas Jeub
% Email: lucasjeub@gmail.com
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
%       implemented in MATLAB," https://github.com/MultilayerGM/MultilayerGM-MATLAB (2016-2019).

weights=zeros(nc,1);
active=find(rand(nc,1)<=q);

weights(active)=gamrnd(theta,1,length(active),1);
weights=weights./sum(weights);

end
