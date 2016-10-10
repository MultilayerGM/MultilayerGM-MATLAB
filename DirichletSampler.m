function weights=DirichletSampler(theta,q,nc)
% Sample weights for community null-distriubtion from a Dirichlet
% distribution
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

% Version: v1.0-alpha1
% Date: Mon 10 Oct 2016 16:12:35 EDT
% Author: Lucas Jeub
% Email: ljeub@iu.edu

weights=zeros(nc,1);
active=find(rand(nc,1)<=q);

weights(active)=gamrnd(theta,1,length(active),1);
weights=weights./sum(weights);

end

