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

% Version: 
% Date: 
% Author: 
% Email: 

weights=zeros(nc,1);
active=find(rand(nc,1)<=q);

weights(active)=gamrnd(theta,1,length(active),1);
weights=weights./sum(weights);

end

