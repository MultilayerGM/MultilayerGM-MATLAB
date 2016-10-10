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


% Version: v1.0-alpha1
% Date: Mon 10 Oct 2016 16:12:35 EDT
% Author: Lucas Jeub
% Email: ljeub@iu.edu

P(1:n_layers,1:n_layers)=p/(n_layers-1);
P=P-diag(diag(P));

end
