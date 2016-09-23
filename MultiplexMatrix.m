function L=MultiplexMatrix(n_layers,p)
% Generate a uniform multiplex transition matrix
%
% Input: 
%   
%   n_layers: number of layers
%
%   p: probability for a state node to copy its community assignment from
%       corresponding state nodes
%
% Output:
%
%   L: transition matrix for use with PartitionGenerator
%
% Note that p<=1

L(1:n_layers,1:n_layers)=p/(n_layers-1);
L=L-diag(diag(L));

end
