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


% Version: 
% Date: 
% Author: 
% Email: 

P(1:n_layers,1:n_layers)=p/(n_layers-1);
P=P-diag(diag(P));

end
