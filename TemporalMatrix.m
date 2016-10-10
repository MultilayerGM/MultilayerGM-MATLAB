function L=TemporalMatrix(n_layers,p)
% Generate a uniform temporal transition matrix
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
%   L: transition matrix for use with PartitionGenerator
%
% Note that p<=1

% Version: 
% Date: 
% Author: 
% Email: 

L=zeros(n_layers);

for i=1:n_layers-1
    L(i,i+1)=p;
end

end
