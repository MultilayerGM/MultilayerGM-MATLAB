function L=BlockMultiplexMatrix(n_blocks,n_layers,p_in,p_out)
% Generate a block-multiplex transition matrix
%
% Input: 
%   
%   n_blocks: number of equal-sized blocks of layers
%
%   n_layers: number of layers in each block
%
%   p_in: probability for a state node to copy its community assignment from
%       corresponding state nodes in the same block
%
%   p_out: [default: 0] probability for a state node to copy its community 
%       assignment from corresponding state nodes in a different block 
%
% Output:
%
%   L: transition matrix for use with PartitionGenerator
%
% Note that p_in+p_out<=1

% Version: 
% Date: 
% Author: 
% Email: 

if nargin<4
    p_out=0;
end

if p_in+p_out>1
    error('MultilayerBenchmark:BlockMultiplexMatrix:probability',...
        'Sum of copying probabilities>1');
end

A=(ones(n_layers)-eye(n_layers))*p_in/(n_layers-1);
B=ones(n_layers)*p_out/(n_layers*(n_blocks-1));

L=kron(eye(n_blocks),A)+kron(ones(n_blocks)-eye(n_blocks),B);

end
