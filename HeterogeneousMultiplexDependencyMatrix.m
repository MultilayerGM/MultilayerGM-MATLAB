function P=HeterogeneousMultiplexDependencyMatrix(n_layers,n_blocks,p_in,p_out)
% Generate a heterogeneous multiplex dependency matrix
%
% Input: 
%   
%   n_layers: number of layers
%
%   n_blocks: number of blocks
%
%   p_in: probability for a state node to copy its community assignment from
%       corresponding state nodes in other layers in the same block
%
%   p_out: probability for a state node to copy its community assignment
%       from corresponding state nodes in other layers in different blocks 
%
% Output:
%
%   P: dependency matrix for use with PartitionGenerator
%
% Note that p_in+p_out<=1
%
%
% Version: 2.0.0
% Date: Thu 11 Jul 2019 15:24:15 CEST
% Author: Lucas Jeub
% Email: lucasjeub@gmail.com
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
%
block_size=floor(n_layers/n_blocks);
S=ones(n_layers,1)*n_blocks;
for i=1:n_blocks
    S((1:block_size)+(i-1)*block_size)=i;
end
P=zeros(n_layers,n_layers);
for i=1:n_blocks
    ind=S==i;
    P(ind,ind)=p_in/(sum(ind)-1);
    P(~ind,ind)=p_out/sum(~ind);
end
P=P-diag(diag(P));

end
