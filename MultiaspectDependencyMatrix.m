function P=MultiaspectDependencyMatrix(layers, p, type)
% Generate a block-multiplex dependency matrix
%
% Input:
%
%   layers: array such that layers(i) is the number of layers in the ith aspect
%
%   p: vector of probabilities where p(i) is the probability for a state node to
%      copy its community assignment from corresponding state nodes in the ith
%      aspect
%
%   type: string specifying the type of dependency for each aspect. The
%         possible types are 'c' (or 'm') for categorical (multiplex) dependency
%         and 'o' (or 't') for ordinal (temporal) dependency.
%
% Output:
%
%   P: Dependency matrix for use with PartitionGenerator
%
% Note that sum(p)<=1
%
% Version: 2.0.0
% Date: Thu 11 Jul 2019 15:17:48 CEST
% Author: Lucas Jeub
% Email: ljeub@iu.edu
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

  if sum(p)>1
      error('MultilayerBenchmark:MultiplexTemporalMatrix:probability',...
          'Sum of copying probabilities>1');
  end

  P = zeros(prod(layers));
  for a = 1:numel(layers)
      P = P + kron(kron(eye(prod(layers(a+1:end))),coupling(layers(a), type(a), p(a))),...
                   eye(prod(layers(1:a-1))));
  end

  end

  function C = coupling(size, type, p)
      switch type
          case {'m', 'c'}
              C = ones(size)-eye(size);
              C = C .* p/(size-1);
          case {'t', 'o'}
              C = diag(ones(size-1,1),1);
              C = C * p;
          otherwise
              error('unknown coupling type %s', type)
      end
  end
