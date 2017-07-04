function x=PowerlawSampler(n,t,x_min,x_max)
% Sample values from a truncated powerlaw distribution
%
% Input: 
%
%   n: number of samples
%
%   t: exponent
%
%   x_min: minimum cut-off
%
%   x_max: maximum cut-off
%
% Output:
%
%   x: sampled values
%
% Version: 1.0.1
% Date: Tue  4 Jul 2017 16:38:06 BST
% Author: Lucas Jeub
% Email: ljeub@iu.edu
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
%       implemented in MATLAB," https://github.com/MultilayerBenchmark/MultilayerBenchmark (2016).


y=rand(n,1);
if t~=-1
    x=((x_max^(t+1)-x_min^(t+1))*y+x_min^(t+1)).^(1/(t+1));
else
    x=x_min*(x_max/x_min).^y;
end

end
