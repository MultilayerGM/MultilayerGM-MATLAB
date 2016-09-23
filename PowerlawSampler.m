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


y=rand(n,1);
if t~=-1
    x=((x_max^(t+1)-x_min^(t+1))*y+x_min^(t+1)).^(1/(t+1));
else
    x=x_min*(x_max/x_min).^y;
end

end
