function ndx = subarray2ind(shape,suba)
%SUBARRAY2IND Linear index from multiple subscripts.
%   Same as SUB2IND but takes array as input rather than multiple arguments.

% Version: 
% Date: 
% Author: 
% Email: 

shape = double(shape);
dims = length(shape);
if size(suba,2)~=dims
   error('MultilayerBenchmark:subarray2ind:DimensionMismatch','The subscript dimension does not match the array size.');
end


if any(min(suba(:,1)) < 1) || any(max(suba(:,1)) > shape(1))
    %Verify subscripts are within range
    error('MultilayerBenchmark:subarray2ind:IndexOutOfRange','Subscript out of range.');
end

ndx = double(suba(:,1));
    %Compute linear indices
    k = cumprod(shape);
    for i = 2:dims
        %%Input checking
        if (any(min(suba(:,i)) < 1)) || (any(max(suba(:,i)) > shape(i)))
            %Verify subscripts are within range
            error('MultilayerBenchmark:subarray2ind:IndexOutOfRange','Subscript out of range.');
        end
        ndx = ndx + (double(suba(:,i))-1)*k(i-1);
    end
end
