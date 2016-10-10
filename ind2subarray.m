function [sub] = ind2subarray(shape,ind)
%IND2SUBARRAY Subscripts from linear index. Same as IND2SUB but
%   returns results as an array instead of separate output arguments.

% Version: v1.0-alpha1
% Date: Mon 10 Oct 2016 16:12:35 EDT
% Author: Lucas Jeub
% Email: ljeub@iu.edu

shape = double(shape);
dims = length(shape);
sub=zeros(length(ind),length(shape));

k = cumprod(shape);
if ind>k(end)
    error('MultilayerBenchmark:ind2subarray:IndexOutOfRange','Index out of range.');
end
for i = dims:-1:2,
    subi = rem(ind-1, k(i-1)) + 1;
    subj = (ind - subi)/k(i-1) + 1;
    sub(:,i) = double(subj);
    ind = subi;
end
sub(:,1)=ind;

end

