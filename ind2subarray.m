function [sub] = ind2subarray(siz,ndx)
%IND2SUBARRAY Subscripts from linear index. Same as IND2SUB but
%   returns results as an array instead of separate output arguments.

siz = double(siz);
lensiz = length(siz);
sub=zeros(length(ndx),length(siz));

k = cumprod(siz);
if ndx>k(end)
    error('MultilayerBenchmark:ind2subarray:IndexOutOfRange','Index out of range.');
end
for i = lensiz:-1:2,
    vi = rem(ndx-1, k(i-1)) + 1;
    vj = (ndx - vi)/k(i-1) + 1;
    sub(:,i) = double(vj);
    ndx = vi;
end
sub(:,1)=ndx;

end

