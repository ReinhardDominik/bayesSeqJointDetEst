function linInds = dim2lin(indices, siz, dim)
% Given an array MAT of size SIZ=[n_1 n_2 ... n_d], and an index array
% INDICES of dimension [n_1 n_2 ... n_dim-1 ... n_d] which "sits on top" of
% the array MAT and indices elements along dimension DIM, the function
% computes the corresponding linear indices of these elements.

sLow    = prod(siz(1:dim-1));
sHigh   = prod(siz(dim+1:end));
linInds = bsxfun(@plus, (1:sLow)', (0:sHigh-1)*sLow*siz(dim));
linInds = linInds(:) + (indices(:)-1)*sLow;

end