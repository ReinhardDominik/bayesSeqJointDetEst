function logSum = logOfSum_array(logVals, dim)
% Numerically stable computation of the logarithm of a sum, where all
% terms are represented in the log-domain.
%
% logVals: logarithmic representation of the values to be added
%
% dim: dimension along the summation should be performed
%
% logSum: logarithmic representation of the sum array


%%% First, permute the array such that we can perform all operations along
%%% the first dimension. This is because we will read out specific
%%% elements of the array in a later step where we must guarantee that one
%%% element is removed in each sub-array along the first dimension.
permInd         = 1:ndims(logVals);
permInd(1)      = dim;
permInd(dim)    = 1;
logValsPerm     = permute(logVals, permInd);

%%% Size of the array and reduced size when we remove the max-elements %%%
s               = size(logValsPerm);
s_reduced       = s;
s_reduced(1)    = s(1)-1;

%%% Locate the max-elements and get their linear indices %%%
[maxVals, inds] = max(logValsPerm, [], 1);
linInds         = dim2lin(inds, s, 1);
mask            = true(numel(logValsPerm),1);
mask(linInds)   = false;

%%% Remove the max-elements %%%
others      = reshape(logValsPerm(mask), s_reduced);

%%% Compute the output value using the permuted arrays %%%
logSumPerm  = maxVals + log1p(sum(exp(bsxfun(@minus, others, maxVals)), 1));

%%% Undo the permutation %%%
logSum = ipermute(logSumPerm, permInd);

end