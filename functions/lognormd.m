function logdist = lognormd(x,dim)
% Normalizes an of logarithmic values into log-probabilities.

    logdist = bsxfun(@minus, x, logOfSum_array(x,dim));

end