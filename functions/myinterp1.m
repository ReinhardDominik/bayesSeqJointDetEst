function interpVal = myinterp1(x,v,xq)
%MYINTERP1 is a custom interpolation function. It uses linear interpolation
%if the query point is inside the grid and nearest interpolation otherwise
%   INPUT:
%       x       grid on x-axis
%       v       function value on grid
%       xq      query points
%
%   OUTPUT:
%   interpVal   vector of the same size as xq containing the interpolated values

    % check dimensions of x and v
    if sum(size(x) ~= size(v)) & sum(size(x) ~= size(v.'))
        error('dimension mismatch');
    end
    
    % initialize output variable
    interpVal = nan(size(xq));
    % get indices of elements outside the grid
    idxOutside=xq > max(x) || xq < min(x);
    % use nearest interpolation for values outside the grid...
    interpVal(idxOutside) = interp1(x,v,xq(idxOutside),'nearest','extrap'); 
    % ... and linear interpolation for the other ones
    interpVal(~idxOutside) = interp1(x,v,xq(~idxOutside),'linear','extrap');    

end
