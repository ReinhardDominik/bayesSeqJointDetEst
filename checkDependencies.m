function [ out ] = checkDependencies()
%CHECKDEPENDENCIES checks whether needed toolboxes and libraries are installed and available
%   See README.md for further details

    str='';
    
    if ~exist('cvx','file') %checking for cvx
        str=sprintf('%s\n%s',str,'cvx not installed!');
    end

    if license('test','signal_toolbox') ~= 1 %checking for signal provessing toolbox
        str=sprintf('%s\n%s',str,'signal processing toolbox is missing!');
    end
    
    if license('test','statistics_toolbox') ~= 1 % checking for statistics toolbox
        str=sprintf('%s\n%s',str,'statistics toolbox is missing!');
    end
    
    % display message
    if ~isempty(str)
       error('Dependencies are not fulfilled \n%s', str);
    end

    
    out=true;
    fprintf('everything is ok!\n');
end

