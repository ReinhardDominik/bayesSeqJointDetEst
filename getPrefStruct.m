function [ prefStruct ] = getPrefStruct( varargin )
%GETPREFSTRUCT returns a struct containing preferences, e.g. cvx solver, used by the methods

    prefStruct=struct();

    % get all available cvx solver
    [~,available_solver]=cvx_solver();
    
    checkCVXSolver = @(x) assert(any(strcmp(available_solver,x)),[x ' is not a installed cvx solver']);
    
    % define input parser
    p = inputParser;
    p.addParameter('cvx_solver','gurobi', checkCVXSolver); 
    p.addParameter('cvx_quiet',false,@islogical);
    
    % parse inputs
    parse(p,varargin{:})
    

    % set input parser fields to prefStruct
    fields=fieldnames(p.Results);    
    for i = 1:numel(fields)
        prefStruct.(fields{i}) = p.Results.(fields{i});
    end
    
end

