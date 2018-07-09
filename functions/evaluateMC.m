function [detErr, estErr, tau, detErr_W, estErr_W, tau_W] = evaluateMC( MC, varargin )
%EVALUATEMC evaluates the Monte Carlo results
%   This function calculates the type-I and type-II error, the MSE under
%   both hypotheses, the average run-length under both hypotheses as well
%   as the overall average run-length. If the fields for the two-step
%   procedure exist, the performance metric for the two-step procedure is
%   calculated as well.
%
%   INPUT:
%       MC          struct containing Monte Carlo results
%                   See GENERALTEST1D.MONTECARLO
%       varargin
%           {1}     prints performance metric to stdout if true
%                   default: true
%
%   OUTPUT:
%       detErr      type-I and type-II error probabilities
%       estErr      MSE und H0 and H1
%       tau         array containing the average run-lenghts
%           1st element: ARL under H0
%           2nd element: ARL under H1
%           3rd element: overall ARL
%       only if two-step procedure is evaluated as well;
%       nan arrays otherwise:
%       detErr_W    type-I and type-II error probabilities
%       estErr_W    MSE und H0 and H1
%       tau_W       array containing the average run-lenghts
%           1st element: ARL under H0
%           2nd element: ARL under H1
%           3rd element: overall ARL

%   check verbose flag
    if nargin == 1
        verbose=true;
    else
        verbose = varargin{1};
    end

    % initialize output variables
    detErr=nan(2,1);
    estErr=nan(2,1);
    tau=nan(3,1);
    detErr_W=nan(2,1);
    estErr_W=nan(2,1);
    tau_W=nan(3,1);
    
    % calculate error probabilities
    detErr(1) = sum(MC.detection(MC.Htrue == 0) == 1 )/sum(MC.Htrue==0);
    detErr(2) = sum(MC.detection(MC.Htrue == 1) == 0 )/sum(MC.Htrue==1);
    
    % calculate MSE
    estErr(1) = mean( (MC.thetaTrue(MC.Htrue==0 & MC.detection == 0) - MC.thetaHat(MC.Htrue==0 & MC.detection == 0) ).^2 );
    estErr(2) = mean( (MC.thetaTrue(MC.Htrue==1 & MC.detection == 1) - MC.thetaHat(MC.Htrue==1 & MC.detection == 1) ).^2 );
    
    % calculate the ARLs
    tau(1) = mean(MC.tau( MC.Htrue == 0) );
    tau(2) = mean(MC.tau( MC.Htrue == 1) );
    tau(3) = mean(MC.tau);
    
    % print everything if enabled
    if verbose
        fprintf('\n\t%s\n\n','==Results for optimal test==');
        fprintf('type-I error: \t%.3f\n', detErr(1) );
        fprintf('type-II error: \t%.3f\n', detErr(2) );
        fprintf('MSE H0: \t%.3f\n', estErr(1) );
        fprintf('MSE H1: \t%.3f\n', estErr(2) );
        
        fprintf('ARL H0: \t%.3f\n', tau(1) );
        fprintf('ARL H1: \t%.3f\n', tau(2) );
        fprintf('ARL: \t\t%.3f\n', tau(3) );        
    end

    % check for fields of the two-step procedure
    if isfield(MC,'tau_W') && isfield(MC,'thetaHat_W') && isfield(MC,'detection_W') && isfield(MC,'varPost_W')
        
        % calculate error probabilities
        detErr_W(1) = sum(MC.detection_W(MC.Htrue == 0) == 1 )/sum(MC.Htrue==0);
        detErr_W(2) = sum(MC.detection_W(MC.Htrue == 1) == 0 )/sum(MC.Htrue==1);
        
        % calculate MSE
        estErr_W(1) = mean( (MC.thetaTrue(MC.Htrue==0 & MC.detection_W == 0) - MC.thetaHat_W(MC.Htrue==0 & MC.detection_W == 0) ).^2 );
        estErr_W(2) = mean( (MC.thetaTrue(MC.Htrue==1 & MC.detection_W == 1) - MC.thetaHat_W(MC.Htrue==1 & MC.detection_W == 1) ).^2 );
        
        % calculate ARLs
        tau_W(1) = mean(MC.tau_W( MC.Htrue == 0) );
        tau_W(2) = mean(MC.tau_W( MC.Htrue == 1) );
        tau_W(3) = mean(MC.tau_W);

        % print everything if enabled
        if verbose    
            fprintf('\n\t%s\n\n','==Results for two-step procedure==');
            fprintf('type-I error: \t%.3f\n', detErr_W(1) );
            fprintf('type-II error: \t%.3f\n', detErr_W(2) );
            fprintf('MSE H0: \t%.3f\n', estErr_W(1) );
            fprintf('MSE H1: \t%.3f\n', estErr_W(2) );

            fprintf('ARL H0: \t%.3f\n', tau_W(1) );
            fprintf('ARL H1: \t%.3f\n', tau_W(2) );
            fprintf('ARL: \t\t%.3f\n', tau_W(3) );        
        end
    end
    
end

