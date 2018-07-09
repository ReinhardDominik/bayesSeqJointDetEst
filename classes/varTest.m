classdef varTest < generalTest1D
    %VARTEST implements the shift-in-variance test based on the abstract class generalTest1D
    %   The constructor and all relevant distribution related functions are implemented.
    %   This class is used to generate the results presented in section 7.2 of 
    %
    %       D. Reinhard, M. Fauss, and A. M. Zoubir, “Bayesian Sequential Joint Detection and Estimation,” Sequential Analysis, 2018. (submitted).
    %
    %   Signal model:
    %
    %       H_0:    p(x | theta) ~ N(0, sigma_0^2)  mu_0 ~ p(sigma_0^2 | H_0)
    %       H_1:    p(x | theta) ~ N(0, sigma_0^2)  mu_1 ~ p(sigma_0^2 | H_1)
    %
    %   varTest Properties:
    %       llh         -   logarithmic version of the likelihood; normalization factor missing
    %       lh          -   likelihood; normalization factor missing
    %
    %   varTest Methods:
    %       varTest          -   constructor
    %       generateData     -   generates samples of the data
    %       calcLogPosterior -   calculates logarithmic version of the posterior density
    %       calcPosterior    -   calculates the posterior density
    %       calcPostPred     -   calculates the posterior predictive
    %       statUpdate       -   transition kernel of the sufficient statistic
    %       
    %       
    %   See also GENERALTEST1D
    
% ----------------------------------------------
%           Constant properties
% ----------------------------------------------

    properties(Constant)
        llh = @(sSq,sigmaSquare,n) -n*log((sqrt(2*pi*sigmaSquare))) -n*sSq./(2*sigmaSquare) ; % logarithmic version of the conditional distribution p( x | sigma^2) ; IMPORTANT: normalization constant is missing
        lh = @(sSq,sigmaSquare,n) (sqrt(2*pi*sigmaSquare)).^(-n).* exp(-n*sSq./(2*sigmaSquare)) ; %conditional distribution p( x | sigma^2) ; IMPORTANT: normalization constant is missing
    end

% ----------------------------------------------
%           Implemented public non-static Methods
% ----------------------------------------------    
    methods
        function obj=varTest(Ntest, thetaGrid, priorGrid, minStat, maxStat, Nstat,x, piVec, constr,reg,priorGen,varargin)
        %varTEST constructs an object of the class VARTest
        %   This constructor stores all input variables in properties of the object. A plausibility check for some of the input variables is done, before they are stored.
        %
        %   INPUT:
        %       Ntest       maximum number of samples used by the sequential scheme
        %       thetaGrid   grid on the paramter space; Ntheta number of elements
        %       priorGrid   discretized version of the prior on theta under both hypothesis; dimensions Ntheta x 2
        %       minStat     Vector containing the minimum values of the grid of the sufficient statistic for each n; repmat is applied if a scalar is passed
        %       maxStat     Vector containing the maximum values of the grid of the sufficient statistic for each n; repmat is applied if a scalar is passed
        %       Nstat       number of grid points used for the sufficient statistic
        %       x           grid on the observation space
        %       piVec       two-element vector containing prior probabilities of the two hypotheses
        %       constr      four-element vector containin the constraints on the error probabilities (first two elements) and on the estimation quality (last two elements)
        %       reg         non-negative, scalar regularization constant used by the linear programming
        %       priorGen    two element cell array containing function handles for sampling from prior distributions under H_0 and H_1, respectively
        %       varargin
        %           {1}     prefstruct, default GETPREFSTRUCT()
        %
        %   OUTPUT:
        %       obj         the constructed object

            % loading default pref struct if none is passed
            if nargin < 12
               obj.prefStruct=getPrefStruct();
            else
                obj.prefStruct = varargin{1};
            end

            obj.Ntest = Ntest;              % storing maximum number of samples
            obj.thetaGrid = thetaGrid(:).'; % reshaping and storing grid of parameter
            obj.dt = thetaGrid(2)-thetaGrid(1);   % aux variable, distance on parameter grid

            % check whether sizes of grid of theta and prior are equal
            if ~isequal(size(priorGrid),[2, numel(thetaGrid)]) && ~isequal(size(priorGrid),[numel(thetaGrid), 2])
               error('size of thetaGrid and corresponding prior do not match'); 
            end


            obj.priorGrid = reshape(priorGrid,numel(thetaGrid),2);   % storing discretized version of prior
            % store max value of statistic, repmat if needed
            if numel(maxStat) == 1
               maxStat=repmat(maxStat,Ntest+1,1); 
            end
            obj.maxStat = maxStat;
            % store MIN value of statistic, repmat if needed
            if numel(minStat) == 1
               minStat=repmat(minStat,Ntest+1,1); 
            end
            obj.minStat = minStat;

            obj.Nstat = Nstat;  % number of grid points for sufficient statistic
            obj.x = x(:);       % reshape and store grid on the observation space
            obj.dx = x(2)-x(1); % distiance on this grid

            % check if we have four constraints...
            if numel(constr) ~= 4
               error('provide a constraint vector which contains exactly 4 elements'); 
            end
            obj.constr=constr(:); % ... and store it
            
            % check range of regularization term
            if reg < 0
               error('regularization parameter must be positive'); 
            elseif reg > 1
               warning('regularization parameter may be too large');
            end
            obj.reg=reg;
            
            % piVec must have exactly two elements
            if numel(piVec) ~= 2
               error('piVec must have exactly two elements');
            end
            % reshape and store it
            obj.piVec = piVec(:);

            % normalize if necessary
            if sum(piVec) ~= 1
               warning('prior of hypotheses do not sum up to one...normalizing');
               obj.piVec = obj.piVec/sum(obj.piVec);
            end
            
            % store distance on grid of sufficient statistic for every n
            obj.dstat=zeros(obj.Ntest+1,1);
            for kk=1:obj.Ntest+1
                stat=linspace(obj.minStat(kk),obj.maxStat(kk),Nstat);
                obj.dstat(kk)=stat(2)-stat(1);
            end
            clear stat;

            % storing function handles for sampling from prior distribution
            obj.priorGen=priorGen;
                       
        end
        
        function data=generateData(obj,theta)
        %GENERATEDATA generate Ntest samples for a given mean theta
            data=randn(obj.Ntest,1)*sqrt(theta);
        end
        
        function postLog = calcLogPosterior(obj,n, varargin)
        %CALCLOGPOSTERIOR calculates the logarithmic version of the posterior distribution p(theta | stat)
        %
        %   INPUT:
        %       obj         the object itself
        %       n           time index
        %       varargin    
        %           {1}     scalar or vector containing sufficient statistic for which the posterior should be computed
        %                   default: grid at time n
        %
        %   OUTPUT:
        %       postLog     logarithmic version of the posterior; size Nstat x Ntheta
        %

            % getting sufficient statistic used later on
            if nargin == 3
                stat=varargin{1}; 
            else
                stat=obj.getStatVec(n);
            end
            lPrior = log(obj.priorGrid*obj.piVec).';   % logarithmic version of the prior p(H_i)*p(theta | H_i)
            
            postLog = bsxfun(@plus, obj.llh(stat,obj.thetaGrid,n), lPrior);% unnormalized version of the posterior
            
            normConst = logOfSum_array(postLog,2);  % getting normalization constant
            postLog = bsxfun(@minus,postLog,normConst) - log(obj.dt);   % and normalize
            
        end
        
        
        function post = calcPosterior(obj,n, stat)
        %CALCPOSTERIOR  calculates the posterior distribution p(theta | stat)
        %
        %   INPUT:
        %       obj     the object itself
        %       n       time instant
        %       stat    scalar or vector of the sufficient statistic which is used for calculating the posterior
        %
        %   OUTPUT:
        %       post    matrix containing the posterior distribution; size Nstat x Ntheta

            prior = (obj.priorGrid*obj.piVec).';        % prior p(H_i)*p(theta | H_i)
            post = bsxfun(@times, obj.lh(stat,obj.thetaGrid,n), prior); % getting unnormalized posterior
            normConst = sum(post,2);    % getting normalization constant
            post = post/(normConst*obj.dt);   % and normalizing
        end
        
        
        function postPred = calcPostPred(obj,n)
        %CALCPOSTPRED calculates the posterior predictive for a given time n
        %   
        %   INPUT:
        %       obj         the object itself
        %       n           time instant
        %
        %   OUTPUT:
        %       postPred    posterior predictive

            post = exp(obj.calcLogPosterior(n));    % calculate posterior

            % getting unnormalized likelihood                       
            lh = normpdf(obj.x,0,sqrt(obj.thetaGrid)).'; %#ok<PROPLC>
            % normalize it
            lh = lh./(sum(lh,2)*obj.dx);   %#ok<PROPLC>
            % getting unnormalized posterior predctive
            postPred = post*lh*obj.dt;%#ok<PROPLC>
            % und normalize posterior predictive
            postPred = bsxfun(@rdivide, postPred,(sum(postPred,2)*obj.dx));
            
        end   

        
    end
    
    methods (Static)
        function statNew = statUpdate(stat,n,x)
        %STATNEW is the transition kernel of the sufficient statistic from
        %time n to n+1 given statistic stat and new observation x
            statNew = bsxfun(@plus, (n-1)*stat/n, x.^2/n);
        end
    end
    
end

