classdef generalTest1D < handle
    %GENERALTEST1D is the abstract class for a sequential joint detection and estimation test using a one dimensional test statistic
    %   This class provides properties needed for a one dimensional test in the joint detection and estimation case.
    %   An implementation of all methods which are independent of the likelihood and the random parameter are implemented directly.
    %   For all other methods exist as abstract methods.
    %
    %   For further information see:
    %   D. Reinhard, M. Fauss, and A. M. Zoubir, “Bayesian Sequential Joint Detection and Estimation,” Sequential Analysis, 2018. (submitted).
    %
    %
    %   generalTest1D Properties:
    %        Ntest      -    maximum number of samples
    %        thetaGrid  -    grid of the unkown parameter
    %        priorGrid  -    gridded version of the prior, priorGrid(:,i) is the prior under hypothesis H_i
    %        Nstat      -    number of grid points of the sufficient statistic
    %        maxStat    -    maximum value of the sufficient statistic (Ntest+1) x 1 vector
    %        minStat    -    minimum value of the sufficient statistic (Ntest+1) x 1 vector
    %        x          -    grid of the observation space
    %        piVec      -    vector containing a priori probabilities: P(H_0) and P(H_1)
    %        constr     -    4-element vector containing detection and estimation constraints
    %        reg        -    scalar regularization constant
    %        D          -    cost functions for continuing (obtained by LP)
    %        G          -    cost functions for stopping
    %        rho2       -    overall cost function (obtained by BI)
    %        D2         -    cost function for continuing (obtained by BI)
    %        rho        -    overall cost function (obtained by LP)
    %        l          -    Lagrange multiplier
    %        dN         -    regularization term
    %        EN         -    expected run-length
    %        priorGen   -    2-element cell array containing function handles to sample from prior
    %        decR       -    decision rule, same dimensions as cost function
    %        prefStruct -    struct containing preferences
    %
    %
    %   generalTest1D Methods:
    %       designTest              -   designs the test using linear programming followed by backward induction
    %       designTestLP            -   designs the test using linear programming
    %       designTestBI            -   designs the test using backward induction
    %       calcStoppingCostFcts    -   calculates the cost functions for stopping
    %       monteCarlo              -   performs a monte carlos simulation for an already designed test
    %       performTest             -   performs a single run of the monte carlo simulation
    %       performTest             -   performs an SPRT followed by an MMSE; single run
    %       calcLR                  -   calculates the likelihood ratio
    %       isDesigned              -   checks whether test is already designed
    %       hasFixedGrid            -   checks whether a fixed grid for the sufficient statistic is used for all n
    %       info                    -   displays some information about the test
    %       showRegions2D           -   visualizes the regions of the optimal test
    %       showRegions2DWald       -   visualizes the regions of the SPRT
    %       calcAllLLR              -   calculates the log-likelihood ratio for all values of the sufficient statistic
    %       plotPrior               -   plots the prior on theta
    %       plotPriors              -   plots the priors under both hypotheses
    %       calcStoppingCost        -   calculates the cost functions for stopping
    %       calcPostProbVar         -   calculates some posterior quantities
    %       getStatVec              -   returns the grid on the sufficient statistic at time n
    %

 
% ----------------------------------------------
%           Properties with protected SetAccess
% ----------------------------------------------
    properties(SetAccess=protected)
        Ntest;      % maximum number of samples
        thetaGrid;  % grid of the unkown parameter
        priorGrid;  % gridded version of the prior, priorGrid(:,i) is the prior under hypothesis H_i
        Nstat;      % number of grid points of the sufficient statistic
        maxStat;    % maximum value of the sufficient statistic (Ntest+1) x 1 vector
        minStat;    % minimum value of the sufficient statistic (Ntest+1) x 1 vector
        x;          % grid of the observation space
        piVec;      % vector containing a priori probabilities: P(H_0) and P(H_1)
        constr;     % 4-element vector containing detection and estimation constraints
        reg;        % scalar regularization constant

        % cost functions
        % for all cost functions: 
        %   first dimention: sufficient statistic, second dimension: time (starting at n=0)
        D;          % cost functions for continuing
        G;          % cost functions for stopping
        % cost functions calculated by backward induction
        rho2;
        D2;

        rho;        % overall cost function
        l;          % Lagrange multiplier
        dN;         % regularization term
        EN;         % expected run-length
        priorGen;   % 2-element cell array containing function handles to sample from prior
        decR;       % decision rule, same dimensions as cost function
        prefStruct; % struct containing preferences
    end

% ----------------------------------------------
%           Properties with protected GetAccess
% ----------------------------------------------    
    properties(GetAccess=protected)
       dt;      % distance on the grid of the parameter
       dx;      % distance on the grid of the observations
       dstat;   % distance on the grid of the sufficient statistic
       
       designed=false;   % flag whether test is already designed
       MCfinished=false; % flag whether the monte carlo simulation already finished
       
    end
    

% ----------------------------------------------
%           Implemented public non-static Methods
% ----------------------------------------------    
    methods
        

        % ---------------------------------------------
        %       Functions related to test design
        % ---------------------------------------------

        function designTest(obj)
        % DESIGNTEST runs the test design
        %   First, the test is designed via Linear Programming, then the cost functions are re-calculated by their recursive definition.
        %   Warnings are displayed when the two cost functions differ too much.
        %   Results are then stored in object properties.
        %
        %   Input:
        %       nothing (only the object itself)
        %
        %   Output:
        %       nothing

            obj.designTestLP();         % design test via linear programming
            obj.calcStoppingCostFcts(); % calculate the cost functions for stopping
            obj.designTestBI();         % re-calculate other cost functions by their recursive definition
            % calculate stopping regions for LP and BI cost functions
            M=obj.G>obj.D;
            M2=obj.G>obj.D2;
            % count percentage of elements where both regions differ
            percGridDiff=sum(M(:)~=M2(:))/numel(obj.G)*100;
            % calculate maximum difference of the rho functions
            maxRhoDiff=max(obj.rho2(:)-obj.rho(:));
            % warn if regions differ
            if percGridDiff > 0.01
               warning(['stopping region differs in ' sprintf('%.2f', (percGridDiff)) '% of the grid points...try to increase the regularization parameter']); 
            % warn if rhos differ too much
            elseif maxRhoDiff > 1
                warning(['LP result and BI result of the rho function differ up to ' sprintf('%.2f',maxRhoDiff) '...try to increase the regularization parameter']);
            end
            obj.designed=true; % now the test is designed

        end
                        
        function designTestLP(obj)
        % DESIGNTESTLP designs the test for the sequential joint detection and estimation problem via Linear Programming
        %
        %   Results are then stored in object properties.
        %   Input:
        %       nothing (only the object itself)
        %
        %   Output:
        %       nothing


            % getting vector of sufficient statistic for n=0
            stat0=obj.getStatVec(0);
            % getting target index for maximization
            id_t = round(numel(stat0)/2);

            %% cvx_block
            cvx_begin
                % selecting cvx_solver stored in prefStruct
                eval(sprintf('cvx_solver %s', obj.prefStruct.cvx_solver));
                cvx_quiet(obj.prefStruct.cvx_quiet);
                % declare cvx variables
                variable rho(obj.Nstat,obj.Ntest+1) nonnegative; % rho functions
                variable l(4,1)  nonnegative; % Lagrange multipliers
                expression D(obj.Nstat,obj.Ntest+1); % costs for continuing

                % objective
                maximize( rho(id_t,1) -obj.piVec(1)*l(1)*obj.constr(1) -obj.piVec(2)*l(2)*obj.constr(2) -obj.piVec(1)*l(3)*obj.constr(3) -obj.piVec(2)*l(4)*obj.constr(4)  +obj.reg*sum(rho*obj.dstat)/(obj.Ntest+1) ); %#ok<PROP,CPROP,*NODEF>
                fprintf('Constraints for n=');
                str='';
                % constraints for all times n=0,...,Ntest
                for n=0:obj.Ntest
                    fprintf(repmat('\b', [1, length(str)]));
                    str=sprintf('%i/%i',n,obj.Ntest);
                    fprintf('%s',str);
                    ii=n+1;

                    % % costs for stopping
                    [pH,v]=obj.calcPostProbVar(n);          % calculating posterior probabilities and variances
                    [g0,g1]=obj.calcStoppingCost(pH,v,l);   % calculate cost functions for stopping and deciding in favor of H_0 and H_1, respectively
                    rho(:,ii) <= g0; %#ok<VUNUS,PROP>       % constraint for H_0
                    rho(:,ii) <= g1; %#ok<VUNUS,PROP>       % constraint for H_1

                    % % costs for continuing
                    if n<obj.Ntest
                        % getting grid of sufficient statistic for current and next time instant
                        statGridCur = obj.getStatVec(n);
                        statGridNext = obj.getStatVec(n+1);

                        % look ahead
                        statNew = obj.statUpdate(statGridCur,n+1,obj.x.');                          % getting sufficient statistic for next step
                        idxStatNew = interp1(statGridNext,1:obj.Nstat,statNew,'nearest','extrap');  % mapping it back to the grid

                        % calculating posterior predictive
                        postPred = obj.calcPostPred(n);

                        % calculate cost for continuing and applying constraint
                        D(:,ii) = 1 + sum(reshape(rho(idxStatNew,ii+1),obj.Nstat,numel(obj.x)).*postPred,2)*obj.dx; %#ok<AGROW,PROP>
                        rho(:,ii) <= D(:,ii); %#ok<VUNUS,PROP>
                    end
                end
            cvx_end
            
            % % assign everything to object properties
            % overall cost function
            obj.rho=rho; %#ok<PROP>
            % cost for continuing at n=Ntest is set to infty
            D(:,end)=inf(size(D(:,end))); %#ok<PROP>
            % storing cost function for continuing
            obj.D=D; %#ok<PROP>
            % storing lagrange multiplier
            obj.l=l; %#ok<CPROP>

            % calculate & store regularization term
            obj.dN = obj.reg*sum(rho*obj.dstat)/(obj.Ntest+1); %#ok<PROP>
            % calculate & store expected run-length
            obj.EN=cvx_optval - obj.dN;

            % warn if ratio of regularization term and expected run-length is too large
            if obj.dN/obj.EN > 0.01
               warning('run-length differs more than 1% from the optimal value. Try to reduce the regularization parameter'); 
            end
        end
        
        
        function calcStoppingCostFcts(obj)
        % CALCSTOPPINGCOSTFCTS calculates the cost function for stopping and the stopping rule of the sequential test
        %   
        %   Results are then stored in object properties.
        %   Input:
        %       nothing (only the object itself)
        %
        %   Output:
        %       nothing

            
            % temporary variables for cost functions for stopping and the decision rule
            tmpG=nan(obj.Nstat,obj.Ntest+1);
            tmpDecR=nan(obj.Nstat,obj.Ntest+1);
            
            
            fprintf('Calculating cost functions for stopping ');
            str='';
            
            % iterating over time
            for n=0:obj.Ntest
                fprintf(repmat('\b', [1, length(str)]));
                str=sprintf('%i/%i',n,obj.Ntest);
                fprintf('%s',str);
                ii=n+1;
                [pH,v]=obj.calcPostProbVar(n);                  % calculating posterior probabilities and variances
                [g0,g1,g]=obj.calcStoppingCost(pH,v,obj.l);     % calculate cost functions for stopping and deciding in favor of H_0 and H_1, respectively
                tmpG(:,ii)=g;                                   % store cost function in temporary variable
                tmpDecR(:,ii)=g0>g1;                            % store decision rule in temporary variable
            end
            % assign both temporary variables to object properties
            obj.G=tmpG;
            obj.decR=tmpDecR;
            fprintf('\t...done!\n');
        end
       

        function designTestBI(obj)
        % DESIGNTESTBI designs the test for the sequential joint detection and estimation problem via Backward Indution
        %   For the design via Backward Induction, the cost coefficients have to be calculated before. If they are not calculated before, the methods ends with an error.
        %   The calculation of the cost functions for stopping is also required beforehand.
        %
        %   Results are then stored in object properties.
        %   Input:
        %       nothing (only the object itself)
        %
        %   Output:
        %       nothing

            % error of cost coefficients are not calculated yet.
            if isempty(obj.l)
                error('The cost coefficients are not calculated. Run designTestLP() before!');
            end
            % error of cost functions for stopping are not calculated yet.
            if isempty(obj.G)
                error('The cost functions for stopping are not calculated. Run calcStoppingCostFcts() before!');
            end

            %% temporary "output" variables
            tmpD=nan(obj.Nstat,obj.Ntest+1);        % cost functions for continuing
            tmpD(:,end)=inf(obj.Nstat,1);           % cost for continuing at n=Ntest is set to infty
            tmpRho=nan(obj.Nstat,obj.Ntest+1);      % overall cost function
            tmpRho(:,end)=obj.G(:,end);             % overall cost function at n=Ntest is set to cost for stopping at n=Ntest
            
            fprintf('Design test using backward induction ');
            str='';
            % calculating cost functions for n=Ntest-1,...,0
            for n=obj.Ntest-1:-1:0
                fprintf(repmat('\b', [1, length(str)]));
                str=sprintf('%i/%i',n,obj.Ntest);
                fprintf('%s',str);

                ii=n+1;
               
                % get grid of sufficient statistic for current and next time instant
                statGridCur = obj.getStatVec(n);
                statGridNext = obj.getStatVec(n+1);
                
                % look ahead step
                statNew = obj.statUpdate(statGridCur,n+1,obj.x.');                          % new value of the sufficient statustic
                idxStatNew = interp1(statGridNext,1:obj.Nstat,statNew,'nearest','extrap');  % map it on the grid
                
                % calculating posterior predictive
                postPred = obj.calcPostPred(n);

                % calculate cost function for continuing the test
                tmpD(:,ii) = 1 + sum(reshape(tmpRho(idxStatNew,ii+1),obj.Nstat,numel(obj.x)).*postPred,2)*obj.dx;
                % calculate overall cost function
                tmpRho(:,ii)=min(tmpD(:,ii),obj.G(:,ii));
               
            end

            % store temporary "output" variables in object properties
            obj.D2=tmpD;
            obj.rho2=tmpRho;
            
            fprintf('\t...done!\n');
        end
     

        % ---------------------------------------------------
        %       Functions related to Monte Carlo Simulation
        % ---------------------------------------------------
        
        function MCres = monteCarlo(obj,Nruns,varargin)
        % MONTECARLO performs a monte carlo simulation once the test is designed
        %   After the is designed, a Monte Carlo simulation can be performed to evaluate the empiral performance measures. If the test is not designed, the function will end with an error.
        %   For the evaluation of the cost functions, the cost functions calculated by Backward Induction are used. A linear Interpolation is used.
        %
        %   Input:
        %       obj         the object itself
        %       Nruns       The number of Monte Carlo runs
        %       varargin
        %           {1}     boolean whether to use the two-stage procedure for benchmarking, default is TRUE
        %           {2}     if set to true, the final sufficient statistic is also stored in output struct, default is false
        %
        %   Output:
        %       MCres       struct containing information about the Monte Carlo simulation
        %           HTrue       Array containing the true hypothesis; dimensions Nruns x 1
        %           thetaTrue   Array containing the true parameter values; dimensions Nruns x 1
        %           theteaHat   Array containing the estimated parameter values; dimensions Nruns x 1
        %           detection   Array containing the decisions; dimensions Nruns x 1
        %           tau         Array containing the run-length; dimensions Nruns x 1
        %           varPost     Array containing the posterior variances under the estimated hypothesis when the test stops; dimensions Nruns x 1
        %         only if varargin{1} is true
        %            thetaHat_W     estimated paramter using two-step procedure; dimensions Nruns x 1
        %            detection_W    decision using two-step procedure; dimensions Nruns x 1
        %            tau_W          run-length using two-step procedure; dimensions Nruns x 1
        %            varPost_W      posterior variance under the estimated hypothesis using two-step procedure; dimensions Nruns x 1
        %         only if varargin{2} is true
        %            stat           value of sufficient statistic when the optimal test stops; dimensions Nruns x 1

    
            % errors if something goes wrong
            if ~obj.designed
                error('Test is not designed yet. Please design the test first!');
            end
            if nargout == 0
              error('running a MC simulation without storing the output makes no sense'); 
            end

            % parsing the varargin
            if nargin == 2
                compare=true;
            else
                compare=varargin{1};
            end
            if nargin == 3
                statFlg=false;
            else
                statFlg=varargin{2};
            end
            
            % generating output struct
            MCres = struct();
            MCres.thetaTrue=nan(Nruns,1);
            MCres.thetaHat=nan(Nruns,1);
            MCres.detection = nan(Nruns,1);
            MCres.tau = nan(Nruns,1);
            MCres.varPost = nan(Nruns,1);
            if compare
               MCres.thetaHat_W=MCres.thetaHat; 
               MCres.detection_W = nan(Nruns,1);
               MCres.tau_W = nan(Nruns,1);
               MCres.varPost_W = nan(Nruns,1);
            end
            if statFlg
               MCres.stat=nan(Nruns,1); 
            end
            
            % generating random true hypothesis and true paramter values
            MCres.Htrue=rand(Nruns,1) > obj.piVec(1);
            MCres.thetaTrue(MCres.Htrue) = obj.priorGen{2}(sum(MCres.Htrue));
            MCres.thetaTrue(~MCres.Htrue) = obj.priorGen{1}(sum(~MCres.Htrue));


            fprintf('Starting Monte Carlo Simulation with %.2e runs\n', Nruns);
            str='';
            tStart=tic;
            for kk=1:Nruns
                % just to print the elapsed and remaining time
                if mod(kk,round(Nruns/1000)) == 0
                    tElapsed=toc(tStart);
                    perc=kk/Nruns*100;
                    tRemain = tElapsed*(100/perc - 1);
                    fprintf(repmat('\b', [1, length(str)-1]));
                    str = sprintf('Progress: %.2f%%%% - %s elapsed - approx %s remaining', perc,duration(0,0,tElapsed), duration(0,0,tRemain));
                    fprintf(str);
                end
               data=obj.generateData(MCres.thetaTrue(kk,1));            % generate the data
               [det, est, tau, varPost,stat] = obj.performTest(data);   % perform a single run
               % and store the results in the output struct
               MCres.detection(kk)=det;
               MCres.thetaHat(kk) = est;
               MCres.tau(kk) = tau;
               MCres.varPost(kk) = varPost;
           
               if compare   % run two-step procedure if enabled
                    [det_W, est_W, tau_W, varPost_W] = obj.performTestWald(data);  % one run
                    % storing the results
                    MCres.detection_W(kk)=det_W;
                    MCres.thetaHat_W(kk) = est_W;
                    MCres.tau_W(kk) = tau_W;
                    MCres.varPost_W(kk) = varPost_W;
                   
               end
               if statFlg % store value of the sufficient statistic if enabled
                  MCres.stat(kk)=stat; 
               end
            end
           fprintf('...done!\n');
        end
        

        function [det, est, tau, varPost,stat] = performTest(obj, data)
        %PERFORMTEST    performs a single run of the optimal test once the test is designed.
        %
        %   INPUT:
        %       obj         the object it self
        %       data        the data used for the test; dimensions Ntest x 1
        %
        %   OUTPUT:
        %       det         detection result
        %       est         parameter estimate
        %       tau         run-length
        %       varPost     posterior variance under detected hypothesis
        %       stat        sufficient statistic when the test stopped


            % errors if test is not designed yet
            if ~obj.designed
                error('Test is not designed yet. Please design the test first!');
            end

            stat=0; % initial value of the sufficient statistic
            for n=1:numel(data)
                ii=n+1;
                stat = obj.statUpdate(stat,n,data(n));  % update the sufficient statistic
                statVec = obj.getStatVec(n);            % getting vector of sufficient statistic for current time instant
                g=myinterp1(statVec,obj.G(:,ii),stat);  % evaluate cost function for stopping
                d=myinterp1(statVec,obj.D2(:,ii),stat); % evaluate cost function for continuing
                if g < d    % evaluate stopping rule
                   break; 
                end
            end
            tau=n;  % store run-length
            [pH, v, m] = obj.calcPostProbVar(tau,stat); % calculate posterior probabilities, posterior variance and posterior mean
            [g0,g1] = obj.calcStoppingCost(pH,v,obj.l); % calculate cost for deciding in favor of H_0 and H_1, respectively
            det = g0>g1;                                % evaluate decision rule
            est=m(:,det+1);                             % get estimate of the parameter
            varPost=v(:,det+1);                         % get posterior variance
      end
        

      function [det, est, tau, varPost] = performTestWald(obj, data)
      %PERFORMTESTWALD    performs a single run of the two-step procedure
        %
        %   INPUT:
        %       obj         the object it self
        %       data        the data used for the test; dimensions Ntest x 1
        %
        %   OUTPUT:
        %       det         detection result
        %       est         parameter estimate
        %       tau         run-length
        %       varPost     posterior variance under detected hypothesis
        %       stat        sufficient statistic when the test stopped


           % calculation of the Wald thresholds
           A=log((1-obj.constr(2))/(obj.constr(1)));
           B=log((obj.constr(2))/(1-obj.constr(1)));

           stat=0;  %inital value of the sufficient statistic
           for n=1:numel(data)
                stat = obj.statUpdate(stat,n,data(n));      % update the sufficient statistic
                lr=obj.calcLR(stat,n);                      % calculate the likelihood ratio
                
                llr=log(lr);                                % converting the likelihood ratio to the log domain
                
                if llr > A || llr < B                       % evaluate the stopping rule
                   break; 
                end
                
           end

            tau=n;          % assigning run-length to output
            if llr >= A     % decide in favor of H_1 if upper threshold crossed
                det = 1;
            elseif llr <= B % decide in favor of H_0 if lower threshold crossed
                det = 0;
            else            % if no decision made so far: decide in favor of H_0 if log likelihood ratio is positive
                det = llr > 0;
            end

            [~,v,m]=obj.calcPostProbVar(tau,stat);  % calculate posterior variance and posterior mean
            est=m(:,det+1);                         % store estimate in output variable
            varPost=v(:,det+1);                     % store posterior variance in output variable  
        end
        
        
        
        function LR = calcLR(obj,stat,n)
        %CALCLR calculates the likelihood ratio for a given sufficient statistic stat and time n
        %   This function calculates the likelihood ratio for given sufficient statistic and time.
        %   To calculate the LR, the posterior probabilities of the hypotheses are claculated first.
        %
        %   INPUT:
        %       stat        vector/scalar containing the sufficient statistic
        %       n           time index
        %
        %   OUTPUT:
        %       LR          likelihood ratio; scalar or vector depending on stat
        %


            Nstat=numel(stat); %#ok<PROPLC>         % getting number of elemts in sufficient statistic vector
            idx0=(obj.priorGrid(:,1)>0);            % getting indices for which theta corresponds to H_0
            idx1=(obj.priorGrid(:,2)>0);            % getting indices for which theta corresponds to H_1
            
            post=obj.calcPosterior(n,stat);         % calculate posterior density p(theta | stat)
            
            pH=nan(Nstat,2); %#ok<PROPLC>           % inizialize vector of posterior probabilities p(H_i | stat)
            pH(:,1) = sum(post(:,idx0),2)*obj.dt;   % calculate p(H_0 | stat)
            pH(:,2) = sum(post(:,idx1),2)*obj.dt;   % calculate p(H_1 | stat)
            
            % calculate likelihood ratio
            LR=obj.piVec(1)/obj.piVec(2).*pH(:,2)./pH(:,1);
        end
        
        function LLR = calcLLR(obj,stat,n)
        %CALCLLR calcujlates the log-likelihood ratio for a given vector stat and a given time n
        %
        %   INPUT:
        %       obje        the object itself
        %       stat        scalar or vector of sufficient statistic
        %       n           time instant
        %   
        %   OUTPUT:
        %       LLR         scalar or vector of log-likelihood ratios


            [pH] = obj.calcPostProbVar(n,stat);     % posterior probabilities of distributions
            LLR = log(obj.piVec(1)/obj.piVec(2).*pH(:,2)./pH(:,1)); % calculating the log-likelihood ratio
        end



        % ---------------------------------------------------
        %       auxiliary/misc functions
        % ---------------------------------------------------
      
        function out = isDesigned(obj)
        %ISDESIGNED returns true when the test is designed and false otherwise
            out = ~( isempty(obj.D) | isempty(obj.G) |isempty(obj.D2) |isempty(obj.l) );
        end
        
        function flg = hasFixedGrid(obj)
        %HASFIXEDGRID returns true when the test uses the same grid of the sufficient statistic for all n
           flg = all( obj.maxStat == obj.maxStat(1)) & all (obj.minStat == obj.minStat(1)); 
        end
        
        function info(obj)
        %INFO displays some information about the object
             fprintf('####grids###\n');
             fprintf('sample space sampled from %.2f to %.2f using %i grid points\n', min(obj.x), max(obj.x), numel(obj.x));
             if obj.hasFixedGrid
                 fprintf('sufficient statistic space sampled from %.2f to %.2f using %i grid points\n', obj.minStat(1), obj.maxStat(1), obj.Nstat);
             else
                 fprintf('sufficient statistic sampled using non constant grid with %i grid points\n', obj.Nstat);
             end
             fprintf('parameter space sampled from %.2f to %.2f using %i grid points\n', min(obj.thetaGrid), max(obj.thetaGrid), numel(obj.thetaGrid));
             
             fprintf('\n');
             fprintf('####constraints###\n');
             fprintf('target errors: %.2e %.2e %.2e %.2e\n', obj.constr);
             fprintf('maximum number of samples: %i\n', obj.Ntest);
             fprintf('priors of the two hypotheses: %.2f %.2f\n', obj.piVec);
             
             
             fprintf('\n');
             if ~obj.isDesigned()
                fprintf('test not designed yet\n'); 
             else
                fprintf('####design results###\n');
                fprintf('optimal cost coefficients: %.2f %.2f %.2f %.2f\n', obj.l);
                fprintf('expected run-length: %.2f\n', obj.EN);
                fprintf('regularization term: %.2e\n', obj.dN);
             end
             
         end

        % ---------------------------------------------------
        %       visualization functions
        % ---------------------------------------------------

         function finalMask= showRegions2D(obj)
            %SHOWREGIONS2D shows the different regions of the optimal test using a heatmap plot
            %   This functions plots the region in which the test has to continue, the region in which the test stops and decides in favor of H_0 and the region in which the test stops and decides in favor of H_1
            %   The x-axis of the plot denotes the time and the y-axis is the sufficient statistic. When the grid of the sufficient statistic is the not same for all n, the y-axis is index of the sufficient statistic vector at time n.
            %   If the test is not designed yet, the function ends with an error.
            %
            %   INPUT:
            %       just the object itself
            %
            %   OUTPUT:
            %       finalMask   a mask indicating the different regions of the test; dimensions Nstat x (Ntest+1)
            %                   0:  test has to continue
            %                   1:  test stops and decides in favor of H_1
            %                  -1:  test stops and decides in favor of H_1

            % error if test not designed
            if ~obj.isDesigned()
                error('test is not designed'); 
            end
            finalMask=nan(size(obj.G));                         % initialize output matrix
            stopMask=obj.G<obj.D2;                              % mask of the stopping region 
            finalMask(~stopMask)=0;
            finalMask(stopMask & obj.decR) = 1;                 % region in which the test stops and decides in favor of H_1
            finalMask(stopMask & ~obj.decR) = -1;               % region in which the test stops and decides in favor of H_0
            figure;                                             % creating the figure
            colormap(flag(3)); lcolorbar({'H0','cont','H1'});   % creating colormap and colorbar
            % y-data
            if obj.hasFixedGrid()                               % grid of sufficient statistic for fixed grid
               ydata=obj.getStatVec(0);
            else
                ydata=1:obj.Nstat;                              % index otherwise
            end

            imagesc('XData',0:obj.Ntest,'YData',ydata, 'CData', finalMask); %plotting the regions

            xlabel('n');                    % xlabel is the time
            xlim([0 obj.Ntest]);            % setting the correct x limits
            ylim([min(ydata), max(ydata)]); % same for y axis

            % label y axis as sufficient statistic or corresponding index
            if obj.hasFixedGrid()
                ylabel('stat');
            else
                ylabel('stat index');
            end
         end
        
        function finalMask=showRegions2DWald(obj)
           %SHOWREGIONS2DWALD shows the different regions of the Wald test using a heatmap plot
            %   This functions plots the region in which the test has to continue, the region in which the test stops and decides in favor of H_0 and the region in which the test stops and decides in favor of H_1
            %   The x-axis of the plot denotes the time and the y-axis is the sufficient statistic. When the grid of the sufficient statistic is the not same for all n, the y-axis is index of the sufficient statistic vector at time n.
            %
            %   INPUT:
            %       just the object itself
            %
            %   OUTPUT:
            %       finalMask   a mask indicating the different regions of the test; dimensions Nstat x (Ntest+1)
            %                   0:  test has to continue
            %                   1:  test stops and decides in favor of H_1
            %                  -1:  test stops and decides in favor of H_1

            % calculation of the Wald thresholds
            A=log((1-obj.constr(2))/(obj.constr(1)));
            B=log((obj.constr(2))/(1-obj.constr(1)));
            disp([A B]) % display the Wald thresholds

            finalMask=zeros(size(obj.G));   % initialization of the output matrix
            llr=obj.calcAllLLR();           % calculate the log likelihood ratio for all n and all values of the sufficient statistic
            finalMask(llr>A)=1;             % stop and decide in favor of H_1
            finalMask(llr<B)=-1;            % stop and decide in favor of H_0

            figure;                         % create the figure
            colormap(flag(3)); lcolorbar({'H0','cont','H1'});   % set colormap of colorbar labels
            % getting data for y-axis
            if obj.hasFixedGrid()
               ydata=obj.getStatVec(0);     % sufficient statistic for constant grid
            else
                ydata=1:obj.Nstat;          % index otherweise
            end

            imagesc('XData',0:obj.Ntest,'YData',ydata, 'CData', finalMask);     % show the regions

            xlabel('n');                        % time is the x-label
            xlim([0 obj.Ntest]);                % setting the limits of the x-axis
            ylim([min(ydata), max(ydata)]);     % and the one of the y-axis

            % label of y-axis
            if obj.hasFixedGrid()
                ylabel('stat');                 % statistic for constant grid
            else
                ylabel('stat index');           % index otherweise
            end
            
         end

        function LLRmtx = calcAllLLR(obj)
        %CALCALLLLR calculates the log-likelihood ratio for all n and all values of the sufficient statistic
        %
        %   INPUT:
        %       just the object itself
        %
        %   OUTPUT:
        %       LLRmtx  matrix containing the values of the log-likelihood ratio; dimensions Nstat x (Ntest+1)

            LLRmtx=nan(obj.Nstat,obj.Ntest);    % initialize output matrix
            for n=0:obj.Ntest                   % loopf over time
               disp(n);                         % display time index
               stat=obj.getStatVec(n);          % get stat vector of time n
               llr=log(obj.calcLR(stat,n));     % calculate the log-likelihood ratio
               LLRmtx(:,n+1)=(llr);             % assign to output
            end
        end
         
         
        function plotPrior(obj)
        %PLOTPRIOR plots the prior p(theta)
        %
        %   INPUT:
        %       just the object itself
        %   
        %   OUTPUT:
        %       nothing
        %
             figure; 
             plot(obj.thetaGrid,obj.priorGrid*obj.piVec);
             xlabel('\theta');
        end
         
        function plotPriors(obj)
        %PLOTPRIOR plots the two priors p(theta | H_0) and p(theta | H_1)
        %
        %   INPUT:
        %       just the object itself
        %   
        %   OUTPUT:
        %       nothing
        %
             figure; 
             plot(obj.thetaGrid,obj.priorGrid);
             xlabel('\theta');
             legend({'$p(\theta \mid H_0)$','$p(\theta\mid H_1)$'},'Interpreter','latex');
        end
         
         


         
    end
% ----------------------------------------------
%           Implemented static Methods
% ----------------------------------------------     
    methods (Static)
        function [g0,g1,g] = calcStoppingCost(pH, v,l)
        %calcSTOPPINGCOST calculates the cost for stopping and deciding in favor of H_0, the cost for stopping and deciding in favor of H_1 and the overall cost for stopping
        %   
        %   INPUT:
        %       pH      posterior probabilities of hypothesis; first dimension: sufficient statistic, second dimension: hypothesis
        %       v       posterior variances under both hypothesis; first dimension: sufficient statistic, second dimension: hypothesis
        %       l       cost coefficients, 4 element vector
        %
        %   OUTPUT:
        %       g0      cost for stopping and deciding in favor of H_0; same dimensions as pH(:,1)
        %       g1      cost for stopping and deciding in favor of H_1; same dimensions as pH(:,1)
        %       g       overall cost for stopping the test

            g0=l(2)*pH(:,2) + l(3)*pH(:,1).*v(:,1); % stopping and deciding in favor of H_0
            g1=l(1)*pH(:,1) + l(4)*pH(:,2).*v(:,2); % stopping and deciding in favor of H_1
            g=min(g0,g1);                           % overall cost for stopping
        end 
    end



% ----------------------------------------------
%           Abstract Methods
% ----------------------------------------------        
    methods (Abstract)
        statUpdate(stat,n,x)                % transition kernel of the sufficient statistic
        generateData(obj, theta)            % drawing samples when a specific theta is given
        calcPostPred(obj,n)                 % calculate the posterior predictive distribution
        calcLogPosterior(obj,n, varargin)   % calculate the logarithmic version of the posterior distribution p(theta | stat)
    end

% ----------------------------------------------
%           Private Methods
% ----------------------------------------------            
    methods(Access=protected)


        function [pH, v, m] = calcPostProbVar(obj,n,varargin)
        % CALCPOSTPROBVAR calculates the posterior probabilities of the hypothesis, the posterior variances and the posterior mean
        %
        %
        %   INPUT:
        %       obj         the object itself
        %       n           the time instant
        %       varargin
        %           {1}     scalar/vector containing values of the sufficient statistic which should be used; default: use whole grid
        %
        %   OUTPUT:
        %       pH          vector containing posterior probabilities of both hypotheses; dimensions Nstat x 2
        %       v           vector containing posterior variances under both hypotheses; dimensions Nstat x 2
        %       m           vector containing posterior m under both hypotheses; dimensions Nstat x 2

            % getting indices of theta grid for which H_0 and H_1 are true
            idx0=(obj.priorGrid(:,1)>0);
            idx1=(obj.priorGrid(:,2)>0);
            
            % gettin size of the stat vector
            Nstat=obj.Nstat;  %#ok<PROPLC>
            % calculate logarithmic version of the posterior distirbution
            if nargin == 3
                postLog=obj.calcLogPosterior(n,varargin{1});
                Nstat=numel(varargin{1});  %#ok<PROPLC>
            else
                postLog=obj.calcLogPosterior(n);
            end
            
            % initialize output variables
            pH=nan(Nstat,2); %#ok<PROPLC>
            m=nan(size(pH));
            v=nan(size(pH));


            % calculate the posterior probabilities of both hypotheses
            pH(:,1) = sum(exp(postLog(:,idx0)),2)*obj.dt;
            pH(:,2) = sum(exp(postLog(:,idx1)),2)*obj.dt;

            
            % calculating the posterior distribution under H_0
            post0L = (postLog);
            post0L(:,~idx0)=-inf;
            post0L = post0L - logOfSum_array(post0L,2) - log(obj.dt);
            post0 = exp(post0L);
    
            m(:,1) = sum(obj.thetaGrid.*post0,2)*obj.dt;    % calculate posterior mean under H_0
            v(:,1) = sum((bsxfun(@plus,obj.thetaGrid,-m(:,1))).^2.*post0,2)*obj.dt; %calculate posterior variance under H_0

            % calculating the posterior distribution under H_1
            post1L = (postLog);
            post1L(:,~idx1)=-inf;
            post1L = post1L - logOfSum_array(post1L,2) - log(obj.dt);
            post1 = exp(post1L);

            m(:,2) = sum(obj.thetaGrid.*post1,2)*obj.dt;    % calculate posterior mean under H_1
            v(:,2) = sum((bsxfun(@plus,obj.thetaGrid,-m(:,2))).^2.*post1,2)*obj.dt; % calculate posterior variance under H_1
        end

        
        function stat=getStatVec(obj,n)
            %GETSTATVEC returns the vector of the sufficient statistic for
            %time instant n
            %
            %   INPUT:
            %       obj     the object itself
            %       n       time instant
            %
            %   OUTPUT:
            %       stat    vector of sufficient statistic; vector containing Nstat elements
           stat=linspace(obj.minStat(n+1),obj.maxStat(n+1),obj.Nstat).';
        end
    end
    
end

