% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % This script reproduces the results of the shift-in-variance problem with low estimation constraint under H_1 presented in Section 7.2 and Appendix H of
% %     D. Reinhard, M. Fauss, and A. M. Zoubir, “Bayesian Sequential Joint Detection and Estimation,” 2018
% %
% % The script designs the test and runs a Monte Carlo simulation afterwards.
% %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%% cleanup
clearvars;
close all;


%% Parameters %%

% parameter space
thetaMin=0.01;
thetaMax=60;
Ntheta=9001;

% sample space
xMin = -20;
xMax = 20;
Nx = 6000;

% space of sufficient statistic
statMin = 0;
statMax = 25;
Nstat = 2100;

% constraints of the test
constr = [0.05 0.05 0.025 0.3];
Ntest = 100;

% parameters of the priors
% unfiform
a0=0.1;
b0=1;
% gamma
a1=1.7;
b1=0.5;
gammaShift = 1.3; % shift of the gamma distribution of H_1

% priors of the hypotheses
piVec=[0.5 0.5];

% Monte Carlo Setting
Nmc=1e7;

% misc
reg = 5e-5;

%% preprocessing
% building grids
thetaGrid=linspace(thetaMin,thetaMax,Ntheta);
dt=thetaGrid(2)-thetaGrid(1);
x=linspace(xMin,xMax,Nx);
stat=linspace(statMin,statMax,Nstat);

% discretized priors
prior=nan(Ntheta,2);
prior(:,1)=unifpdf(thetaGrid,a0,b0);
prior(:,2)=gampdf(thetaGrid - gammaShift ,a1,b1);
normConst=sum(prior,1)*dt;
if any(abs(1-normConst) > 0.01)
    warning('priors whould integrate up to one...normalizing them');
    disp(normConst)
    prior=prior./normConst;
end
clear normConst;

% function handles generating the priors
priorGen=cell(2,1);
priorGen{1} = @(N) a0 + (b0-a0)*rand(N,1);
priorGen{2} = @(N) gamrnd(a1,b1,N,1) + gammaShift;

% basename for storing results
st=dbstack; % getting stack trace
simulationName = st(end).name; % getting function/script name of the most high level function
clear st

testObj=varTest(Ntest, thetaGrid, prior, statMin,statMax, Nstat,x, piVec, constr,reg,priorGen);

testObj.designTest();

save(['./results/' simulationName '.mat'],'testObj');
fprintf('\n%s\n',['test saved in ./results/' simulationName '.mat']);

MCres = testObj.monteCarlo(Nmc,true,true);

save(['./results/' simulationName '_MC.mat'],'MCres');
fprintf('\n%s\n',['MC results saved in ./results/' simulationName '_MC.mat']);

evaluateMC(MCres);
