% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % This script reproduces the results of the shift-in-mean problem with relaxed estimation constraint under H_0 presented in Section 7.1 of
% %     D. Reinhard, M. Fauss, and A. M. Zoubir, “Bayesian Sequential Joint Detection and Estimation,” 2018
% %
% % The script designs the test and runs a Monte Carlo simulation afterwards.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


%% cleanup
clearvars;
close all;


%% Parameters %%

% parameter space
thetaMin=-12;
thetaMax=12;
Ntheta=4800;

% sample space
xMin=-15;
xMax=15;
Nx=6000;

% space of sufficient statistic
statMin=-8;
statMax=8;
Nstat=1600;

% parameter of likelihood
sigma=2;

% constraints of the test
constr=[0.05 0.025 0.35 0.9];
Ntest=100;

% parameters of the priors
a0=1.7;
b0=1;
a1=1.7;
b1=1;

% priors of the hypotheses
pi_vec=[0.5 0.5];

% Monte Carlo Setting
Nmc=1e7;

% misc
reg = 5e-4;

%% preprocessing
% building grids
thetaGrid=linspace(thetaMin,thetaMax,Ntheta);
dt=thetaGrid(2)-thetaGrid(1);
x=linspace(xMin,xMax,Nx);
stat=linspace(statMin,statMax,Nstat);

% discretized priors
prior=nan(Ntheta,2);
prior(:,1)=gampdf(-thetaGrid,a0,b0);
prior(:,2)=gampdf(thetaGrid,a1,b1);
normConst=sum(prior,1)*dt;
if any(abs(1-normConst) > 0.01)
    warning('priors whould integrate up to one...normalizing them');
    prior=prior./normConst;
end
clear normConst;

% function handles sampling from the priors
priorGen=cell(2,1);
priorGen{1} = @(N) -gamrnd(a0,b0,N,1);
priorGen{2} = @(N) gamrnd(a1,b1,N,1);

% basename for storing results
st=dbstack; % getting stack trace
simulationName = st(end).name; % getting function/script name of the most high level function
clear st

testObj=meanTest(Ntest, thetaGrid, prior, statMin, statMax, Nstat,x, pi_vec, constr,reg,priorGen, sigma.^2);

testObj.designTest();

save(['./results/' simulationName '.mat'],'testObj');
fprintf('\n%s\n',['test saved in ./results/' simulationName '.mat']);

MCres = testObj.monteCarlo(Nmc,true,true);

save(['./results/' simulationName '_MC.mat'],'MCres');
fprintf('\n%s\n',['MC results saved in ./results/' simulationName '_MC.mat']);

evaluateMC(MC);
