%% Non-Tensorial ESSTV for steady state parameter searching
% This script determines the parameters for a catalytic cracking reaction.
% The dynamic data set is first generated using MATLAB ODE solver, then fit
% using the parallel tempering algorithm.
%
% Dependencies: ParTemp.m in ../ directory
clear; clc
addpath('../')
addpath(genpath(pwd))

%% Determining dynamic data
% STEADY STATE
filename = 'SS.xlsx'; % Create an excel file with the first column 
% containing the shear rate and the second column containing shear stress
% for the steady state experimental data

DATA = xlsread(filename);

exp.TimeSS = [40; 30; 20; 16; 10; 6; 4; 2; 2; 2; 2; 2; 2; 2; 2; 2];
exp.ShearSS = DATA(:,1);
exp.StressSS = DATA(:,2);

% -------------------------------------------------------------------------
%% Fitting parameters
Mu0_c0=  0.0112415108323193; % Zero shear viscosity
MuINF_c0= 0.00404648370673836; % Infinite shear viscosity
TauC0= 0.0737203279463336; % RBC deformation time constant
tr10=  1.46599951771715;
tr20=  0.138452190122515;
MuR0=  0.0507902954111678; % Rouleaux viscosity
Sigy0= 0.00311483540501976; % Yield stress
TauLAM0=  0.418053358885443; % Overall rouleaux rebuild time constant
Gc0= 1.14959488638368; % RBC elastic modulus
GR0= 0.239756660459452; % Rouleaux elastic modulus

stochSOLN = [Mu0_c0,
    MuINF_c0,
    TauC0,
    tr10,
    tr20,
    MuR0,
    Sigy0];%,TauLAM0,Gc0,GR0];

% Minimization objective as a function of parameters
objective = @(parVec) SteadyState_OBJ(parVec, exp);

% Parallel tempering function call
qoptim = ParTemp(objective, ... % function to be minimized
                stochSOLN, ...  % Starting value (must be a row vector)
                1, ...          % EBHot (unity here)
                stochSOLN*0,'MIN', ... % Minimum value constraint
                stochSOLN*100,'MAX'); % Maximum value constraint


%% Plotting


[obj,pred] = SteadyState_OBJ(qoptim,exp);
save("CALCULATIONS\SS_partemp.mat","qoptim","pred","exp")


figure(1);
loglog(exp.ShearSS, exp.StressSS, 'ro',exp.ShearSS, pred.SigmaR,'b^', ...
    exp.ShearSS, pred.SigmaVISC,'g>',exp.ShearSS, pred.SigmaTOT,'k-', ...
    'MarkerSize',8,'LineWidth',2.);
set(gca,'FontSize',14,'FontWeight','bold','linewidth',2, 'FontName','Times');
legend('DATA','rouleaux stress','visc stress','tot stress','Location','NorthWest');
xlabel('Shear Rate (1/s)');
ylabel('Stress (Pa)');
xlim([min(exp.ShearSS) max(exp.ShearSS)]);
ylim([.0001 5]);

figure(2);
semilogx(exp.ShearSS, pred.LambdaSS, 'k-.','MarkerSize',8,'LineWidth',2.);
set(gca,'FontSize',14,'FontWeight','bold','linewidth',2, 'FontName','Times');
xlabel('Shear Rate (1/s)');
ylabel('Lambda');
xlim([min(exp.ShearSS) max(exp.ShearSS)]);
ylim([0 1]);