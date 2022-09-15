%% t-ESSTV
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

clear DATA
%% STEPS DATA:

filename = 'Steps.xlsx'; % Creat an excel file that contains 6 tabs in the 
% following order: ["Down 5.0 -> 0.1","Down 10.0 -> 0.1","Down 20.0 -> 0.1",
                  % "Up 0.1 -> 5.0","Up 0.1 -> 10.0","Up 0.1 -> 20.0"]

names = ["Down_5p0_0p1","Down_10p0_0p1","Down_20p0_0p1", ...
    "Up_0p1_5p0","Up_0p1_10p0","Up_0p1_20p0"];

for i=1:length(names)
    DATA{i} = xlsread(filename,i);
end

exp.tD1=DATA{1}(:,1);  exp.sigD1=DATA{1}(:,2); exp.gamD1=DATA{1}(:,4);
exp.tD2=DATA{2}(:,1);  exp.sigD2=DATA{2}(:,2); exp.gamD2=DATA{2}(:,4);
exp.tD3=DATA{3}(:,1);  exp.sigD3=DATA{3}(:,2); exp.gamD3=DATA{3}(:,4);

exp.tU1=DATA{4}(:,1);  exp.sigU1=DATA{4}(:,2); exp.gamU1=DATA{4}(:,4);

% exp.tU2=DATA{5}(:,1);  exp.sigU2=DATA{5}(:,2)-.005; exp.gamU2=DATA{5}(:,4);
% exp.tU3=DATA{6}(:,1);  exp.sigU3=DATA{6}(:,2)-.01; exp.gamU3=DATA{6}(:,4);

exp.tU2=DATA{5}(:,1);  exp.sigU2=DATA{5}(:,2); exp.gamU2=DATA{5}(:,4);
exp.tU3=DATA{6}(:,1);  exp.sigU3=DATA{6}(:,2); exp.gamU3=DATA{6}(:,4);

% Step experiments (in order)
exp.gam_a = [0.1,0.1,0.1,5,10,20]; % Initial shear rates (1/s)
exp.gam_b = [5,10,20,0.1,0.1,0.1]; % Target shear rate (1/s)

% Times, stresses, shear, and norms in proper order as is written above
exp.times = {exp.tU1,exp.tU2,exp.tU3,exp.tD1,exp.tD2,exp.tD3};
exp.shears = {exp.gamU1,exp.gamU2,exp.gamU3,exp.gamD1,exp.gamD2,exp.gamD3};
exp.stresses = {exp.sigU1,exp.sigU2,exp.sigU3,exp.sigD1,exp.sigD2,exp.sigD3};
exp.NormVal = [max(exp.sigU1),max(exp.sigU2),max(exp.sigU3), ...
    max(exp.sigD1),max(exp.sigD2),max(exp.sigD3)];

clear DATA

% -------------------------------------------------------------------------
%% Fitting parameters

Mu0_c=  0.0112415108323193; % Zero shear viscosity
MuINF_c= 0.00404648370673836; % Infinite shear viscosity
TauC= 0.0737203279463336; % RBC deformation time constant
Sigy= 0.00311483540501976; % Yield stress
tr1=  1.46599951771715;
tr2=  0.138452190122515;
MuR=  0.0507902954111678; % Rouleaux viscosity
TauLAM=  1.18053358885443; % Overall rouleaux rebuild time constant
Gc= 1.14959488638368; % RBC elastic modulus
GR= 0.239756660459452; % Rouleaux elastic modulus

SSqoptim = reshape(cell2mat(struct2cell(load("CALCULATIONS\SSqoptim.mat"))),[1,7]);
start = reshape(cell2mat(struct2cell(load("CALCULATIONS\tESSTV_Stochastic.mat").best_par)),[1,10]);


% Minimization objective as a function of parameters
objective = @(parVec) tESSTV_OBJ(parVec, exp);

% Parallel tempering function call
qoptim = ParTemp(objective, ... % function to be minimized
                start, ...  % Starting value (must be a row vector)
                1, ...          % EBHot (unity here)
                start*0.001,'MIN', ... % Minimum value constraint
                start*10^3,'MAX'); % Maximum value constraint

%% Save variables

save("CALCULATIONS\tESSTV_partemp.mat","qoptim")
save("DATA\tESSTV_partemp.mat","exp")

%% Plotting
calc = load("CALCULATIONS\tESSTV_partemp.mat");
DATA = load("DATA\tESSTV_partemp.mat");

qoptim = calc.qoptim;
exp = DATA.exp;

[obj,pred] = tESSTV_OBJ(qoptim,exp);

stressTOT = pred.stressTOT;
stressElas = pred.stressElas;
stressVisc = pred.stressVisc;

figure(1);
loglog(exp.ShearSS, exp.StressSS, ...
    'ro',exp.ShearSS, stressTOT,'k-.', ...
    exp.ShearSS, stressVisc,'b-.', ...
    exp.ShearSS, stressElas,'g-.', ...
    'MarkerSize',8,'LineWidth',2.);
set(gca,'FontSize',14,'FontWeight','bold','linewidth',2, 'FontName','Times');
legend('DATA', ...
    'tot stress', ...
    'visc stress', ...
    'rouleaux stress', ...
    'Location', ...
    'NorthWest');
xlabel('Shear Rate (1/s)');
ylabel('Stress (Pa)');
xlim([min(exp.ShearSS) max(exp.ShearSS)]);
ylim([.0001 5]);

figure(2);
semilogx(exp.ShearSS, pred.lambdaSS, 'k-.','MarkerSize',8,'LineWidth',2.);
set(gca,'FontSize',14,'FontWeight','bold','linewidth',2, 'FontName','Times');
xlabel('Shear Rate (1/s)');
ylabel('Lambda');
xlim([min(exp.ShearSS) max(exp.ShearSS)]);
ylim([0 1]);

figure(3);
semilogx(exp.tU1, exp.sigU1,'ro', ...
    exp.tU2, exp.sigU2,'bd', ...
    exp.tU3, exp.sigU3,'g*', ...
    pred.t{1}, pred.s{1},'r-.', ...
    pred.t{2}, pred.s{2},'b-.', ...
    pred.t{3}, pred.s{3},'g-.', ...
    'MarkerSize',6,'LineWidth',2);
set(gca,'FontSize',14,'FontWeight','bold','linewidth',2, 'FontName','Times');
xlabel('Time (s)');
ylabel(' Stress (Pa)');
legend('5(1/s)', ...
    '10(1/s)', ...
    '20(1/s)');
title('Step Ups from 0.1(1/s) to');

figure(4);
semilogx( pred.t{1}, pred.lam{1},'r-.', ...
    pred.t{2}, pred.lam{2},'b-.', ...
    pred.t{3}, pred.lam{3},'g-.', ...
    'MarkerSize',6,'LineWidth',2);
set(gca,'FontSize',14,'FontWeight','bold','linewidth',2, 'FontName','Times');
xlabel('Time (s)');
ylabel(' Lambda');
legend('5(1/s)', ...
    '10(1/s)', ...
    '20(1/s)');
title('Step Ups from 0.25(1/s) to');

figure(5);
semilogx(exp.tD1, exp.sigD1,'ro', ...
    exp.tD2, exp.sigD2,'bd', ...
    exp.tD3, exp.sigD3,'g*', ...
    pred.t{4}, pred.s{4},'r-.', ...
    pred.t{5}, pred.s{5},'b-.', ...
    pred.t{6}, pred.s{6},'g-.', ...
    'MarkerSize',6,'LineWidth',2);
set(gca,'FontSize',14,'FontWeight','bold','linewidth',2, 'FontName','Times');
xlabel('Time (s)');
ylabel(' Stress (Pa)');
legend('5(1/s)', ...
    '10(1/s)', ...
    '20(1/s)');
title('Step Downs from 5.0(1/s) to');

figure(6);
semilogx( pred.t{4}, pred.lam{4},'r-.', ...
    pred.t{5}, pred.lam{5}, 'b-.', ...
    pred.t{6}, pred.lam{6},'g-.', ...
    'MarkerSize',6,'LineWidth',2);
set(gca,'FontSize',14,'FontWeight','bold','linewidth',2, 'FontName','Times');
xlabel('Time (s)');
ylabel(' Lambda');

figure(7);
title('Step up from 0.1 to 5s-1');
semilogx(pred.t{3}, pred.lam{3}, ...
    pred.t{3},pred.Strain_es{3}, ...
    pred.t{3}, pred.Shear_ps{3},'MarkerSize',6,'LineWidth',2);
set(gca,'FontSize',14,'FontWeight','bold','linewidth',2, 'FontName','Times');
legend('lambda', ...
    'GammaE', ...
    'ShearP');
xlabel('Time (s)');

figure(8);
title('Step down from 5 to 0.1s-1');
semilogx(pred.t{6}, pred.lam{6}, ...
    pred.t{6}, pred.Strain_es{6}, ...
    pred.t{6}, pred.Shear_ps{6},'MarkerSize',6,'LineWidth',2);
set(gca,'FontSize',14,'FontWeight','bold','linewidth',2, 'FontName','Times');
legend('lambda', ...
    'GammaE', ...
    'ShearP');
xlabel('Time (s)');