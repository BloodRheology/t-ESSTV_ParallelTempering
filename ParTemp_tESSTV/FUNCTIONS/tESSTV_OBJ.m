function [obj,pred] = tESSTV_OBJ(parVec,exp)

% These are constants in our current model
at=1;
d= 1/2;
m= 1.5;
GammaINF=1;

par.mu0 = parVec(1);
par.muinf = parVec(2);
par.tauC = parVec(3);

par.tr1 = parVec(4);
par.tr2 = parVec(5);
par.muR = parVec(6);
par.sigy0 = parVec(7);

par.taulam = parVec(8);
par.GR = parVec(9);
par.GC = parVec(10);

% Organize experimental data
times = exp.times;
shears = exp.shears;
stresses = exp.stresses;


% Steps shear rates
gam_a = exp.gam_a; % 1/s initial shear rates
gam_b = exp.gam_b; % 1/s target shear rates

% Steady State analytical solution
[ERROR_SS,SSpred] = SteadyState_OBJ(parVec,exp);

pred.stressTOT = SSpred.SigmaTOT;
pred.stressElas = SSpred.SigmaR;
pred.stressVisc = SSpred.SigmaVISC;
pred.lambdaSS = SSpred.LambdaSS;

%% Transience
IC_TIME_VECT=[30;30;30;3;2;2];

% Caculating IC for transient
parfor ii=1:length(gam_a)
    SHEAR=gam_a(ii);
    ZZZZ = 0;
    Lambda00=1;
    Gamma0=GammaINF*Lambda00;

    StressVE0=0.001;
    Stress110=0;

    StressRyx0=0.001;
    StressRxx0=0;
    StressRyy0=0;
    StressRzz0=0;

    y0=[Lambda00 Gamma0 StressVE0 Stress110 StressRyx0 StressRxx0 StressRyy0 StressRzz0];
    tmax=IC_TIME_VECT(ii);
    t0=[0 tmax];

    consVec = [at,d,m,GammaINF,ZZZZ,SHEAR];
    dataVec = {times{ii},shears{ii}};
    
    tstart = tic;
    options = odeset('RelTol',1e-6,'AbsTol',1e-6,'Stats','off','Events',@(t,X)myEvent(t,X,tstart));
    fun = @(t,y) mHAWB_ODE(t,y,par,consVec,dataVec);
    [T, Y] = ode23s(fun,t0,y0,options);

    % Replace with steady state values
    Time = T;
    Lambda = Y(:,1);
    Gamma = Y(:,2);
    StressVE = Y(:,3);
    Stress111 = Y(:,4);
    StressR = Y(:,5);

    StressRxx = Y(:,6);
    StressRyy = Y(:,7);
    StressRzz = Y(:,8);

    LAMBDA_IC(ii,1)=(Lambda(end,1));
    GAMMA_IC(ii,1)=Gamma(end,1);
    STRESS_ICyx(ii,1)=(StressVE(end,1));
    STRESS_ICxx(ii,1)=(Stress111(end,1));
    RSTRESS_ICyx(ii,1)=StressR(end,1);
    RSTRESS_ICxx(ii,1)=StressRxx(end,1);
    RSTRESS_ICyy(ii,1)=StressRyy(end,1);
    RSTRESS_ICzz(ii,1)=StressRzz(end,1);
end



ERROR = zeros(length(gam_a),1);
ERROR_b = zeros(length(gam_a),1);
t = cell(length(gam_a),1);
sig_b = cell(length(gam_a),1);
lam_b = cell(length(gam_a),1);
sig_xx = cell(length(gam_a),1);
sig_yy = cell(length(gam_a),1);
sig_zz = cell(length(gam_a),1);
Strain_es = cell(length(gam_a),1);
Shear_ps = cell(length(gam_a),1);

parfor ii=1:length(gam_b)

    ZZZZ=1;

    SHEAR=gam_b(ii);
    Lambda00=LAMBDA_IC(ii,1);
    Gamma0=GAMMA_IC(ii,1);
    StressVE0=STRESS_ICyx(ii,1);
    Stress110=STRESS_ICxx(ii,1);
    StressRyx0=RSTRESS_ICyx(ii,1);
    StressRxx0=RSTRESS_ICxx(ii,1);
    StressRyy0=RSTRESS_ICyy(ii,1);
    StressRzz0=RSTRESS_ICzz(ii,1);

    ShearVECT=shears{ii};
    TimeVECT=times{ii};

    t0=times{ii};
    StressSOLN=stresses{ii};
    y0=[Lambda00 Gamma0 StressVE0 Stress110 StressRyx0 StressRxx0 StressRyy0 StressRzz0];

    consVec = [at,d,m,GammaINF,ZZZZ,SHEAR];
    dataVec = {TimeVECT,ShearVECT};
    
    tstart = tic;
    options = odeset('RelTol',1e-6,'AbsTol',1e-6,'Stats','off','Events',@(t,X)myEvent(t,X,tstart));
    fun = @(t,y) mHAWB_ODE(t,y,par,consVec,dataVec);
    [T, Y] = ode23s(fun,t0,y0,options);


    Time=T;
    Lambda=Y(:,1);
    Gamma=Y(:,2);
    StressVE=Y(:,3);
    Stress111=Y(:,4);
    StressR=Y(:,5);
    StressRxx=Y(:,6);
    StressRyy=Y(:,7);
    StressRzz=Y(:,8);

    StressTOT = StressR + StressVE;
    StressTOTXX = StressRxx + Stress111;
    Strain_e=Gamma;

   
    % Rogers et al addition:
    ShearREAL = [];
    for iiii=1:length(Lambda)
        if Time(iiii,1)<2
            ShearREAL(iiii,1)=interp1(t0, ShearVECT, Time(iiii,1));
        else
            ShearREAL(iiii,1)=SHEAR;
        end
        if ShearREAL(iiii,1)<0
            ShearREAL(iiii,1)=SHEAR;
        end

    end

    GMax0=min( (GammaINF.*(Lambda)), GammaINF);
    Shear_p=ShearREAL./(2-min((Gamma),GMax0)./GMax0);


    t{ii} = Time;
    sig_b{ii} = StressTOT;
    lam_b{ii} = Lambda;
    sig_xx{ii} = StressTOTXX;
    sig_yy{ii} = StressRyy;
    sig_zz{ii} = StressRzz;
    Strain_es{ii} = Strain_e;
    Shear_ps{ii} = Shear_p;

    for abc=1:length(Time)
        ERROR_b(ii) = ERROR_b(ii) + (sig_b{ii}(abc,1) - StressSOLN(abc,1))^2;
        if Time(abc,1)<2
            mult=5;%3.5;
            ERROR(ii) = ERROR(ii) + mult*(sig_b{ii}(abc,1) - StressSOLN(abc,1))^2;
        else
            mult=1;
            ERROR(ii) = ERROR(ii) + mult*(sig_b{ii}(abc,1)-StressSOLN(abc,1))^2;
        end
    end
    ERROR(ii) = (sqrt(ERROR(ii))/exp.NormVal(ii))/length(sig_b{ii});
    ERROR_b(ii) = sqrt(ERROR_b(ii))/length(sig_b{ii});


end

ERROR_TRANS = sum(ERROR)/6;
ERROR_TRANS_COMP = sum(ERROR_b)/6;

ERROR_TOT=ERROR_TRANS+ERROR_SS;

obj = ERROR_TOT;

pred.t = t;
pred.s = sig_b;
pred.lam = lam_b;
pred.sig_xx = sig_xx;
pred.sig_yy = sig_yy;
pred.sig_zz = sig_zz;
pred.Strain_es = Strain_es;
pred.Shear_ps = Shear_ps;

pred.ErrorSS = ERROR_SS;
pred.ErrorTRANS = ERROR_TRANS;

end