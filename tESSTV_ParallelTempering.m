%PAR_SIM_ANNEAL_SIAM_CAT_CRACK   26JUN18reversible reaction...using
%
clear; clc; format longG;  pause on;

global   Mu0_c MuINF_c TauC tr1 tr2 at d m MuR Sigy0 Gc GR TauLAM TauG G0 SHEAR GammaINF ALPHA tmax TimeVECT ShearVECT ZZZZ
Z=cputime;

% STEADY STATE
filename = 'SS.xlsx'; % Create an excel file with the first column
% containing the shear rate and the second column containing shear stress
% for the steady state experimental data

DATA = xlsread(filename);

ExpShearSS=DATA(:,1);  ExpStressSS=DATA(:,2);

SS_TIME_VECT=[40; 30; 20; 16; 10; 6; 4; 2; 2; 2; 2; 2; 2; 2; 2; 2];
clear DATA
%% STEP and TRIANGLE DATA:

filename = 'Steps.xlsx'; % Creat an excel file that contains 6 tabs in the
% following order: ["Down 5.0 -> 0.1","Down 10.0 -> 0.1","Down 20.0 -> 0.1",
% "Up 0.1 -> 5.0","Up 0.1 -> 10.0","Up 0.1 -> 20.0"]

names = ["Down_5p0_0p1","Down_10p0_0p1","Down_20p0_0p1", ...
    "Up_0p1_5p0","Up_0p1_10p0","Up_0p1_20p0"];

for i=1:length(names)
    DATA{i} = xlsread(filename,i);
end

TimeD1=DATA{1}(:,1);  StressD1=DATA{1}(:,2); tspan_SOLN3=TimeD1; ShearD1=DATA{1}(:,4);
TimeD2=DATA{2}(:,1);  StressD2=DATA{2}(:,2); tspan_SOLN4=TimeD2; ShearD2=DATA{2}(:,4);
TimeD3=DATA{3}(:,1);  StressD3=DATA{3}(:,2); tspan_SOLN8=TimeD3; ShearD3=DATA{3}(:,4);

TimeU1=DATA{4}(:,1);  StressU1=DATA{4}(:,2); tspan_SOLN1=TimeU1; ShearU1=DATA{4}(:,4);
TimeU2=DATA{5}(:,1);  StressU2=DATA{5}(:,2)-.005; tspan_SOLN2=TimeU2; ShearU2=DATA{5}(:,4);
TimeU3=DATA{6}(:,1);  StressU3=DATA{6}(:,2)-.01; tspan_SOLN7=TimeU3; ShearU3=DATA{6}(:,4);
%%
no=length(StressU1);

NormVal1=max(StressU1); NormVal2=max(StressU2); NormVal3=max(StressU3);
NormVal4=max(StressD1); NormVal5=max(StressD2); NormVal6=max(StressD3);

SHEAR1a=0.1;    SHEAR1b=5;
SHEAR2a=0.1;    SHEAR2b=10;
SHEAR7a=0.1;    SHEAR7b=20;

SHEAR3a=5;        SHEAR3b=.1;
SHEAR4a=10;       SHEAR4b=.1;
SHEAR8a=20;       SHEAR8b=.1;
%%

NB=10;
TotTime=1;   EPS_LIMIT=.01;  STOPPER=0;  DUMMY_VAR=0;  No_N_Ex_COLD=5;   CONC_CHAR=.5;  hCRIT=.0075;
EbMAX=1;  EbMIN=1/100000; E_RATIO=EbMIN/EbMAX;   NCHUNK_MIN=25;  NCHUNK_MAX=175;  NCHUNK_RATIO=NCHUNK_MAX/NCHUNK_MIN; hmin=.001;  hmax=.01;  hRATIO=hmin/hmax;   Eps_TEST=1.2e-8; h_TEST=0.0005; ERROR_MAX=0.001;
for i=1:NB   %INITIALIZING EVERYTING
    if i==1

        Temp(1,i)=EbMAX/EbMAX;
        NCHUNK(1,i)=NCHUNK_MIN;
        Epsilon(1,i)=(ERROR_MAX*(Temp(1,i)/Temp(1,1)))/CONC_CHAR;
        Error(1,i)=100000000;
        Guess(1,i)=2.01342182702437*rand;         % GR
        Guess(2,i)=5.69711803882554e-05*rand;     % Gc    
        Guess(3,i)=0.00892264342441323*rand;      % TauLAM
        Guess(4,i)=0.244077207493339*rand;     % tr1
        Guess(5,i)=0.132087656312376*rand;     % tr2
        Guess(6,i)=0.431143585306016*rand;        % Mu0_c
        Guess(7,i)=0.00475528117049142*rand;     % MuINF_c
        Guess(8,i)=57.0172848812029*rand;       % TauC
        Guess(9,i)=0.0187072577413047*rand;        % MuR
        Guess(10,i)=3.80017755774228e-06*rand;     % Sigy0


    else

        Temp(1,i)=Temp(1,i-1)*(E_RATIO)^(1/(NB-1)) ;%Temp(1, i-1)/4;
        NCHUNK(1,i)=floor(NCHUNK(1,i-1)*(NCHUNK_RATIO)^(1/(NB-1)));
        Epsilon(1,i)=(ERROR_MAX*(Temp(1,i)/Temp(1,1)))/CONC_CHAR;

        Error(1,i)=100000;
        Guess(1,i)=Guess(1,1);     %n
        Guess(2,i)=Guess(2,i);
        Guess(3,i)=Guess(3,i);
        Guess(4,i)=Guess(4,i);     
        Guess(5,i)=Guess(5,i);
        Guess(6,i)=Guess(6,i);
        Guess(7,i)=Guess(7,i);     
        Guess(8,i)=Guess(8,i);
        Guess(9,i)=Guess(9,i);
        Guess(9,i)=Guess(10,i);


    end
end

STOP_CRIT=Temp(1,NB);

TotTime=1;

at=1;
d= 1/2;
m= 1.5;
GammaINF=1;

%Temp(1,NB)=.0000000000000000000000000000000000001;
ANNEALLIMIT=10000;  ERROR_TOT_BEST=10000;   ZM=500;
%NCHUNK=[ ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM];

LIMIT=10;
ICOUNT=0;

GR_BEST = Guess(1,1);
Gc_BEST = Guess(2,1);
TauLAM_BEST = Guess(3,1);
tr1_BEST = Guess(4,1);
tr2_BEST = Guess(5,1);
Mu0_c_BEST = Guess(6,1);
MuINF_c_BEST = Guess(7,1);
TauC_BEST = Guess(8,1);
MuR_BEST = Guess(9,1);
Sigy0_BEST = Guess(10,1);

ERROR_SS_BEST = 10000;
ERROR_TRANS_BEST = 10000;


for i=1:LIMIT



    if i<=LIMIT
        ICOUNT=ICOUNT+1;



        if rem(ICOUNT, 1)==0
            i
            ERROR_SS_BEST
            ERROR_TRANS_BEST
            ERROR_TOT_BEST
            GR_BEST
            Gc_BEST
            TauLAM_BEST
        end

        for iii=1:NB-1           
            if rem(ICOUNT, NCHUNK(iii))==0   %switching the best values to the coldest Temp
                if iii==1
                    option=1;
                elseif iii>1 & iii<NB
                    option=2;
                end
                if option==1
                    PROB=exp( (1/Temp(1,2)-1/Temp(1,1))* (Error(1,2)-Error(1,1)) ); %Addition modification of sometimes exchanging bad parameters to the right...
                    X=rand;
                    %pause
                    %if k+1<NB
                    k=iii;
                    % This is the conditional statement to accept or reject
                    % but what are we accepting or rejecting?
                    if Error(1,k)<Error(1,k+1)  | rand < PROB
                        DUMMY=Error;  DUMMYG=Guess;
                        Error(1,k)=DUMMY(1,k+1);
                        Error(1,k+1)=DUMMY(1,k);

                        Guess(1,k)=DUMMYG(1,k+1);
                        Guess(2,k)=DUMMYG(2,k+1);
                        Guess(3,k)=DUMMYG(3,k+1);
                        Guess(4,k)=DUMMYG(4,k+1);
                        Guess(5,k)=DUMMYG(5,k+1);
                        Guess(6,k)=DUMMYG(6,k+1);
                        Guess(7,k)=DUMMYG(7,k+1);
                        Guess(8,k)=DUMMYG(8,k+1);
                        Guess(9,k)=DUMMYG(9,k+1);
                        Guess(10,k)=DUMMYG(10,k+1);


                        Guess(1,k+1)=DUMMYG(1,k);
                        Guess(2,k+1)=DUMMYG(2,k);
                        Guess(3,k+1)=DUMMYG(3,k);
                        Guess(4,k+1)=DUMMYG(4,k);
                        Guess(5,k+1)=DUMMYG(5,k);
                        Guess(6,k+1)=DUMMYG(6,k);
                        Guess(7,k+1)=DUMMYG(7,k);
                        Guess(8,k+1)=DUMMYG(8,k);
                        Guess(9,k+1)=DUMMYG(9,k);
                        Guess(10,k+1)=DUMMYG(10,k);

                    end
                elseif option==2
                    PROB=exp( (1/Temp(1,iii+1)-1/Temp(1,iii))* (Error(1,iii+1)-Error(1,iii)));
                    X=rand;
                    k=iii;
                    if Error(1,k)<Error(1,k+1)  | rand < PROB
                        DUMMY=Error;  DUMMYG=Guess;
                        Error(1,k)=DUMMY(1,k+1);
                        Error(1,k+1)=DUMMY(1,k);

                        Guess(1,k)=DUMMYG(1,k+1);
                        Guess(2,k)=DUMMYG(2,k+1);
                        Guess(3,k)=DUMMYG(3,k+1);
                        Guess(4,k)=DUMMYG(4,k+1);
                        Guess(5,k)=DUMMYG(5,k+1);
                        Guess(6,k)=DUMMYG(6,k+1);
                        Guess(7,k)=DUMMYG(7,k+1);
                        Guess(8,k)=DUMMYG(8,k+1);
                        Guess(9,k)=DUMMYG(9,k+1);
                        Guess(10,k)=DUMMYG(10,k+1);


                        Guess(1,k+1)=DUMMYG(1,k);
                        Guess(2,k+1)=DUMMYG(2,k);
                        Guess(3,k+1)=DUMMYG(3,k);
                        Guess(4,k+1)=DUMMYG(4,k);
                        Guess(5,k+1)=DUMMYG(5,k);
                        Guess(6,k+1)=DUMMYG(6,k);
                        Guess(7,k+1)=DUMMYG(7,k);
                        Guess(8,k+1)=DUMMYG(8,k);
                        Guess(9,k+1)=DUMMYG(9,k);
                        Guess(10,k+1)=DUMMYG(10,k);

                    end
                end
            end
        end


        for j=1:NB
            GUESS1=Guess(1,j);
            GUESS2=Guess(2,j);
            GUESS3=Guess(3,j);
            GUESS4=Guess(4,j);
            GUESS5=Guess(5,j);
            GUESS6=Guess(6,j);
            GUESS7=Guess(7,j);
            GUESS8=Guess(8,j);
            GUESS9=Guess(9,j);
            GUESS10=Guess(10,j);


       

            if ICOUNT<=1/2*LIMIT
                MULT=.15;
            elseif ICOUNT>1/2*LIMIT & ICOUNT <=3/4*LIMIT
                MULT=.1;
            else
                MULT=.05;
            end



            if ICOUNT==1  %no change to parameter values during first time through
                GR=GUESS1;
                Gc=GUESS2;
                TauLAM=GUESS3;
                tr1=GUESS4;
                tr2=GUESS5;
                Mu0_c=GUESS6;
                MuINF_c=GUESS7;
                TauC=GUESS8;
                MuR=GUESS9;
                Sigy0=GUESS10;
            else

%                 GR=(sqrt(GUESS1)+MULT*(.5-rand))^2;   %first three are transient par
%                 Gc=(sqrt(GUESS2)+MULT*(.5-rand))^2;
%                 TauLAM=(sqrt(GUESS3)+MULT*(.5-rand))^2;
% 
%                 tr1=(sqrt(GUESS4)+MULT*.025*(.5-rand))^2;   % these are 'steady state' par.
%                 tr2=(sqrt(GUESS5)+MULT*.015*(.5-rand))^2;
%                 Mu0_c=(sqrt(GUESS6)+MULT*.025*(.5-rand))^2;
%                 MuINF_c=(sqrt(GUESS7)+MULT*.01*(.5-rand))^2;
%                 TauC=(sqrt(GUESS8)+MULT*.01*(.5-rand))^2;
%                 MuR=(sqrt(GUESS9)+MULT*.025*(.5-rand))^2;
%                 Sigy0=(sqrt(GUESS10)+MULT*.025*(.5-rand))^2;

                GR=log(exp(MULT*randn)*exp(GUESS1)); % GR
                Gc=log(exp(MULT*randn)*exp(GUESS2)); % Gc
                TauLAM=log(exp(MULT*randn)*exp(GUESS3)); % TauLAM
                
                tr1=log(exp(MULT*0.025*randn)*exp(GUESS4)); % tr1
                tr2=log(exp(MULT*0.015*randn)*exp(GUESS5)); % tr2
                Mu0_c=log(exp(MULT*0.025*randn)*exp(GUESS6)); % Mu0_c
                MuINF_c=log(exp(MULT*0.025*randn)*exp(GUESS7)); % MuINF_c
                TauC=log(exp(MULT*0.01*randn)*exp(GUESS8)); % TauC
                MuR=log(exp(MULT*0.025*randn)*exp(GUESS9)); % MuR
                Sigy0=log(exp(MULT*0.025*randn)*exp(GUESS10)); % Sigy0
            end


            %Steady State IC:
            SHEAR=300;
            ZZZZ=0;
            Lambda00=1;%(tr2*SHEAR_START^d+  1)/(tr1*SHEAR_START^at +  tr2*SHEAR_START^d+  1);
            Gamma0=GammaINF*Lambda00;
            StressVE0=0.001;%((Mu0_c-MuINF_c)/(1+TauC*SHEAR_START)+MuINF_c)*SHEAR_START;
            Stress110=0;%2*((Mu0_c-MuINF_c)/(1+TauC*SHEAR5)+MuINF_c)*SHEAR_START^2/Gc;
            StressRyx0=0.001;%(Lambda00*Sigy0 + MuR*SHEAR1a*Lambda00^m);
            StressRxx0=0;  StressRyy0=0;  StressRzz0=0; ERROR_SS=0; ERROR_SS_COMP=0;

            y0=[Lambda00 Gamma0 StressVE0 Stress110 StressRyx0 StressRxx0 StressRyy0 StressRzz0];
            t0=[0 3];
            Time1=0; StressTOT1=0; Lambda1=0; Gamma1=0; GMax1=0;  T=0; Y=0; Shear1=0; ShearP1=0; StressR1=0; StressER1=0; StressVER1=0; StressVE1=0;
            tstart = tic;
            options = odeset('RelTol',1e-6,'AbsTol',[1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6],'Stats','off','Events',@(t,X)myEvent(t,X,tstart));

            
            [T, Y] = ode23s(@mHAWB_ODE,t0,y0,options);


            Z=length(T);
            Time1=T; Lambda00=Y(Z,1); Gamma0=Y(Z,2);  StressVE0=Y(Z,3); Stress110=Y(Z,4); StressRyx0=Y(Z,5);
            StressRxx0=Y(Z,6);  StressRyy0=Y(Z,7);  StressRzz0=Y(Z,8);

            for ii=1:length(ExpShearSS)
                SHEAR=ExpShearSS(ii,1);
                y0=[Lambda00 Gamma0 StressVE0 Stress110 StressRyx0 StressRxx0 StressRyy0 StressRzz0];
                tmax=SS_TIME_VECT(ii,1);
                t0=[0 tmax];
                Time1=0; StressTOT1=0; Lambda1=0; Gamma1=0; GMax1=0;  T=0; Y=0; Shear1=0; ShearP1=0; StressR1=0; StressER1=0; StressVER1=0; StressVE1=0;
                StressRxx=0; StressRyy=0; StressRzz=0;
                tstart = tic;
                options = odeset('RelTol',1e-6,'AbsTol',[1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6],'Stats','off','Events',@(t,X)myEvent(t,X,tstart));
                [T, Y] = ode23s(@mHAWB_ODE,t0,y0,options);
                Time1=T; Lambda1=Y(:,1); Gamma1=Y(:,2);  StressVE1=Y(:,3); Stress111=Y(:,4); StressR1=Y(:,5);
                StressRxx=Y(:,6);  StressRyy=Y(:,7);  StressRzz=Y(:,8);
                StressTOT1=StressR1+ StressVE1;
                StressXX1=Stress111+StressRxx;
                StressYY1=StressRyy;
                StressN1=StressXX1-StressYY1;

                Z=length(StressTOT1);
                STRESSb(ii,1)=(StressTOT1(Z,1));
                ESTRESSb(ii,1)=(StressR1(Z,1));
                VSTRESSb(ii,1)=(StressVE1(Z,1));
                LAMBDAb(ii,1)=(Lambda1(Z,1));
                STRESSN1b(ii,1)=(StressN1(Z,1));
                % ERROR_SS=ERROR_SS + ((STRESSb(i,1)-ExpStressSS(i,1))/1)^2;
                ERROR_SS=ERROR_SS+((STRESSb(ii,1)-ExpStressSS(ii,1))/ExpStressSS(ii,1))^2;
                % ERROR_SS=ERROR_SS + ((STRESSb(i,1)-ExpStressSS(i,1))/ExpStressSS(i,1))^2;
            end
            ERROR_SS=(sqrt(ERROR_SS)/length(ExpStressSS));
            % ERROR_SS_COMP=(sqrt(ERROR_SS_COMP)/length(ExpStressSS));

            V_VECT(1,1)=SHEAR1a;
            V_VECT(2,1)=SHEAR3a;  V_VECT(3,1)=SHEAR4a;  V_VECT(4,1)=SHEAR8a;
            IC_TIME_VECT=[30;3;2;2];
            % Caculating IC for transient
            for ii=1:length(V_VECT)
                SHEAR=V_VECT(ii,1);
                y0=[Lambda00 Gamma0 StressVE0 Stress110 StressRyx0 StressRxx0 StressRyy0 StressRzz0];
                tmax=IC_TIME_VECT(ii,1);
                t0=[0 tmax];
                Time1=0; StressTOT1=0; Lambda1=0; Gamma1=0; GMax1=0;  T=0; Y=0; Shear1=0; ShearP1=0; StressR1=0; StressER1=0; StressVER1=0; StressVE1=0;
                tstart = tic;
                options = odeset('RelTol',1e-6,'AbsTol',[1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6],'Stats','off','Events',@(t,X)myEvent(t,X,tstart));
                [T, Y] = ode23s(@mHAWB_ODE,t0,y0,options);
                Time1=T; Lambda1=Y(:,1); Gamma1=Y(:,2);  StressVE1=Y(:,3); Stress111=Y(:,4); StressR1=Y(:,5);
                StressRxx=Y(:,6);  StressRyy=Y(:,7);  StressRzz=Y(:,8);
                StressTOT1=StressR1+ StressVE1;

                Z=length(StressTOT1);
                LAMBDA_IC(ii,1)=(Lambda1(Z,1));
                GAMMA_IC(ii,1)=Gamma1(Z,1);
                STRESS_ICyx(ii,1)=(StressVE1(Z,1));
                STRESS_ICxx(ii,1)=(Stress111(Z,1));
                RSTRESS_ICyx(ii,1)=StressR1(Z,1);
                RSTRESS_ICxx(ii,1)=StressRxx(Z,1);
                RSTRESS_ICyy(ii,1)=StressRyy(Z,1);
                RSTRESS_ICzz(ii,1)=StressRzz(Z,1);
            end

            %Transient
            for ii=1:6
                ZZZZ=1;
                if ii==1
                    SHEAR=SHEAR1b;
                    Lambda00=LAMBDA_IC(1,1);
                    Gamma0=GAMMA_IC(1,1);
                    StressVE0=STRESS_ICyx(1,1);
                    Stress110=STRESS_ICxx(1,1);
                    StressRyx0=RSTRESS_ICyx(1,1);
                    StressRxx0=RSTRESS_ICxx(1,1);
                    StressRyy0=RSTRESS_ICyy(1,1);
                    StressRzz0=RSTRESS_ICzz(1,1);
                    ShearVECT=ShearU1; TimeVECT=tspan_SOLN1;

                    t0=TimeU1; StressSOLN1=StressU1;
                    y0=[Lambda00 Gamma0 StressVE0 Stress110 StressRyx0 StressRxx0 StressRyy0 StressRzz0];

                elseif ii==2
                    SHEAR=SHEAR2b;
                    Lambda00=LAMBDA_IC(1,1);
                    Gamma0=GAMMA_IC(1,1);
                    StressVE0=STRESS_ICyx(1,1);
                    Stress110=STRESS_ICxx(1,1);
                    StressRyx0=RSTRESS_ICyx(1,1);
                    StressRxx0=RSTRESS_ICxx(1,1);
                    StressRyy0=RSTRESS_ICyy(1,1);
                    StressRzz0=RSTRESS_ICzz(1,1);
                    ShearVECT=ShearU2; TimeVECT=tspan_SOLN2;
                    t0=TimeU2; StressSOLN2=StressU2;
                    y0=[Lambda00 Gamma0 StressVE0 Stress110 StressRyx0 StressRxx0 StressRyy0 StressRzz0];

                elseif ii==3
                    SHEAR=SHEAR7b;
                    Lambda00=LAMBDA_IC(1,1);
                    Gamma0=GAMMA_IC(1,1);
                    StressVE0=STRESS_ICyx(1,1);
                    Stress110=STRESS_ICxx(1,1);
                    StressRyx0=RSTRESS_ICyx(1,1);
                    StressRxx0=RSTRESS_ICxx(1,1);
                    StressRyy0=RSTRESS_ICyy(1,1);
                    StressRzz0=RSTRESS_ICzz(1,1);
                    ShearVECT=ShearU3; TimeVECT=tspan_SOLN7;
                    t0=TimeU3; StressSOLN3=StressU3;
                    y0=[Lambda00 Gamma0 StressVE0 Stress110 StressRyx0 StressRxx0 StressRyy0 StressRzz0];

                elseif ii==4
                    SHEAR=SHEAR3b;
                    Lambda00=LAMBDA_IC(2,1);
                    Gamma0=GAMMA_IC(2,1);
                    StressVE0=STRESS_ICyx(2,1);
                    Stress110=STRESS_ICxx(2,1);
                    StressRyx0=RSTRESS_ICyx(2,1);
                    StressRxx0=RSTRESS_ICxx(2,1);
                    StressRyy0=RSTRESS_ICyy(2,1);
                    StressRzz0=RSTRESS_ICzz(2,1);
                    ShearVECT=ShearD1; TimeVECT=tspan_SOLN3;
                    t0=TimeD1; StressSOLN4=StressD1;
                    y0=[Lambda00 Gamma0 StressVE0 Stress110 StressRyx0 StressRxx0 StressRyy0 StressRzz0];

                elseif ii==5
                    SHEAR=SHEAR4b;   %changed from 3 to 2
                    Lambda00=LAMBDA_IC(3,1);
                    Gamma0=GAMMA_IC(3,1);
                    StressVE0=STRESS_ICyx(3,1);
                    Stress110=STRESS_ICxx(3,1);
                    StressRyx0=RSTRESS_ICyx(3,1);
                    StressRxx0=RSTRESS_ICxx(3,1);
                    StressRyy0=RSTRESS_ICyy(3,1);
                    StressRzz0=RSTRESS_ICzz(3,1);
                    ShearVECT=ShearD2; TimeVECT=tspan_SOLN4;
                    t0=TimeD2; StressSOLN5=StressD2;
                    y0=[Lambda00 Gamma0 StressVE0 Stress110 StressRyx0 StressRxx0 StressRyy0 StressRzz0];

                elseif ii==6
                    SHEAR=SHEAR8b;
                    Lambda00=LAMBDA_IC(4,1);
                    Gamma0=GAMMA_IC(4,1);
                    StressVE0=STRESS_ICyx(4,1);
                    Stress110=STRESS_ICxx(4,1);
                    StressRyx0=RSTRESS_ICyx(4,1);
                    StressRxx0=RSTRESS_ICxx(4,1);
                    StressRyy0=RSTRESS_ICyy(4,1);
                    StressRzz0=RSTRESS_ICzz(4,1);
                    ShearVECT=ShearD3; TimeVECT=tspan_SOLN8;
                    t0=TimeD3; StressSOLN6=StressD3;
                    y0=[Lambda00 Gamma0 StressVE0 Stress110 StressRyx0 StressRxx0 StressRyy0 StressRzz0];
                end


                Time1=0; StressTOT1=0; Lambda1=0; Gamma1=0; GMax1=0;  T=0; Y=0; Shear1=0; ShearP1=0; StressR1=0; StressER1=0; StressVER1=0; StressVE1=0;
                tstart = tic;
                options = odeset('RelTol',1e-6,'AbsTol',[1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6],'Stats','off','Events',@(t,X)myEvent(t,X,tstart));
                [T, Y] = ode23s(@mHAWB_ODE,t0,y0,options);
                Time1=T; Lambda1=Y(:,1); Gamma1=Y(:,2);  StressVE1=Y(:,3); Stress111=Y(:,4); StressR1=Y(:,5);
                StressRxx=Y(:,6);  StressRyy=Y(:,7);  StressRzz=Y(:,8);
                StressTOT1=StressR1+ StressVE1;
                StressTOTXX=StressRxx+Stress111;
                Strain_e=Gamma1;

                %Rogers et al addition:
                ShearREAL = [];
                for iiii=1:length(Lambda1)
                    if Time1(iiii,1)<2
                        ShearREAL(iiii,1)=interp1(TimeVECT, ShearVECT, Time1(iiii,1));
                    else
                        ShearREAL(iiii,1)=SHEAR;
                    end
                    if ShearREAL(iiii,1)<0
                        ShearREAL(iiii,1)=SHEAR;
                    end

                end
                GMax0=min( (GammaINF.*(Lambda1)), GammaINF);
                Shear_p=ShearREAL./(2-min((Gamma1),GMax0)./GMax0);

                if ii==1
                    t1=Time1;
                    s1b=StressTOT1;
                    lam1b=Lambda1;
                    s1xx=StressTOTXX;
                    s1yy=StressRyy;
                    s1zz=StressRzz;
                    Strain_e1= Strain_e;
                    Shear_p1=Shear_p;
                    ERROR1=0; ERROR1b=0;
                    for abc=1:length(Time1)
                        ERROR1b=ERROR1b +(s1b(abc,1)-StressSOLN1(abc,1))^2;
                        if Time1(abc,1)<2
                            mult=5;%3.5;
                            ERROR1=ERROR1+mult*(s1b(abc,1)-StressSOLN1(abc,1))^2;
                        else
                            mult=1;
                            ERROR1=ERROR1+mult*(s1b(abc,1)-StressSOLN1(abc,1))^2;
                        end
                    end
                    ERROR1=(sqrt(ERROR1)/NormVal1)/length(s1b);
                    ERROR1b=sqrt(ERROR1b)/length(s1b);
                    %ERROR1=norm(s1b-StressSOLN1)/1;%length(s1b);

                elseif ii==2
                    t2=Time1;
                    s2b=StressTOT1;
                    lam2b=Lambda1;
                    s2xx=StressTOTXX;
                    s2yy=StressRyy;
                    s2zz=StressRzz;
                    Strain_e2= Strain_e;
                    Shear_p2=Shear_p;
                    ERROR2=0;ERROR2b=0;
                    for abc=1:length(Time1)
                        ERROR2b=ERROR2b +(s2b(abc,1)-StressSOLN2(abc,1))^2;
                        if Time1(abc,1)<2
                            mult=5;%3.5;
                            ERROR2=ERROR2+mult*(s2b(abc,1)-StressSOLN2(abc,1))^2;
                        else
                            mult=1;
                            ERROR2=ERROR2+mult*(s2b(abc,1)-StressSOLN2(abc,1))^2;
                        end
                    end
                    ERROR2=(sqrt(ERROR2)/NormVal2)/length(s2b);
                    ERROR2b=sqrt(ERROR2b)/length(s2b);
                    %ERROR2=norm(s2b-StressSOLN2)/1;%length(s2b);

                elseif ii==3
                    %Step up from 0.1 to 20
                    t3=Time1;
                    s3b=StressTOT1;
                    Plasma3=StressVE1;
                    Rouleaux3=StressR1;
                    lam3b=Lambda1;
                    s3xx=StressTOTXX;
                    s3yy=StressRyy;
                    s3zz=StressRzz;
                    Strain_e3= Strain_e;
                    Shear_p3=Shear_p;
                    N1_3=StressTOTXX-StressRyy;
                    ERROR3=0; ERROR3b=0;
                    for abc=1:length(Time1)
                        ERROR3b=ERROR3b +(s3b(abc,1)-StressSOLN3(abc,1))^2;
                        if Time1(abc,1)<2
                            mult=10;%3.5;
                            ERROR3=ERROR3+mult*(s3b(abc,1)-StressSOLN3(abc,1))^2;
                        else
                            mult=1;
                            ERROR3=ERROR3+mult*(s3b(abc,1)-StressSOLN3(abc,1))^2;
                        end
                    end
                    ERROR3=(sqrt(ERROR3)/NormVal3)/length(s3b);
                    ERROR3b=sqrt(ERROR3b)/length(s3b);
                    %ERROR3=norm(s3b-StressSOLN3)/1;%length(s3b);

                elseif ii==4
                    t4=Time1;
                    s4b=StressTOT1;
                    lam4b=Lambda1;
                    s4xx=StressTOTXX;
                    s4yy=StressRyy;
                    s4zz=StressRzz;
                    Strain_e4= Strain_e;
                    Shear_p4=Shear_p;
                    ERROR4=0;ERROR4b=0;
                    for abc=1:length(Time1)
                        ERROR4b=ERROR4b +(s4b(abc,1)-StressSOLN4(abc,1))^2;
                        if Time1(abc,1)<2
                            mult=5;%5.5;
                            ERROR4=ERROR4+mult*(s4b(abc,1)-StressSOLN4(abc,1))^2;
                        else
                            mult=1;
                            ERROR4=ERROR4+mult*(s4b(abc,1)-StressSOLN4(abc,1))^2;
                        end
                    end
                    ERROR4=(sqrt(ERROR4)/NormVal4)/length(s4b);
                    ERROR4b=sqrt(ERROR4b)/length(s4b);
                    %ERROR4=norm(s4b-StressSOLN4)/1;%length(s4b);

                elseif ii==5
                    t5=Time1;
                    s5b=StressTOT1;
                    lam5b=Lambda1;
                    s5xx=StressTOTXX;
                    s5yy=StressRyy;
                    s5zz=StressRzz;
                    Strain_e5= Strain_e;
                    Shear_p5=Shear_p;
                    ERROR5=0;ERROR5b=0;
                    for abc=1:length(Time1)
                        ERROR5b=ERROR5b +(s5b(abc,1)-StressSOLN5(abc,1))^2;
                        if Time1(abc,1)<2
                            mult=5;%5.5;
                            ERROR5=ERROR5+mult*(s5b(abc,1)-StressSOLN5(abc,1))^2;
                        else
                            mult=1;
                            ERROR5=ERROR5+mult*(s5b(abc,1)-StressSOLN5(abc,1))^2;
                        end
                    end
                    ERROR5=(sqrt(ERROR5)/NormVal5)/length(s5b);
                    ERROR5b=sqrt(ERROR5b)/length(s5b);
                    %ERROR5=norm(s5b-StressSOLN5)/1;%length(s5b);

                elseif ii==6
                    %step down from 20 to 0.1
                    t6=Time1;
                    s6b=StressTOT1;
                    Plasma6=StressVE1;
                    Rouleaux6=StressR1;
                    lam6b=Lambda1;
                    s6xx=StressTOTXX;
                    s6yy=StressRyy;
                    s6zz=StressRzz;
                    Strain_e6= Strain_e;
                    Shear_p6=Shear_p;
                    N1_6=StressTOTXX-StressRyy;
                    ERROR6=0; ERROR6b=0;
                    for abc=1:length(Time1)
                        ERROR6b=ERROR6b +(s6b(abc,1)-StressSOLN6(abc,1))^2;
                        if Time1(abc,1)<2
                            mult=5;%5.5;
                            ERROR6=ERROR6+mult*(s6b(abc,1)-StressSOLN6(abc,1))^2;
                        else
                            mult=1;
                            ERROR6=ERROR6+mult*(s6b(abc,1)-StressSOLN6(abc,1))^2;
                        end
                    end
                    ERROR6=(sqrt(ERROR6)/NormVal6)/length(s6b);
                    ERROR6b=sqrt(ERROR6b)/length(s6b);
                    %ERROR6=norm(s6b-StressSOLN6)/1;%length(s6b);
                end
            end
            ERROR_TRANS=(ERROR1+ERROR2+ERROR3+ERROR4+ERROR5+ERROR6)/6;
            ERROR_TRANS_COMP=(ERROR1b+ERROR2b+ERROR3b+ERROR4b+ERROR5b+ERROR6b)/6;

            ERROR_TOT=ERROR_TRANS+ERROR_SS;
            %ERROR_TOT=.9*ERROR_TRANS_COMP + .1*ERROR_SS;

%%
            if i>1 & i<=ANNEALLIMIT
                ANNEALP=exp( (Error(1,j)-ERROR_TOT)/Temp(1,j));
            end

%             if ERROR_TOT<ERROR_TOT_BEST
% 
%                 GR0=GR;         GRBEST=GR;
%                 GC0=Gc;         GcBEST=Gc;
%                 TauLAM0=TauLAM; TauLAMBEST=TauLAM;
%                 tr10=tr1;         tr1BEST=tr1;
%                 tr20=tr2;         tr2BEST=tr2;
%                 Mu0_c0=Mu0_c;     Mu0_cBEST=Mu0_c;
%                 MuINF_c0=MuINF_c;     MuINF_cBEST=MuINF_c;
%                 TauC0= TauC;      TauCBEST=TauC;
%                 MuR0=MuR;         MuRBEST=MuR;
%                 Sigy00=Sigy0;     Sigy0BEST=Sigy0;
% 
%                 s1=s1b; lam1=lam1b;  N11=s1xx-s1yy;  N21=s1yy-s1zz;
%                 s2=s2b; lam2=lam2b;  N12=s2xx-s2yy;  N22=s2yy-s2zz;
%                 s3=s3b; lam3=lam3b;  N13=s3xx-s3yy;  N23=s3yy-s3zz;
%                 s4=s4b; lam4=lam4b;  N14=s4xx-s4yy;  N24=s4yy-s4zz;
%                 s5=s5b; lam5=lam5b;  N15=s5xx-s5yy;  N25=s5yy-s5zz;
%                 s6=s6b; lam6=lam6b;  N16=s6xx-s6yy;  N26=s6yy-s6zz;
% 
%                 STRESS=STRESSb;
%                 ESTRESS=ESTRESSb;
%                 VSTRESS=VSTRESSb;
%                 LAMBDA=LAMBDAb;
% 
%                 ERROR_SS_BEST=ERROR_SS;
%                 ERROR_TRANS_BEST=ERROR_TRANS;
%                 ERROR_TOT_BEST=ERROR_TOT;
% 
%                 ERROR_SS_COMP_BEST=ERROR_SS_COMP;
%                 ERROR_TRANS_COMP_BEST=ERROR_TRANS_COMP;
%                 
%                 iBEST=i;
% 
%             end

           
            if i>1 & i<=ANNEALLIMIT
                X=rand;
                if ERROR_TOT_BEST<Error(1,j)
                    if ERROR_TOT<ERROR_TOT_BEST
                        GR_BEST = GR;
                        Gc_BEST = Gc;
                        TauLAM_BEST = TauLAM;
                        tr1_BEST = tr1;
                        tr2_BEST = tr2;
                        Mu0_c_BEST = Mu0_c;
                        MuINF_c_BEST = MuINF_c;
                        TauC_BEST = TauC;
                        MuR_BEST = MuR;
                        Sigy0_BEST = Sigy0;    

                        s1=s1b; lam1=lam1b;  N11=s1xx-s1yy;  N21=s1yy-s1zz;
                        s2=s2b; lam2=lam2b;  N12=s2xx-s2yy;  N22=s2yy-s2zz;
                        s3=s3b; lam3=lam3b;  N13=s3xx-s3yy;  N23=s3yy-s3zz;
                        s4=s4b; lam4=lam4b;  N14=s4xx-s4yy;  N24=s4yy-s4zz;
                        s5=s5b; lam5=lam5b;  N15=s5xx-s5yy;  N25=s5yy-s5zz;
                        s6=s6b; lam6=lam6b;  N16=s6xx-s6yy;  N26=s6yy-s6zz;

                        STRESS=STRESSb;
                        ESTRESS=ESTRESSb;
                        VSTRESS=VSTRESSb;
                        LAMBDA=LAMBDAb;

                        ERROR_SS_BEST=ERROR_SS;
                        ERROR_TRANS_BEST=ERROR_TRANS;
                        ERROR_TOT_BEST=ERROR_TOT;

                        ERROR_SS_COMP_BEST=ERROR_SS_COMP;
                        ERROR_TRANS_COMP_BEST=ERROR_TRANS_COMP;
                        iBEST=i;
                    end
                    Guess(1,j)=GR;
                    Guess(2,j)=Gc;
                    Guess(3,j)=TauLAM;
                    Guess(4,j)=tr1;
                    Guess(5,j)=tr2;
                    Guess(6,j)=Mu0_c;
                    Guess(7,j)=MuINF_c;
                    Guess(8,j)=TauC;
                    Guess(9,j)=MuR;
                    Guess(10,j)=Sigy0;


                    Error(1,j)=ERROR_TOT;
                elseif ERROR_TOT>=Error(1,j) & X< ANNEALP
                    Guess(1,j)=GR;
                    Guess(2,j)=Gc;
                    Guess(3,j)=TauLAM;
                    Guess(4,j)=tr1;
                    Guess(5,j)=tr2;
                    Guess(6,j)=Mu0_c;
                    Guess(7,j)=MuINF_c;
                    Guess(8,j)=TauC;
                    Guess(9,j)=MuR;
                    Guess(10,j)=Sigy0;


                    Error(1,j)=ERROR_TOT;
                elseif ERROR_TOT>=Error(1,j) & X>=ANNEALP
                    %Error(1,j)=Error(1,j);
                end
            end
        end

        RESULTS(i,1)=GR_BEST;
        RESULTS(i,2)=Gc_BEST;
        RESULTS(i,3)=TauLAM_BEST;

        RESULTS(i,4)=tr1_BEST;
        RESULTS(i,5)=tr2_BEST;
        RESULTS(i,6)=Mu0_c_BEST;
        RESULTS(i,7)=MuINF_c_BEST;
        RESULTS(i,8)=TauC_BEST;
        RESULTS(i,9)=MuR_BEST;
        RESULTS(i,10)=Sigy0_BEST;

        RESULTS(i,11)=ERROR_TOT_BEST;
        RESULTS(i,12)=i;



    end

end

RunTime=cputime-Z

MetaData = ["Zero Shear Viscosity";
    "Infinite Shear Viscosity";
    "RBC Deformation Time Constant";
    "Yield Stress";
    "tr1";
    "tr2";
    "Rouleaux Viscosity";
    "Overall Structure Rebuild Time Constant";
    "RBC Elastic Modulus";
    "Rouleaux Elastic Modulus"];
Best_Fits = [Mu0_c_BEST;
            MuINF_c_BEST;
            TauC_BEST;
            Sigy0_BEST;
            tr1_BEST;
            tr2_BEST;
            MuR_BEST;
            TauLAM_BEST;
            Gc_BEST;
            GR_BEST];

% ADD IN UNITS
T = table(MetaData,Best_Fits)
filename = 'Metadata.xlsx';
writetable(T,filename,'Sheet','Rheology Parameters')


%% Plotting
figure(1);
loglog(ExpShearSS, ExpStressSS, 'ro',ExpShearSS, STRESS,'k-.',ExpShearSS, VSTRESS,'b-.',ExpShearSS, ESTRESS,'g-.','MarkerSize',8,'LineWidth',2.);
set(gca,'FontSize',14,'FontWeight','bold','linewidth',2, 'FontName','Times');
legend('DATA','tot stress','visc stress','rouleaux stress','Location','NorthWest');
xlabel('Shear Rate (1/s)');
ylabel('Stress (Pa)');
xlim([min(ExpShearSS) max(ExpShearSS)]);
ylim([.0001 5]);
