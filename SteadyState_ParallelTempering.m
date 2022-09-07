%ARMSTRONG_ethixo_HAWB_SS_26JAN21
%FOR: CDTs O'Malley and Bailey
clear; clc; format longG; pause on;
%Note that this mfile is to be used for fitting steady state HAWB
%parameters to data taken at UDEL on ARESG2 and data taken at USMA on
%DHR-3.

Z=cputime;

%SS DATA:
%This is where you read in the SS data to be fit
%Shear rate in first col./ stress in the second col.

filename = 'SS.xlsx'; % Create an excel file with the first column
% containing the shear rate and the second column containing shear stress
% for the steady state experimental data

DATA = xlsread(filename);

ExpShearSS=DATA(:,1);  ExpStressSS=DATA(:,2);  ExpViscSS=ExpStressSS./ExpShearSS; N=length(ExpShearSS);

%Initial guesses from paper:
Mu0_c=  0.145468906131364;
MuINF_c=  0.00316706382891868;
TauC= .351797011418833;
Sigy0=0.20498850301406;
tr1=  0.70951493191527;
tr2=  0.111778944141652;
a=1.0644714293727;
MuR= 0.0216741984005847;
m=1.5;
d=1/2;

%% INITIALIZING BOLTZMANN ENERGY LEVELS

NB=10;
EPS_LIMIT=.01;  STOPPER=0;  DUMMY_VAR=0;  No_N_Ex_COLD=5;   CONC_CHAR=.5;  hCRIT=.0075;
EbMAX=1;  EbMIN=1/100000; E_RATIO=EbMIN/EbMAX;   NCHUNK_MIN=25;  NCHUNK_MAX=175;
NCHUNK_RATIO=NCHUNK_MAX/NCHUNK_MIN; hmin=.001;  hmax=.01;  hRATIO=hmin/hmax;   Eps_TEST=1.2e-8;
h_TEST=0.0005; ERROR_MAX=0.001;
for i=1:NB   %INITIALIZING EVERYTHING
    if i==1

        Temp(1,i)=EbMAX/EbMAX;
        NCHUNK(1,i)=NCHUNK_MIN;
        Epsilon(1,i)=(ERROR_MAX*(Temp(1,i)/Temp(1,1)))/CONC_CHAR;
        Error(1,i)=100000;
        Guess(1,i)=Mu0_c; % Mu0_c     %n
        Guess(2,i)=MuINF_c; % MuINF_c
        Guess(3,i)=TauC;  % TauC
        Guess(4,i)=Sigy0;  % Sigy0
        Guess(5,i)=tr1;  % tr1
        Guess(6,i)=tr2;  % tr2
        Guess(7,i)=MuR;  % MuR

    else

        Temp(1,i)=Temp(1,i-1)*(E_RATIO)^(1/(NB-1)) ;%Temp(1, i-1)/4;
        NCHUNK(1,i)=floor(NCHUNK(1,i-1)*(NCHUNK_RATIO)^(1/(NB-1)));
        Epsilon(1,i)=(ERROR_MAX*(Temp(1,i)/Temp(1,1)))/CONC_CHAR;

        Error(1,i)=10000;
        Guess(1,i)=Guess(1,1);  % Mu0_c   %n
        Guess(2,i)=Guess(2,1);  % MuINF_c
        Guess(3,i)=Guess(3,1);  % TauC
        Guess(4,i)=Guess(4,1);  % Sigy0
        Guess(5,i)=Guess(5,1);  % tr1
        Guess(6,i)=Guess(6,1);  % tr2
        Guess(7,i)=Guess(7,1);  % MuR


    end
end

STOP_CRIT=Temp(1,NB);



%Temp(1,NB)=.0000000000000000000000000000000000001;
ANNEALLIMIT= 50000;  ERRORBEST=10000000000;   ZM=500;
%NCHUNK=[ ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM ZM];

LIMIT=100000;
ICOUNT=0;

for i=1:LIMIT



    if i<=LIMIT
        ICOUNT=ICOUNT+1;



        if rem(ICOUNT, 500)==0
            i
            ERRORBEST
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
                    if Error(1,k)<Error(1,k+1)  | rand < PROB
                        DUMMY=Error;  DUMMYG=Guess;
                        Error(1,k)=DUMMY(1,k+1);
                        Error(1,k+1)=DUMMY(1,k);

                        for nn = 1:length(Guess(:,1))
                            Guess(nn,k)=DUMMYG(nn,k+1);
                            Guess(nn,k+1)=DUMMYG(nn,k);
                        end

                    end
                elseif option==2
                    PROB=exp( (1/Temp(1,iii+1)-1/Temp(1,iii))* (Error(1,iii+1)-Error(1,iii)));
                    X=rand;
                    k=iii;
                    if Error(1,k)<Error(1,k+1)  | rand < PROB
                        DUMMY=Error;  DUMMYG=Guess;
                        Error(1,k)=DUMMY(1,k+1);
                        Error(1,k+1)=DUMMY(1,k);

                        for nn = 1:length(Guess(:,1))
                            Guess(nn,k)=DUMMYG(nn,k+1);
                            Guess(nn,k+1)=DUMMYG(nn,k);
                        end

                    end
                end
            end
        end


        for j=1:NB
            % Establish the GUESS vector for this iteration
            GUESS = zeros(length(Guess(:,1)),1);
            for nn = 1:length(Guess(:,1))
                GUESS(nn) = Guess(nn,j);
            end
            
            % Establish a mulitiplier which decreases with iterations
            if ICOUNT<=LIMIT*0.5
                MULT=0.5;
            elseif i>0.5*LIMIT && i<3/4*LIMIT
                MULT=0.15;
            else
                MULT=0.05;
            end

            for nn = 1:length(Guess(:,1))
                % Matt chose to make the multiplier smaller for some of the
                % parameters
                if nn == 1 || 3 || 5
                    GUESS(nn) = (sqrt(GUESS(nn))+MULT*0.1*(0.5-rand))^2;
                else
                    GUESS(nn) = (sqrt(GUESS(nn))+MULT*(0.5-rand))^2;
                end
            end


            Mu0_cNEW=GUESS(1);
            MuINF_cNEW=GUESS(2);
            Sigy0NEW=GUESS(3);
            tr1NEW=GUESS(4);
            tr2NEW=GUESS(5);
            TauCNEW=GUESS(6);
            MuRNEW=GUESS(7);
            aNEW=1;
            dNEW=1/2;


            % Steady State Equations
            LambdaSS=(tr2NEW.*ExpShearSS.^dNEW+1)./(tr1NEW.*ExpShearSS.^aNEW+tr2NEW.*ExpShearSS.^dNEW+1);   %Structure
            SigmaR=LambdaSS.*Sigy0NEW + MuRNEW.*LambdaSS.^m.*ExpShearSS;   %stress contribution from rouloux
            SigmaVISC=((Mu0_cNEW-MuINF_cNEW)./(1+TauCNEW.*ExpShearSS) +MuINF_cNEW).*ExpShearSS;  %Stress contribution from viscosity
            VISC_VE=((Mu0_cNEW-MuINF_cNEW)./(1+TauCNEW.*ExpShearSS) +MuINF_cNEW);
            SigmaTOT=SigmaVISC+SigmaR;  %total stress

            ERROR=0;

            for jj = 1:N
                ERROR = ERROR + ((SigmaTOT(j,1)-ExpStressSS(j,1))/ExpStressSS(j,1))^2;
            end

            s=sqrt(ERROR)/N;

            CURRENTERROR=s;

            if i>1 && i<=ANNEALLIMIT
                ANNEALP=exp((Error(1,j)-CURRENTERROR)/Temp(1,j));
            end

            if CURRENTERROR<ERRORBEST
                ERRORBEST=CURRENTERROR;

                Mu0_cBEST=GUESS(1);
                MuINF_cBEST=GUESS(2);
                Sigy0BEST=GUESS(3);
                tr1BEST=GUESS(4);
                tr2BEST=GUESS(5);
                TauCBEST=GUESS(6);
                MuRBEST=GUESS(7);

                iBEST=i;

                LambdaSS_BEST = LambdaSS;
                SigmaR_BEST = SigmaR;
                SigmaVISC_BEST = SigmaVISC;
                VISC_VE_BEST = VISC_VE;
                SigmaTOT_BEST = SigmaTOT;
            end
            if i>1 && i<=ANNEALLIMIT
                X=rand;
                if CURRENTERROR<Error(1,j)
                    if CURRENTERROR<ERRORBEST
                        ERRORBEST=CURRENTERROR;

                        Mu0_cBEST=GUESS(1);
                        MuINF_cBEST=GUESS(2);
                        Sigy0BEST=GUESS(3);
                        tr1BEST=GUESS(4);
                        tr2BEST=GUESS(5);
                        TauCBEST=GUESS(6);
                        MuRBEST=GUESS(7);

                        iBEST=i;

                        LambdaSS_BEST = LambdaSS;
                        SigmaR_BEST = SigmaR;
                        SigmaVISC_BEST = SigmaVISC;
                        VISC_VE_BEST = VISC_VE;
                        SigmaTOT_BEST = SigmaTOT;
                    end

                    for nn = 1:length(GUESS)
                        Guess(nn,j)=GUESS(nn);
                    end

                    Error(1,j)=CURRENTERROR;
                elseif CURRENTERROR>=Error(1,j) && X< ANNEALP

                    for nn = 1:length(GUESS)
                        Guess(nn,j)=GUESS(nn);
                    end

                    Error(1,j)=CURRENTERROR;
                elseif CURRENTERROR>=Error(1,j) && X>=ANNEALP
                    %Error(1,j)=Error(1,j);

                end
            end
        end

        RESULTS(i,1)=Mu0_cBEST;
        RESULTS(i,2)=MuINF_cBEST;
        RESULTS(i,3)=Sigy0BEST;
        RESULTS(i,4)=tr1BEST;
        RESULTS(i,5)=tr2BEST;
        RESULTS(i,6)=TauCBEST;
        RESULTS(i,7)=MuRBEST;
        RESULTS(i,8)=ERRORBEST;
        RESULTS(i,9)=i;
    end
end

RunTime=cputime-Z

Mu0_cBEST
MuINF_cBEST
Sigy0BEST
tr1BEST
tr2BEST
TauCBEST
MuRBEST
ERRORBEST
i

LambdaSS_BEST
SigmaR_BEST
SigmaVISC_BEST
VISC_VE_BEST
SigmaTOT_BEST


figure(1);
loglog(ExpShearSS, ExpStressSS, 'ro',ExpShearSS, SigmaR_BEST,'b^',ExpShearSS, SigmaVISC_BEST,'g>',ExpShearSS, SigmaTOT_BEST,'k-','MarkerSize',8,'LineWidth',2.);
set(gca,'FontSize',14,'FontWeight','bold','linewidth',2, 'FontName','Times');
legend('DATA','rouleoux stress','visc stress','tot stress','Location','NorthWest');
xlabel('Shear Rate (1/s)');
ylabel('Stress (Pa)');
xlim([min(ExpShearSS) max(ExpShearSS)]);
ylim([.0001 5]);

figure(2);
semilogx(ExpShearSS, LambdaSS_BEST, 'k-.','MarkerSize',8,'LineWidth',2.);
set(gca,'FontSize',14,'FontWeight','bold','linewidth',2, 'FontName','Times');
xlabel('Shear Rate (1/s)');
ylabel('Lambda');
xlim([min(ExpShearSS) max(ExpShearSS)]);
ylim([0 1]);
