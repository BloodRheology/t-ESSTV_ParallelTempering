function [obj,pred] = SteadyState_OBJ(parVec,exp)

at=1;
d= 1/2;
m= 1.5;

par.mu0 = parVec(1);
par.muinf = parVec(2);
par.tauC = parVec(3);

par.tr1 = parVec(4);
par.tr2 = parVec(5);
par.muR = parVec(6);
par.sigy0 = parVec(7);

% par.taulam = parVec(8);
% par.GR = parVec(9);
% par.GC = parVec(10);

N = length(exp.ShearSS);

LambdaSS=(par.tr1.*exp.ShearSS.^d+1)./(par.tr1.*exp.ShearSS.^at+par.tr2.*exp.ShearSS.^d+1);   %Structure
SigmaR=LambdaSS.*par.sigy0 + par.muR.*LambdaSS.^m.*exp.ShearSS;   %stress contribution from rouleaux
SigmaVISC=((par.mu0-par.muinf)./(1+par.tauC.*exp.ShearSS) + par.muinf).*exp.ShearSS;  %Stress contribution from viscosity
VISC_VE=((par.mu0 - par.muinf)./(1+par.tauC.*exp.ShearSS) + par.muinf);
SigmaTOT=SigmaVISC+SigmaR;  %total stress


ERROR = sum(((SigmaTOT-exp.StressSS)./exp.StressSS).^2);

s=sqrt(ERROR)/N;

% if LambdaSS(1) > 1
%     s = 1000;
% end


obj = s;

pred.LambdaSS = LambdaSS;
pred.SigmaR = SigmaR;
pred.SigmaVISC = SigmaVISC;
pred.VISC_VE = VISC_VE;
pred.SigmaTOT = SigmaTOT;
end