function [best] = stochasticSolver(objective,par0,ITER)

best.TOTerror = 100;
best.SSerror = 1000;
best.TRANSerror = 1000;

fprintf('Prepare for iterations... \n\n')

for i = 1:ITER
    if rem(i,1)==0
        fprintf('Iteration Number: %d\n\n',i)
        fprintf('Steady State Error: %d\n',best.SSerror)
        fprintf('Transient Error: %d\n',best.TRANSerror)
        fprintf('TotalError: %d\n\n',best.TOTerror)
        par0
    end

    if i<=1/2*ITER
        MULT=.15;
    elseif i>1/2*ITER & i <=3/4*ITER
        MULT=.1;
    else
        MULT=.05;
    end

    if i==1  %no change to parameter values during first time through
        par = par0;
    else

        par.mu0 = (sqrt(par0.mu0)+MULT*.025*(.5-rand))^2; % these are 'steady state' par. 
        par.muinf = (sqrt(par0.muinf)+MULT*.01*(.5-rand))^2;
        par.tauC = (sqrt(par0.tauC)+MULT*.01*(.5-rand))^2;

        par.tr1 = (sqrt(par0.tr1)+MULT*.025*(.5-rand))^2;   
        par.tr2 = (sqrt(par0.tr2)+MULT*.015*(.5-rand))^2;
        par.muR = (sqrt(par0.muR)+MULT*.025*(.5-rand))^2;
        par.sigy0 = (sqrt(par0.sigy0)+MULT*.025*(.5-rand))^2;

        par.taulam = (sqrt(par0.taulam)+MULT*(.5-rand))^2;  %these three are transient par
        par.GR = (sqrt(par0.GR)+MULT*(.5-rand))^2;   
        par.GC = (sqrt(par0.GC)+MULT*(.5-rand))^2;
    end

    [obj,pred] = objective(cell2mat(struct2cell(par)));

    if obj < best.TOTerror % Is the total error less than the best so far?
        % Accept new parameters
        par0 = par;

        % Store the parameters, error, and predictions in structure "best"
        best.best_par = par;
        best.TOTerror = obj;
        best.SSerror = pred.ErrorSS;
        best.TRANSerror = pred.ErrorTRANS;
        best.pred = pred;
    end
end

fprintf("\nBest Values... \n")
end