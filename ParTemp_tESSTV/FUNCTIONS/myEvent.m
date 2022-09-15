function [values,isterminal,direction] = myEvent(t,X,tstart)

 values(1) = t;
 %  Don't let integration go for more than 1 seconds.
 values(2) = toc(tstart) < 1;
 isterminal = true(size(values));
 direction = zeros(size(values));
end