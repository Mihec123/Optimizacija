function [X,vr] = popolnoPrirejanjeCeloSt(utezi,Graf)
%Graf je incidenca matrika
c = -utezi;
A = size(Graf);
n = A(1);
m = length(utezi);

b = ones(n,1);

[X,FVAL] = intlinprog(c,[1:m],Graf,b,[],[],zeros(m,1),ones(m,1));
vr = -FVAL; 

end

