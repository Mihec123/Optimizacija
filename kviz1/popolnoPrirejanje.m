function [X,vr] = popolnoPrirejanje(utezi,Graf)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
c = -utezi;
A = size(Graf);
n = A(1);
m = length(utezi);

b = ones(n,1);

[X,FVAL] = linprog(c,Graf,b,[],[],zeros(m,1),ones(m,1));
vr = -FVAL; 

end

