function [X,vr] = popolnoPrirejanje(utezi,Graf)
%funkcija vrne resitev relaksiranega problema popolnega prirejanje grafa Graf z utezmi na povezavah
%Vhod:
% utezi 1xn matrika utezi povezav
% Graf je podan z incidencno matriko mxn ,kjer je n stevilo povezav, m
%       stevilo vozlisc. Graf(i,j) = 1, ce je krajisce povezave j vozlisce
%       i in 0 sicer
%Izhod:
% X matrika 1xn, ki nam pove deleze povezav, ki so v prirejanju

%za dodaten opis glej https://www.imsc.res.in/~meena/matching/lecture5.pdf

c = -utezi;%resujemo max problem
A = size(Graf);
n = A(1);
m = length(utezi);

b = ones(n,1);

[X,FVAL] = linprog(c,Graf,b,[],[],zeros(m,1),ones(m,1));
vr = -FVAL; 

end

