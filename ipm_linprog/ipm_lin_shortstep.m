function [ vrednost,x,y,s, iter, napaka] = ipm_lin_shortstep( c,A,b,x0,y0,sigma,eps,maxit)
% Opis:
%   ipm_lin_shortstep metoda notranjih tock s kratkim korakom za resevanje
%   optimizacijskega problema
%       min c'x
%       p.p. Ax = b
%            x >= 0
%
% Definicija:
%   [ vrednost,x,y,s, iter, napaka] = ipm_lin_shortstep( c,A,b,x0,y0,sigma,natancnost,maxit)
%
% Vhodni  podatki:
%   c minimizacijska funkcija (min c'x) predstavljena z matriko velikosti n x 1
%   A linearni pogoji ki omejujejo x predstavljeni z matriko
%       velikosti m x n
%   b desna stran sistema Ax = b predstavljena z matriko velikosti m x 1
%   x0 zacetni priblizek za resitev opt. problema velikosti n x 1 veljati mora da je x0 v
%       strogi notranjosti  glede na pogoje
%   y0 zacetni priblizek za resitev dualnega opt. problema velikosti m x 1 veljati mora da je y0 v
%       strogi notranjosti  glede na pogoje duala
%   sigma vrednost, ki doloca premik proti srediscni poti default vrednost sigma
%       je sigma = 1 - 0.4/sqrt(n)
%   eps je vrednost, ki nam pove toleranco napake dobljene resitve
%       default vrednost natancnosti je eps = 1e-6
%   maxit pove maksimalno stevilo iteracij, ki jih izvedemo
%       default vrednost maxit je maxit = 100
%
% Izhodni  podateki:
%   vrednost nam pove vrednost optimizacijskega problema min c'x
%   x matrika n x 1 vrednosti x optimalnega problema 
%   y matrika m x 1 vrednosti y dualnega problema
%   s matrika m x 1 vrednosti dopolnilnih spremenljivk dualnega problema
%   iter stevilo iteracij, ki jih je metoda izvedla
%   napaka velikost napake metode merjene kot x'*s

vel = size(A);
n = vel(2);

if nargin < 6
    sigma = 1 - 0.4/sqrt(n);
end
if nargin < 7
    eps = 1e-6;
end
if nargin < 8
    maxit = 100;
end

x=x0;
y=y0;
s=c-A'*y;



iter = 0;
napaka = x'*s;

while (napaka > eps) && (iter < maxit)
    iter = iter +1;
    tau = (x'*s)/n;
    mu = sigma*tau;
    rp = A*x -b;
    rd = A'*y+s-c;
    rc = x.*s- mu;
    dy = (A*inv(diag(s))*diag(x)*A') \ (rp-A*inv(diag(s))*(rc-diag(x)*rd));
    ds = rd - A'*dy;
    dx = inv(diag(s))*(rc-x.*ds);
    x = x-dx;
    y = y-dy;
    s = s-ds;
    napaka = x'*s;
end

vrednost = c'*x;

end

