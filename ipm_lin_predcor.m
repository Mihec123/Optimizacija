function [ vrednost,x,y,s, iter, napaka] = ipm_lin_predcor( c,A,b,x0,y0,sigma,faktor,eps,maxit)
% Opis:
%   ipm_lin_predcor metoda notranjih tock prediktor korektor
%   optimizacijskega problema
%       min c'x
%       p.p. Ax = b
%            x >= 0
%
% Definicija:
%   [ vrednost,x,y,s, iter, napaka] = ipm_lin_predcor( c,A,b,x0,y0,sigma,natancnost,maxit)
%
% Vhodni  podatki:
%   c minimizacijska funkcija (min c'x) predstavljena z matriko velikosti n x 1
%   A linearni pogoji ki omejujejo x predstavljeni z matriko
%       velikosti m x n
%   b desna stran sistema Ax = b predstavljena z matriko velikosti m x 1
%   x0 zacetni priblizek za resitev opt. problema velikosti n x 1 veljati
%   mora da je x0 strogo pozitiven
%   y0 zacetni priblizek za resitev dualnega opt. problema velikosti m x 1
%       veljati mora da je y0 pozitiven
%   sigma vrednost, ki doloca premik proti srediscni poti default vrednost sigma
%       je sigma = 0.5
%   faktor je vrednost med 0<faktor<1, ki pove s kaksnim faktorjem manjsamo iskanje
%       alfe pri newtnovi metodi v smeri x in s. Default value = 0.8
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
    sigma = 0.5;
end
if nargin < 7
    faktor = 0.8;
end
if nargin < 8
    eps = 1e-6;
end
if nargin < 9
    maxit = 100;
end

x=x0;
y = y0;
s = c - A'*y0;

%prvi korak predikotor

rp = b-A*x;
rd = c - A'*y-s;
rc = - x.*s;
X = diag(x);
S = diag(s);

dy1 = (A*inv(S)*X*A') \ (rp-A*inv(S)*(rc-x.*rd));
ds1 = rd - A'*dy1;
dx1 = S \ (rc-X*ds1);

%drugi korak korektor

tau = (x'*s)/n;
mu = sigma*tau;
rcc = mu*ones(n,1) - x.*s - dx1.*ds1;

dy = (A*inv(S)*X*A') \ (rp-A*inv(S)*(rcc-x.*rd));
ds = rd - A'*dy;
dx = S \ (rcc-X*ds);


%pogledamo kako velik korak lahko naredimo v smeri x da smo
%se dopustni
alfap = 1;
napakax = true;
while napakax
    xtemp = x+alfap*dx;      
    if prod(double(xtemp > 0)) > 0
        napakax = false;
    else
        alfap = faktor * alfap;
    end          
end

%pogledamo kako velik korak lahko naredimo v smeri s da smo
%se dopustni
alfad = 1;
napakas = true;
while napakas
    stemp = s+alfad*ds;      
    if prod(double(stemp > 0)) > 0
        napakas = false;
    else
        alfad = faktor * alfad;
    end          
end

x = xtemp;
y = y + alfad*dy;
s = stemp;


%dobili smo dopustne resitve v strogi notranjosti zato lahko pozenemo
%metodo ipm_lin_longstep na x,y

[ vrednost,x,y,s, iter, napaka] = ipm_lin_longstep( c,A,b,x,y,sigma,faktor,eps,maxit);
iter = iter+1; %kot en korak stejemo se predikotor in korektor



end

