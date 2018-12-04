function [ vrednost,x,y,s, iter, napaka] = ipm_lin_longstep( c,A,b,x0,y0,sigma,eps,maxit)
% Opis:
%   ipm_lin_longstep metoda notranjih tock z dolgim korakom za resevanje
%   optimizacijskega problema
%       min c'x
%       p.p. Ax = b
%            x >= 0
%
% Definicija:
%   [ vrednost,x,y,s, iter, napaka] = ipm_lin_longstep( c,A,b,x0,y0,sigma,natancnost,maxit)
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
%       je sigma = 0.5
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
    
    
    %pogledamo kako velik korak lahko naredimo v smeri x da smo
    %se dopustni
    alfap = 1;
    napakax = true;
    while napakax
        xtemp = x-alfap*dx;      
        if prod(double(xtemp > 0)) > 0
            napakax = false;
        else
            alfap = alfap - 0.05;
        end          
    end
    
    %pogledamo kako velik korak lahko naredimo v smeri s da smo
    %se dopustni
    alfad = 1;
    napakas = true;
    while napakas
        stemp = s-alfad*ds;      
        if prod(double(stemp > 0)) > 0
            napakas = false;
        else
            alfad = alfad - 0.05;
        end          
    end
    
    %za korak v y vzamemo kar velikost koraka v s
    x = xtemp;
    y = y - alfad*dy;
    s = stemp;
    
    napaka = x'*s;
end

vrednost = c'*x;

end

