function [ Xk,y,vrednost,iter,err_p,err_d ] = Metoda_robnih_tock_sdp( A,C,b,y0,sigma,itermax,napaka )
% Opis:
%   Metoda_robnih_tock_sdp metoda robnih tock za
%   optimizacijski problem
%       min <C,X>
%       p.p. A(X) = b
%            X >= 0 (X p.s.d.)
%
%       kjer je A(X) = [<A1,X>; ...;<An,X>] in A'(y) = sum (y_i.*Ai)
% Definicija:
%   [ Xk,y,vrednost,iter,err_p,err_d ] = Metoda_robnih_tock_sdp( A,c,b,y,X0,sigma,itermax,napaka )
%
% Vhodni  podatki:
%   c minimizacijska funkcija (min c'*X(:)) predstavljena z matriko
%       velikosti m x 1, kjer je c = C(:) in C neka m1 x m1 matrika torej
%       za m velja m = m1*m1
%   A matrika velikosti n x m kjer vsaka vrstica predstavlja matriko Ai(:)'
%       in velja da je matrika Ai velikosti m1 x m1
%   b desna stran sistema A(X) = b predstavljena z matriko velikosti n x 1
%   X0 zacetni priblizek za resitev opt. problema velikosti m1 x m1 veljati
%       mora da je X0 strogo pozitivno definitna
%   y0 zacetni priblizek za resitev dualnega opt. problema velikosti n x 1
%       veljati mora da je y0 pozitiven
%   sigma vrednost, ki doloca premik proti srediscni poti default vrednost sigma
%       je sigma = 0.5
%   napaka je vrednost, ki nam pove toleranco napake dobljene resitve
%       default vrednost natancnosti je eps = 1e-6
%   maxit pove maksimalno stevilo iteracij, ki jih izvedemo
%       default vrednost maxit je maxit = 500
%
% Izhodni  podateki:
%   vrednost nam pove vrednost optimizacijskega problema min <C,X>
%   X matrika m1 x m1 vrednosti X optimalnega problema 
%   y matrika n x 1 vrednosti y dualnega problema
%   Z matrika m1 x m1 vrednosti dopolnilnih spremenljivk dualnega problema
%   iter stevilo iteracij, ki jih je metoda izvedla
%   napaka velikost napake metode merjene kot X(:)'*Z(:)

if nargin < 6
    sigma = 0.5;
end

if nargin < 7
    itermax = 2000;
end
if nargin < 8
    napaka = 1e-6;
end

iter = 0;
err_p = Inf;
err_d = Inf;

[n,~] = size(C);
Xk = zeros(n);

y = y0;

V1 = chol(A*A','lower');

while (iter < itermax) && ((err_p > napaka) || (err_d > napaka))
    W = reshape(A'*y,[n,n]) - C - 1/sigma*Xk;
    [V,D] = eig(W);
    plus = max(D,0);
    minus = min(D,0);
    Zk = V*plus*inv(V);
    Xk = -sigma*V*minus*inv(V);
    %normalno
    %y = (A*A')\(A*(Zk(:)+C(:))+1/sigma*(A*Xk(:)-b));
    
    
    %nardimo choleskega
    y1 = V1\(A*(Zk(:)+C(:))+1/sigma*(A*Xk(:)-b));
    y = V1'\y1;
    
    iter = iter+1;
    
    err_p = norm(A*Xk(:)-b);
    err_d = norm(Zk-reshape(A'*y,[n,n])+C);
    if mod(iter,10) == 0
        text = strcat('iter: ',int2str(iter));
        text = strcat(text,'  err_p: ');
        text = strcat(text,num2str(err_p));
        text = strcat(text,'  err_d: ');
        text = strcat(text,num2str(err_d));
        disp(text)
    end
end
vrednost = C(:)'*Xk(:);
end

