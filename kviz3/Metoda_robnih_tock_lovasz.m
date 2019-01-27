function [ Xk,y,vrednost,iter,err_p,err_d ] = Metoda_robnih_tock_lovasz( povezave,C,b,sigma,itermax,napaka )
% Opis:
%   Metoda_robnih_tock_lovasz metoda ki racuna vrednost lovazseve theta
%   funkcije
%   optimizacijski problem
%       max <J,X>
%       p.p tr(X) = 1
%           X_ij = 0 za ij in E
%           X s.p.d
% Definicija:
%   [ Xk,y,vrednost,iter,err_p,err_d ] = Metoda_robnih_tock_lovasz( povezave,C,b,sigma,itermax,napaka )
%
% Vhodni  podatki:
%   C minimizacijska m1 x m1 matrika
%   povezave matrika m x 2, ki pove katere tocke so povezane
%   b desna stran sistema A(X) = b predstavljena z matriko velikosti n x 1
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
%   err_p velikost napake primarnega programa
%   err_d velikost napake dualnega programa

if nargin < 4
    sigma = 0.5;
end

if nargin < 5
    itermax = 2100;
end
if nargin < 6
    napaka = 1e-6;
end

disp(itermax)
iter = 0;
err_p = Inf;
err_d = Inf;

[n,~] = size(C);
Xk = zeros(n);

y = rand(length(b),1);
delimo = 2*ones(length(b),1);
delimo(end) =  n;

operatory = operatorAt(y,povezave,n);

while (iter < itermax) && ((err_p > napaka) || (err_d > napaka))
    
    W = operatory - C - 1/sigma*Xk;
    [V,D] = eig(W);
    plus = max(D,0);
    minus = min(D,0);
    Zk = V*plus*V';
    Xk = -sigma*V*minus*V';
    %specialno za stable set
    opA = operatorA(Zk+C,povezave);
    opA1 = operatorA(Xk,povezave);
    y = (opA+1/sigma*(opA1-b))./delimo;
    
    operatory = operatorAt(y,povezave,n);

    
    iter = iter+1;
    
    err_p = norm(opA1-b);
    err_d = norm(Zk-operatory+C);
    if mod(iter,100) == 0
        vrednost = C(:)'*Xk(:);
        text = strcat('iter: ',int2str(iter));
        text = strcat(text,'  err_p: ');
        text = strcat(text,num2str(err_p));
        text = strcat(text,'  err_d: ');
        text = strcat(text,num2str(err_d));
        text = strcat(text,'  vrednost: ');
        text = strcat(text,num2str(vrednost));
        disp(text)
    end
end
vrednost = C(:)'*Xk(:);
end




