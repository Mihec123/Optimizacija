function [ Xk,y,D1,Zk,vrednost,iter,err_p,err_d ] = okrepljena_lagrangeva_metoda_lovazs( povezave,b,C,sigma,itermax,napaka )
% Opis:
%   okrepljena_lagrangeva_metoda_lovazs je okrepljena lagranzeva metoda za
%   dual izbolsane theta funkcije
%       min <J,X>
%       p.p. tr(X) = 1
%       X_ij = 0 za ij in E
%       X => 0
%       X s.p.d
%
% Definicija:
%   [ Xk,y,D1,Zk,vrednost,iter,err_p,err_d ] = okrepljena_lagrangeva_metoda_lovazs( povezave,b,C,sigma,itermax,napaka )
%
% Vhodni  podatki:
%   povezave matrika m x 2, ki pove katere tocke so povezane
%   b matrika m x 1 ki doloca desno stran pogojev v sdp programu
%   C matrika n x n, katero optimiziramo pri sdp programu
%   sigma vrednost, ki doloca premik proti srediscni poti default vrednost sigma
%       je sigma = 0.5
%   napaka je vrednost, ki nam pove toleranco napake dobljene resitve
%       default vrednost natancnosti je eps = 1e-6
%   maxit pove maksimalno stevilo iteracij, ki jih izvedemo
%       default vrednost maxit je maxit = 2000
%
% Izhodni  podateki:
%   vrednost nam pove vrednost optimizacijskega problema min 0.5x'Px + q'x
%   Xk matrika n x n vrednosti x optimalnega problema 
%   iter stevilo iteracij, ki jih je metoda izvedla
%   err_p napaka primarnega problema
%   err_d napaka dualnega problema

if nargin < 6
    itermax = 500;
end
if nargin < 7
    napaka = 1e-6;
end
[n,m] = size(C);
Xk = zeros(n); %lahko eksperimentiras
D1 = zeros(n); %ta uredu, lah bi blo tud kej druzga

iter = 0;
err_p = Inf;
err_d = Inf;

y = rand(length(b),1);
delimo = 2*ones(length(b),1);
delimo(end) =  n;


while (iter < itermax) && ((err_p > napaka) || (err_d > napaka))
    W = operatorAt(y,povezave,n) - C -D1- 1/sigma*Xk;
    %W = operatorAt(y,povezave,n) - C - 1/sigma*Xk;
    [V,D] = eig(W);
    plus = max(D,0);
    minus = min(D,0);
    Zk = V*plus*V';
    %test = operatorAt(y,povezave,n)-C;
    %Zk(test-Zk<=0) = test(test-Zk<=0);
    Xk = -sigma*V*minus*V';
    D1 = max(0,-C-Zk+operatorAt(y,povezave,n)-1/sigma*Xk);
    %specialno za stable set
    y = (operatorA(Zk+C+D1,povezave)+1/sigma*(operatorA(Xk,povezave)-b))./delimo;
    %y = (operatorA(Zk+C,povezave)+1/sigma*(operatorA(Xk,povezave)-b))./delimo;

    
    iter = iter+1;
    
    err_p = norm(operatorA(Xk,povezave)-b);
    err_d = norm(Zk-operatorAt(y,povezave,n)+C+D1);
    %err_d = norm(Zk-operatorAt(y,povezave,n)+C);
    if mod(iter,1) == 0
        vrednost = b'*y;
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
vrednost = b'*y;
end
