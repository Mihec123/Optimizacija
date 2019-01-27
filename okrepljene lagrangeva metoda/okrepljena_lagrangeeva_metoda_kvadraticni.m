function [ x,vrednost,err_p,err_d,iter] = okrepljena_lagrangeeva_metoda_kvadraticni( P,q,A,b,sigma,itermax,napaka )
% Opis:
%   okrepljena_lagrangeeva_metoda_kvadraticni je okrepljena lagranzeva
%   metoda za kvadraticne probleme
%       min 0.5x'Px + q'x
%       p.p. A(X) <= b
%
% Definicija:
%   [ x,vrednost,err_p,err_d,iter] = okrepljena_lagrangeeva_metoda_kvadraticni( P,q,A,b,sigma,itermax,napaka )
%
% Vhodni  podatki:
%   P  pozitivno definitna n x n matrika
%   q  vektor velikosti n x 1
%   A matrika velikosti m x n 
%   b matrika velikosti m x 1
%   sigma vrednost, ki doloca premik proti srediscni poti default vrednost sigma
%       je sigma = 0.5
%   napaka je vrednost, ki nam pove toleranco napake dobljene resitve
%       default vrednost natancnosti je eps = 1e-6
%   maxit pove maksimalno stevilo iteracij, ki jih izvedemo
%       default vrednost maxit je maxit = 2000
%
% Izhodni  podateki:
%   vrednost nam pove vrednost optimizacijskega problema min 0.5x'Px + q'x
%   x matrika n x 1 vrednosti x optimalnega problema 
%   iter stevilo iteracij, ki jih je metoda izvedla
%   err_p napaka primarnega problema
%   err_d napaka dualnega problema

if nargin < 6
    itermax = 2000;
end
if nargin < 7
    napaka = 1e-6;
end
[n,m] = size(A);
lambda = zeros(n,1); %lahko eksperimentiras
s = zeros(length(b),1); %ta uredu, lah bi blo tud kej druzga

iter = 0;
err_p = Inf;
err_d = Inf;

V1 = chol(P+sigma.*A'*A,'lower');

while (iter < itermax) && ((err_p > napaka) || (err_d > napaka))
    
    %normalno
    %x = (P+sigma.*A'*A) \ ((-q-A'*lambda)-sigma.*(A'*(s-b))); %bol z razcepom choleskega
    
    %z razcepom choleskega
    x1 = V1\((-q-A'*lambda)-sigma.*(A'*(s-b)));
    x = V1'\x1;
    
    
    s = max(-A*x+b-lambda/sigma,0);
    lambda = lambda + sigma*(A*x+s-b);
    err_p = abs(max(b-A*x-s));
    err_d = abs(max(P*x+q+A'*lambda));
    iter = iter+1;
    
    
    if mod(iter,10) == 0
        text = strcat('iter: ',int2str(iter));
        text = strcat(text,'  err_p: ');
        text = strcat(text,num2str(err_p));
        text = strcat(text,'  err_d: ');
        text = strcat(text,num2str(err_d));
        disp(text)
    end
    
    
end

vrednost = 0.5*x'*P*x+q'*x;
        




end

