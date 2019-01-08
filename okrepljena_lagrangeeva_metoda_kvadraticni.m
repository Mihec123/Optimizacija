function [ x,vrednost,err_p,err_d,iter] = okrepljena_lagrangeeva_metoda_kvadraticni( P,q,A,b,sigma,itermax,napaka )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin < 6
    itermax = 300;
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

while (iter < itermax) && ((err_p > napaka) || (err_d > napaka))
    
    x = (P+sigma.*A'*A) \ ((-q-A'*lambda)-sigma.*(A'*(s-b))); %bol z razcepom choleskega
    %z razcepom choleskega
    
    %DN
    
    s = max(-A*x+b-lambda/sigma,0);
    lambda = lambda + sigma*(A*x+s-b);
    err_p = abs(max(b-A*x-s));
    err_d = abs(max(P*x+q+A'*lambda));
    iter = iter+1;
    
    
end

vrednost = 0.5*x'*P*x+q'*x;
        




end

