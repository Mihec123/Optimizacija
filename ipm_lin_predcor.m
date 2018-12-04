function [ vrednost,x,y,s, iter, napaka] = ipm_lin_predcor( c,A,b,x0,s0,sigma,maxit)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



eps = 1e-6;
x=x0;
s = s0;
y = A' \(c-s);

vel = size(A);
n = vel(2);

if nargin < 6
    sigma = 0.5;
end


%prvi korak predikotor

rp = b-A*x;
rd = c - A'*y-s;
rc = - x.*s;
X = diag(x);
S = diag(s);

dy1 = (A*inv(S)*X*A') \ (rp-A*inv(S)*(rc-x.*rd));
ds1 = rd - A'*dy1;
dx1 = S \ (rc-X*ds1);

% x = x +dx1;
% y = y+dy1;
% s = s+ds1;

%drugi korak korektor

tau = (x'*s)/n;
mu = sigma*tau;
% rp = b-A*x;
% rd = c - A'*y-s;
rcc = mu*ones(n,1) - x.*s - dx1.*ds1;
% X = diag(x);
% S = diag(s);

dy = (A*inv(S)*X*A') \ (rp-A*inv(S)*(rcc-x.*rd));
ds = rd - A'*dy;
dx = S \ (rcc-X*ds);


alfap = 1;
napakax = true;
while napakax
    xtemp = x+alfap*dx;      
    if prod(double(xtemp > 0)) > 0
        napakax = false;
    else
        alfap = alfap - 0.05;
    end          
end

alfad = 1;
napakas = true;
while napakas
    stemp = s+alfad*ds;      
    if prod(double(stemp > 0)) > 0
        napakas = false;
    else
        alfad = alfad - 0.05;
    end          
end

x = xtemp;
y = y + alfad*dy;
s = stemp;

[ vrednost,x,y,s, iter, napaka] = NotranjeTocke_longstep( c,A,b,x,y)




end

