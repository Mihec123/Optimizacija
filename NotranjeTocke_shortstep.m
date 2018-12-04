function [ vrednsot,x,y,s, iter, napaka] = NotranjeTocke_shortstep( c,A,b,x0,y0,sigma,maxit)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
eps = 1e-6;
x=x0;
y=y0;
s=c-A'*y;

vel = size(A);
n = vel(2);

iter = 0;

if nargin < 6
    sigma = 1 - 0.4/sqrt(n);
end


napaka = x'*s;

while napaka > eps
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

vrednsot = c'*x;

end

