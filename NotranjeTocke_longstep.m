function [ vrednsot,x,y,s, iter, napaka] = NotranjeTocke_longstep( x0,y0,c,A,b,sigma)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

eps = 1e-6;
x=x0;
y=y0;
s=c-A'*y;

vel = size(A);
n = vel(1);

iter = 0;

if nargin < 6
    sigma = 0.5;
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
    
    alfap = 1;
    napakax = true;
    while napakax
        xtemp = x-alfap*dx;
        ytemp = y-alfap*dy;
        stemp = s-alfap*ds;
        
        if prod(double(xtemp > 0)) > 0  && prod(double(stemp > 0)) >0
            napakax = false;
        else
            alfap = alfap - 0.1;
        end          
    end
    x = xtemp;
    y = ytemp;
    s = stemp;
    
    napaka = x'*s;
    napaka
end

vrednsot = c'*x;

end

