function [ vrednsot,x,y,s, iter, napaka] = NotranjeTocke_longstep( c,A,b,x0,y0,sigma,maxit)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

eps = 1e-6;
x=x0;
y=y0;
s=c-A'*y;

vel = size(A);
n = vel(2);

iter = 0;

if nargin < 6
    sigma = 0.5;
end


napaka = x'*s;

while (napaka > eps) && iter < 100 
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
        if prod(double(xtemp > 0)) > 0
            napakax = false;
        else
            alfap = alfap - 0.05;
        end          
    end
    
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
    x = xtemp;
    y = y - alfad*dy;
    s = stemp;
    
    napaka = x'*s;
end

vrednsot = c'*x;

end

