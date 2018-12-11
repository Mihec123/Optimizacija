%vaje 4.12.2018

%nalogi 1,2
A = [1 1 1];
b = 1;
c = [-1 -3 -4];

[ vrednost,x,y,s, iter, napaka] = ipm_lin_shortstep( c',A,b',[1/3 1/3 1/3]',[-5])
[ vrednost,x,y,s, iter, napaka] = ipm_lin_longstep( c',A,b',[1/3 1/3 1/3]',[-5])
[ vrednost,x,y,s, iter, napaka] = ipm_lin_predcor( c',A,b',ones(3,1),ones(1,1))

%naloga 3
randn('state', 0);
rand('state', 0);
n = 500;
m = 400;
c = rand(n,1) + 0.5;
x0 = abs(randn(n,1));
A = abs(randn(m,n));
b = A*x0;

x0 = ones(n,1);
y0 = zeros(m,1);

[ vrednost,x,y,s, iter, napaka] = ipm_lin_predcor( c,A,b,x0,y0);
vrednost
napaka
iter