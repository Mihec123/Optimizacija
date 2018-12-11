%kviz2

c = [2 3 3 4]'; %dim c = n x 1 kar je 4 x 1
A = [1 1 1 1; 1 1 2 1; 0 2 1 2];% dim A = m x n kar je  3 x 4
b = [4 5 5]'; % dim b = m x 1 kar je 3 x 1

x0 = ones(4,1);
y0 = 0.5*ones(3,1);
sigma = 0.5;
eps = 1e-6;
faktor = 0.8;
maxit = 100;

[ vrednost,x,y,s, iter, napaka] = ipm_lin_predcor( c,A,b,x0,y0,sigma,faktor,eps,maxit)