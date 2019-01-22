addpath('./data')
e = dlmread('HoG_6575.dat');
%poracunamo z metodo robnih tock
[ A,b,C ] = convert( e );
X0 = eye(length(C));
y0 = rand(length(b),1);
tic;
[ Xk,y,vrednost,iter,err_p,err_d ] = Metoda_robnih_tock_sdp( A,C,b,y0,0.5);
toc