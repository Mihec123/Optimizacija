addpath('./data')
e = dlmread('HoG_6575.dat');
%poracunamo z metodo robnih tock
[ A,b,C ] = convert_stable( e );
X0 = eye(length(C));
y0 = rand(length(b),1);
% tic;
% [ Xk,y,vrednost,iter,err_p,err_d ] = Metoda_robnih_tock_sdp( A,C,b,y0,0.5);
% toc
% 
% 
% %max_cut problem
% e = [5 6 0;1 2 1; 2 3 1;3 4 1; 4 5 1; 5 1 1; 5 2 1 ];
% [ A,b,C ] = convert_maxCut( e );