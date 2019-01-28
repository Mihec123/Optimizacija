
%1.naloga a primer
% e = dlmread('brock400-1_den75.dat');
% [ A,b,C ] = convert_stable( e );
% sigma = 0.05/e(1);
% 
% [ Xk,y,vrednost,iter,err_p,err_d ] = Metoda_robnih_tock_lovasz( A,C,b,sigma,15000);
% [ Xk,y,D1,Zk,vrednost,iter,err_p,err_d ] = okrepljena_lagrangeva_metoda_lovazs( A,b,C,sigma);

%1.naloga b primer

% e = dlmread('keller4_clq.dat');
% [ A,b,C ] = convert_stable( e );
% sigma = 0.05/e(1);
% %[ Xk,y,vrednost,iter,err_p,err_d ] = Metoda_robnih_tock_lovasz( A,C,b,sigma);
% [ Xk,y,D1,Zk,vrednost,iter,err_p,err_d ] = okrepljena_lagrangeva_metoda_lovazs( A,b,C,sigma);


%2.naloga

% e = dlmread('keller4_clq.dat');
% [ S,vrednost,Xk ] = hevristika( e );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MAX CUT

%3.naloga

e = dlmread('graph_max_cut.dat');
%e = dlmread('brock400-1_den75.dat');
[ b,C ] = convert_maxCut( e );
sigma = 0.05/e(1);

[ vrednost,X,y,Z, iter, napaka] = ipm_sdp_predcor_maxCut( -C,b,eye(length(C)),-300*ones(e(1),1));

%4.naloga 

% e = dlmread('graph_max_cut.dat');
% n = e(1);
% opt = Goemans_Williamson( e,eye(n),-300*ones(n,1));