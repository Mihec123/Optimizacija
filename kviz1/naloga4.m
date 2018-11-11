%naloga4
%graf incidencna matrika wikipedija

% G = [1 0 0 0 1 0; 1 1 0 0 0 1; 0 1 1 0 0 0; 0 0 1 1 0 0; 0 0 0 1 1 1];
% %zmisljene
% utezi = [0.5, 0.8, 0.85, 0.8, 0.7, 0.65];
% 
% %[X,vr] = popolnoPrirejanjeCeloSt(utezi,G)
% 
% [X1,vr1] = popolnoPrirejanje(utezi,G)
% 
% %dvodeln graf
% 
% 
% G1 = [1 1 0; 1 0 0; 0 0 1; 0 1 1];
% utezi1 = [0.85,0.9,0.7];
% 
% %[X2,vr2] = popolnoPrirejanjeCeloSt(utezi1,G1)
% 
% [X3,vr3] = popolnoPrirejanje(utezi1,G1)


%primer iz naloge

G = zeros(14,18);

%prva mnozica
G(1,1:3) = ones(1,3);
G(2,4:5) = ones(1,2);
G(3,6:8) = ones(1,3);
G(4,9:11) = ones(1,3);
G(5,12:13) = ones(1,2);
G(6,14:16) = ones(1,3);
G(7,17:18) = ones(1,2);

%druga mnozica
G(8,[1 4 9]) = ones(1,3);
G(9,[6 12 14]) = ones(1,3);
G(10,[2 7 17]) = ones(1,3);
G(11,[5 10]) = ones(1,2);
G(12,[8 15]) = ones(1,2);
G(13,[3 16]) = ones(1,2);
G(15,[11 13 18]) = ones(1,3);

%utezi

utezi = [57 48 95 87 70 96 74 60 75 81 90 55 75 85 26 60 64 88];

[X4,vr4] = popolnoPrirejanje(utezi,G);
find(X4>0.5) % ker ni racunan celo stevilsko rezultati niso cist 0 al pa 1 ampak 0.000001 ,0.99999991 kej tazga