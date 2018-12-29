%vaje 3 18.12.2018

C = [1 2 3; 2 9 0; 3 0 7];
A1 = -[1 0 1; 0 3 7; 1 7 5];
A2 = -[0 2 8 ; 2 6 0; 8 0 4];
b = [11;9];
c = C(:);
A = zeros(2,9);
A(1,:) = A1(:)';
A(2,:) = A2(:)';

K.s = 3; % povemo nad gledamo nad stozcem psd matrik dimenzije 3x3
%pars.fid = 0;

[x,y,info] = sedumi(A,-b,c,K);
%podal smo v dualni obliki zato potimalna resitev y ne x;


%naloga2

J = ones(5);
j = J(:);
A1 = eye(5);

A2 = zeros(5);
A2(1,2) = 1;
A2(2,1) = 1;

A3 = zeros(5);
A3(2,3) = 1;
A3(3,2) = 1;

A4 = zeros(5);
A4(3,4) = 1;
A4(4,3) = 1;

A5 = zeros(5);
A5(4,5) = 1;
A5(5,4) = 1;

A6 = zeros(5);
A6(5,1) = 1;
A6(1,5) = 1;

A = zeros(6,25);
A(1,:) = A1(:)';
A(2,:) = A2(:)';
A(3,:) = A3(:)';
A(4,:) = A4(:)';
A(5,:) = A5(:)';
A(6,:) = A6(:)';

%vektor b
b = [1; 0; 0 ; 0; 0 ;0];

%dimenzija 5x5
K.s = 5;

%iscemo max zato -j
[x,y,info] = sedumi(A,b,-j,K);

%resitev je j'*x

%naloga 3

load('podatki_max_cut_sdp.mat')

K.s = 50;
[x,y,info] = sedumi(A,b,c,K);


%naloga 4

%resujemo dual

M = reshape(c,[50,50]);
C = M;
c = C(:);
A1 = eye(length(M));
b = 1;
A = A1(:)';
K.s = length(M);
[x,y,info] = sedumi(A,b,c,K);

