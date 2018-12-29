%test ipm_sdp_predcor

%m1,m,n so nastimane kot v opisu funkcije

m1 = 5;
m = m1^2;
n = 6;


C = ones(m1);
c = C(:);

A1 = eye(m1);

A2 = zeros(m1);
A2(1,2) = 1;
A2(2,1) = 1;

A3 = zeros(m1);
A3(2,3) = 1;
A3(3,2) = 1;

A4 = zeros(m1);
A4(3,4) = 1;
A4(4,3) = 1;

A5 = zeros(m1);
A5(4,5) = 1;
A5(5,4) = 1;

A6 = zeros(m1);
A6(5,1) = 1;
A6(1,5) = 1;

A = zeros(n,m);
A(1,:) = A1(:)';
A(2,:) = A2(:)';
A(3,:) = A3(:)';
A(4,:) = A4(:)';
A(5,:) = A5(:)';
A(6,:) = A6(:)';

b = zeros(n,1);
b(1) = 1;

X0 = eye(m1);
y0 = zeros(n,1);
y0(1) = -6;

[ vrednost,X,y,Z, iter, napaka] = ipm_sdp_predcor( -c,A,b,X0,y0); %maxsimiziramo zato -c