function [ a,b,vr ] = separacijskaRavnina( T1,T2 )
%iscemo separacijsko premico y = a*x+b za skupini tock T1 in T2
%ce ta obstaja vrnemo tisto, ki ima najvecjo vertikalno razdaljo do tock
%T1 = {(x1,y1),...,(xn,yn)}, T2 = {(u1,v1),...,(um,vm)}
%Vhod:
% T1 matrika tock nx2 za katere velja T1 > y=a*x+b (lezijo nad premico)
% T2 matrika tock mx2 za katere velja T2 <y = a*x+b (lezijo pod premico)
%Izhod:
% a smerni koeficient premice y = a*x+b
% b odsek na ordinatni osi premice y = a*x+b
% vr razdalja do najblizje tocke

T = size(T1);
n = T(1);
T = size(T2);
m = T(1);

%resujemo max min razdalj
%uvedemo novo spremenljivko z
%kjer velja z <= razdalja v y-osi tock do premic
%problem je potem max z
%p.p premica loci tocke
% z <= razdalja

c = -[1 0 0]; %ker max problem
A = [ones(n,1), T1(:,1) ones(n,1);ones(m,1), -T2(:,1) -ones(m,1)];
b = [T1(:,2);-T2(:,2)];


%premica ima vecjo y vrednost kot T2 in manjso y vrednost kot T1
A = [A; zeros(n,1) T1(:,1) ones(n,1); zeros(m,1) -T2(:,1) -ones(m,1)];
b = [b;T1(:,2);-T2(:,2)];

[X,FVAL,EXITFLAG,OUTPUT] = linprog(c,A,b);
vr = -FVAL;
a = X(2);
b = X(3);

x = linspace(-10,10);
y = a.*x+b;

figure;
hold on;
plot(T1(:,1),T1(:,2),'ro');
plot(T2(:,1),T2(:,2),'ko');
plot(x,y,'b')
hold off;
end

