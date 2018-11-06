function [ a,b,vr ] = separacijskaRavnina( T1,T2 )
%iscemo separacijsko premico y = a*x+b za skupini tock T1 in T2
%ce ta obstaja vrnemo tisto, ki ima najvecjo vertikalno razdaljo do tock
%T1 = {(x1,y1),...,(xn,yn)}, T2 = {(u1,v1),...,(um,vm)}
%smatramo da T1 > y=a*x+b in T2 <y = a*x+b

T = size(T1);
n = T(1);
T = size(T2);
m = T(1);


c = -[1 0 0];
A = [ones(n,1), T1(:,1) ones(n,1);ones(m,1), -T2(:,1) -ones(m,1)];
b = [T1(:,2);-T2(:,2)];

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

