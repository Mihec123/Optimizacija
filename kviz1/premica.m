function [w,w1] = premica(tocke)
%funkcija premica isce premico, ki se najbolje prilega danim tockam glede
%na norm(A*w-b,1)+norm(w,inf), kjer je w=(k,n). Problem resujemo z
%linearnim programiranjem
% Vhod:
% tocke matrika nx2 tock
% Izhod:
% w = (k,n) koeficient premice oblike y = kx+n dobljene z norm(A*w-b,1)+norm(w,inf)
% w1 = (k,n) koeficient premice oblike y = kx+n dobljene z norm(A*w-b,2)+norm(w,2)

A = [tocke(:,1),ones(length(tocke),1)];
b = tocke(:,2);

%iscemo min(norm(A*w-b,1)+norm(w,inf)), kjer je w = (k,n)
%min(A*w-b) = min(sum(e_i)), p.p. -e_i \leq a_i'w-b_i \leq e_i

M = [-A, -eye(length(tocke));A,-eye(length(tocke))];
D = [-b;b];
c = [0 0 ones(1,length(tocke))];

%min(w) = min e p.p. -e \leq w-b_i \leq e

%dodamo novo neznanko e
M = [M zeros(length(tocke)*2,1)];
M = [M;-ones(length(tocke),2) zeros(length(tocke)) -ones(length(tocke),1);ones(length(tocke),2) zeros(length(tocke)) -ones(length(tocke),1)];
D = [D;-b;b];
c = [c 1];

X = linprog(c,M,D);

w = X(1:2);

%2 norma
w1 = [A;eye(2)]\[b;zeros(2,1)];

%se narisemo

x = linspace(0,100);
y = w(1).*x+w(2);
y1 = w1(1).*x+w1(2);

hold on;
plot(tocke(:,1),tocke(:,2),'o');
plot(x,y)
plot(x,y1)
hold off;


end

