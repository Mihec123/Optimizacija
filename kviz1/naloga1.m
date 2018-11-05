%naloga1


%tocke so podane v matriki tocke [x y] kjer je x stolpec x koordinat, y
% stolpec y koordinat
tocke = [10 52; 26 48; 30 45; 34 45; 39 58; 39 44; 42 45; 46 42; 53 42; 43 38];

%iscemo premico oblike y = k*x +n

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

[X,FVAL,EXITFLAG,OUTPUT] = linprog(c,M,D);

k=X(1);
n=X(2);

x = linspace(0,100);
y = k.*x+n;

%2 norma
w = [A;eye(2)]\[b;zeros(2,1)];

y1 = w(1).*x+w(2);

hold on;
plot(tocke(:,1),tocke(:,2),'o');
plot(x,y)
plot(x,y1)
hold off;