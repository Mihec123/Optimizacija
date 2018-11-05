function [X,FVAL,EXITFLAG,OUTPUT,LAMBDA] = separacija( T1,T2 )
%iscemo separacijsko premico y = a*x+b za skupini tock T1 in T2
%ce ta obstaja vrnemo tisto, ki ima najvecjo vertikalno razdaljo do tock
%T1 = {(x1,y1),...,(xn,yn)}, T2 = {(u1,v1),...,(um,vm)}
%smatramo da T1 > y=a*x+b in T2 <y = a*x+b


n = length(T1);
m = length(T2);

c =-[sum(T1(:,1))-sum(T2(:,1)),n-m]; %[(-u1-...-um+x1+...+xn)a,(n-m)b] minus ker resujemo max ne min
A = zeros(n+m,2);
A(1:n,1) = T1(:,1);
A(1:n,2) = ones(n,1);
A(n+1:end,1) = -T2(:,1);
A(n+1:end,2) = -ones(m,1);
b = [T1(:,2);-T2(:,2)];
[X,FVAL,EXITFLAG,OUTPUT,LAMBDA] = linprog(c,A,b);


end

