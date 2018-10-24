b = [52 48 45 45 58 44 45 42 42 38]';
a = [10 26 30 34 39 39 42 46 53 43]';
A = [a ones(length(a),1)];


%prva norma
A1 = zeros(length(A)*2,length(A)+2);
b1 = zeros(length(A)*2,1);
A1(1:2:end,1:end-2) = -eye(length(A));
A1(2:2:end,1:end-2) = -eye(length(A));
b1(1:2:end) = -b;
b1(2:2:end) = b;
A1(1:2:end,end-1) = -A(:,1);
A1(2:2:end,end-1) = A(:,1);
A1(1:2:end,end) = -A(:,2);
A1(2:2:end,end) = A(:,2);
c = ones(1,length(A)+2);
c(end) = 0;
c(end-1) = 0;

x = linprog(c',A1,b1);
k = x(end-1);
n = x(end);
t = linspace(10,53);
v = k*t+n;


%neskoncna norma


c2 = [1 0 0];
A2 = zeros(length(A)*2,3);
b2 = zeros(length(A)*2,1);

A2(:,1)=-ones(length(A)*2,1);
A2(1:2:end,end-1) = -A(:,1);
A2(2:2:end,end-1) = A(:,1);
A2(1:2:end,end) = -A(:,2);
A2(2:2:end,end) = A(:,2);

b2(1:2:end) = -b;
b2(2:2:end) = b;
x1 = linprog(c2,A2,b2);
k2 = x1(2);
n2 = x1(3);

t2 = linspace(10,53);
v2 = k2*t+n2;

hold on;
plot(a,b,'*')
plot(t,v)
plot(t2,v2)
hold off;

