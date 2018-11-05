%naloga2

T1 = [1,1;2,1;1.5,2;5,3];
T2 = [2,5;3,6;1,5;2,7];

[X,FVAL,EXITFLAG,OUTPUT,LAMBDA] =separacija( T2,T1 );

a = X(1);
b = X(2);
x = linspace(-10,10);
y = a.*x+b;

figure;
hold on;
plot(T1(:,1),T1(:,2),'ro');
plot(T2(:,1),T2(:,2),'ko');
plot(x,y,'b')
hold off;