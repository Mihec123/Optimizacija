%naloga 5

c = [-3; -1; -3; 0; 0; 0];
A = [2 1 1 1 0 0; 1 2 3 0 1 0; 2 2 1 0 0 1];
b = [2;5;6];
J = [4,5,6];
% c = [0;0 ;0; 0; 0; 1; 1];
% A = [0 1 2 0 0 1 0;-1 4 2 1 0 0 0; 2 1 3 0 -1 0 1];
% b = [44;110;132];
% J = [4 6 7];
[x,vr,y,st] = simpleksMetoda( c,A,b,J);
x
linprog(c,[],[],A,b,zeros(length(c),1))

%drug primer
c = [-1 -1 0 0 0]';
A = [1 0 2 1 0;0 1 -1 0 1; 1 1 1 0 0];
b = [1;1;2];
[x,vr,y,st] = simpleksMetoda( c,A,b,[],100,2);
x
linprog(c,[],[],A,b,zeros(length(c),1))