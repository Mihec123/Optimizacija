function [ r,s] = najkrog( premice1,premice2 )
%vrnemo krog z najvecjim radijem r, ki je dolocen s politopom dolocen z
%premicami.
%Vhod:
%premice1 matrika n1x2 [k n], kjer smatramo da je premica oblike y<=kx+n
%premice2 matrika n2x2 [k n], kjer smatramo da je premica oblike y>=kx+n
%Izhod:
% r radij najvecjega kroga v tem politopu
% s=(s1,s2) koordinati sredisca kroga


%lotimo se tako da izracunamo kot pravokotnice na vsako premico in potem
%pri vsakem takem kotu pogledamo omejitev te premice
% (s1+rcos(fi),s2+rsin(fi)) <= y=ax+b (za prvi kvadrant)

pravkot1 = -1./premice1(:,1);
pravkot2 = -1./premice2(:,1);
alfa1 = atan(pravkot1(:,1));
alfa2 = atan(pravkot2(:,1));

%iscemo max radij
c = -[1 0 0]; %max [r s1 s2]
%pogoji

A=zeros(length(premice1)+length(premice2),3);
b = zeros(length(premice1)+length(premice2),1);
for i=1:length(premice1)
    k = premice1(i,1);
    if k<0        
        A(i,:) = [sin(alfa1(i))-cos(alfa1(i))*k, -k, 1]; %[r(sin(fi)-k*cos(fi)),-k*s1,s2]
        b(i) = premice1(i,2); %n
    else
        A(i,:) = [sin(pi+alfa1(i))-cos(pi+alfa1(i))*k, -k, 1];
        b(i) = premice1(i,2); %n                
    end
end
n = length(premice1);
for j = 1:length(premice2)
    k = premice2(j,1);
    if k >0
        radtodeg(alfa2(j))
        %radtodeg(alfa2(j))
        A(n+j,:) = [-sin(alfa2(j))+cos(alfa2(j))*k, k, -1]; %[r(k*cos(fi)-sin(fi)),k*s1,-s2] + fi treba mert iz prave strani
        b(n+j) = -premice2(j,2); %n
    else
        A(n+j,:) = [-sin(-(pi-alfa2(j)))+cos(-(pi-alfa2(j)))*k, k, -1]; %[r(k*cos(fi)-sin(fi)),k*s1,-s2] + fi treba mert iz prave strani
        b(n+j) = -premice2(j,2); %n        
    end
    
end

[X,FVAL,EXITFLAG] = linprog(c,A,b);
r = X(1);
s = X(2:3);
s1 = X(2);
s2 = X(3);

%se narisemo (%slika ni ok razmerja osi so poroblem)
x = linspace(s1-r-3,s1+r+3);
figure;
axis([s1-r-3, s1+r+3, s2-r-3, s2+r+3 ])
hold on;
for i = 1:length(premice1)
    y = premice1(i,1).*x + premice1(i,2);
    plot(x,y,'b')
    pbaspect([1 1 1])
end

for j = 1:length(premice2)
    y = premice2(j,1).*x + premice2(j,2);
    plot(x,y,'b')
    pbaspect([1 1 1])
end
%se krog
th = 0:pi/50:2*pi;
xunit = r * cos(th) + s1;
yunit = r * sin(th) + s2;
plot(xunit, yunit,'k');
pbaspect([1 1 1])
end

