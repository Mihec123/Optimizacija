function [x,vr,y,st] = simpleksMetoda( c,A,b,J,maxstep)
% funkcija izvede postopek simpleksne metode problema min c'x p.p Ax =b,
% x>=0
%Vhod: 
% c nx1 matrika kriterijske funkcije
% A mxn matrika
% b mx1 matika
% J mx1 matrika zacetne bazne resitve, ce ta ni podana jo poiscemo z
%   dvofazno simpleksno metodo
% maxstep najvecje stevilo korakov simpleksne metode, ce ga ne podamo
%   izvedemo najvec 100 korakov
%Izhod:
% x nx1 matriak resitev problema
% vr vrednost kriterijske funkcije pri resitvi x
% y resitev dualnega problema
% st stevilo izvedenih korakov simpleksne metode


if nargin < 5
    maxstep = 100;
end
if nargin < 4
    J = dopustnaResitev(c,A,b)
end


%sledimo navodilom naloge5 iz pdf-ja

V = size(A);
m = V(1);
n = V(2);
K = [1:n];
K(J) = 0;
K(K==0) = [];
x = zeros(length(c),1);

st=0;
while st<maxstep
    xj = A(:,J)\b;
    vr = c(J)'*xj;
    y = A(:,J)'\c(J);
    cpr = c(K)'-c(J)'*inv(A(:,J))*A(:,K);
    if cpr >= 0
        x(J) = xj;
        break;
    else
        temp = cpr;
        [val,indeks] = min(temp);
        ro = K(indeks);
        apr = A(:,J)\A(:,ro);
        if apr <= 0
            vr = inf;
            break;
        end

        v = xj./apr;
        temp = v;
        temp(temp<0) = Inf;
        [val,indeks] = min(temp);
        ro1 = J(indeks);

        %posodobimo J in K
        J(J == ro1) = [];
        J = sort([J ro]);
        K = [1:n];
        K(J) = 0;
        K(K==0) = [];
    end
    st = st+1;

end
    


end

