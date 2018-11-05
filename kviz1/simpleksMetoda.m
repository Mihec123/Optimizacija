function [x,vr,y,st] = simpleksMetoda( c,A,b,J,maxstep)
%resujemo simpleksno metodo v standardni obliki

if nargin < 5
    maxstep = 100;
end
if nargin < 4
    J = dopustnaResitev(c,A,b)
end

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

