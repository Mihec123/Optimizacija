function [x,vr,y,st] = simpleksMetoda( c,A,b,J,maxstep)
%resujemo simpleksno metodo v standardni obliki

if nargin < 5
    maxstep = 100;
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
    K
    J
    xj = A(:,J)\b;
    vr = c(J)'*xj;
    y = A(:,J)'\c(J);
    cpr = c(K)'-c(J)'*inv(A(:,J))*A(:,K);
    if cpr > 0
        x(J) = xj;
        break;
    else
        temp = cpr<0;
        temp1 = find(temp > 0);
        ro = K(temp1(1));
        apr = A(:,J)\A(:,ro);
        apr
        xj
        if apr <= 0
            vr = inf;
            break;
        end

        v = xj./apr;
        v
        v(v == Inf) = -1;
        temp = v >0;
        temp1 = find(temp > 0);
        ro1 = J(temp1(1));

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

