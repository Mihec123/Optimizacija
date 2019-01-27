function b = operatorA(X,povezave)
[n,m] = size(povezave);
b = zeros(n+1,1);
for i = 1:n
    b(i) = 2*X(povezave(i,1),povezave(i,2));
end
b(end) = sum(diag(X));
end
