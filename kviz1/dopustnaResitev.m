function J = dopustnaResitev(c,A,b)
%poisce dopustno resitev za linearni program ce ta obstaja


Dim = size(A);
m = Dim(1);
n = Dim(2);

A1 = [A eye(m)];
A1
J1 = [n+1:1:n+m];
J1
c1 = [zeros(1,length(c)) ones(1,m)]';
c1
[x,vr,y,st] = simpleksMetoda( c1,A1,b,J1);
if vr ~=0
    display('ni dopustne resitve')
else
    J = find(x>0)';
    if length(J) < m
        [vr,ind] = min(x);
        J = sort([J ind]);
    end
    
end



end

