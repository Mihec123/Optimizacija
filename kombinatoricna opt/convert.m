function [ A,b,C ] = convert( e )
%pretvori tabelo grafa v podatke za stabilno množico

m1 = e(1,2);
n = e(1,1);
m = (m1+1);

C = ones(n);
b = zeros(m,1);
b(m) = 1;
A = zeros(m,n^2);
Id = eye(n);
A(m,:) = Id(:)';

for k=1:m1
    Ai = zeros(n);
    Ai(e(k+1,1),e(k+1,2))=1;
    Ai(e(k+1,2),e(k+1,1))=1;
    
    A(k,:) = Ai(:)';
end




end

