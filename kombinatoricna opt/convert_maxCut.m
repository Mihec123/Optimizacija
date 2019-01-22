function [ A,b,C ] = convert_maxCut( e )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

m1 = e(1,2);
n = e(1,1);
m = m1+1;

A=zeros(n);

for i = 1:m1
    A(e(i+1,1),e(i+1,2)) = e(i+1,3);
end

C = 1/4*(diag(A*ones(n,1))-A);


end

