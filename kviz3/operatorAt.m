function [Y] = operatorAt( y,povezave,n )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

Y = zeros(n);
for i = 1:length(povezave)
    Y(povezave(i,1),povezave(i,2)) = y(i);
    Y(povezave(i,2),povezave(i,1)) = y(i);
end

Y = Y + diag(y(end)*ones(n,1));



end

