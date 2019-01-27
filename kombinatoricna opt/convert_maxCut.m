function [ A,b,C ] = convert_maxCut( e )
% Opis:
%   metoda convert_stable pripravi matrike za re?evanje optimizacijskega
%   problema max cut
%   max <1/4L,X>
%     p.p diag(X) = e
%         X s.d.p
% Vhodni  podatki:
%   e matrika m+1 x 3 kjer e(1,1)=n pomeni stevilo vozlisc, e(1,2) = m pomeni
%       stevilo povezav, [e(i,1) e(i,2) e(i,3)] pomeni vozlisce
%       e(i,1) je povezano z vozliscem e(i,2) z utezjo e(i,3), i > 1
%
% Izhodni  podateki:
%   A matrika m x n^2, ki doloca pogoje v sdp programu
%   b matrika m x 1 ki doloca desno stran pogojev v sdp programu
%   C matrika n x n, katero optimiziramo pri sdp programu

m = e(1,2);
n = e(1,1);

A = zeros(n,n^2);
b = ones(n,1);

Atemp = zeros(n); %matrika sosednosti
for i = 1:m
    Atemp(e(i+1,1),e(i+1,2)) = e(i+1,3);
    Atemp(e(i+1,2),e(i+1,1)) = e(i+1,3);
end
Atemp

C = 1/4*(diag(Atemp*ones(n,1))-Atemp);

for i = 1:n
    Ai = zeros(n);
    Ai(i,i) = 1;
    A(i,:) = Ai(:);
end


end

