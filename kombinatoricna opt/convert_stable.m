function [ A,b,C ] = convert_stable( e )
% Opis:
%   metoda convert_stable pripravi matrike za re?evanje optimizacijskega
%   problema stabilne mno?ice
%   max <J,X>
%     p.p tr(X) = 1
%         X_ij = 0 za ij in E
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

