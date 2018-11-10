function J = menjavaBaze(c,A,b,J,izstopna)
%funkcija menjavaBaze izvede en korak menjave baze v simpleksni metodi.
%Funkcija je namenjena menjavi nezazeljenih spremenljivk v 1.fazi dvofazne
%simpleksne metode (zamnejamo eno dopolnilno spremenljivko za prvo ustrezno)
%
%Vhod:
% c nx1 matrika kriterijske funkcije
% A mxn matrika
% b mx1 matika
% J mx1 matrika trenutne baze
% izstopna 1x1 matrika z izstopno spremenljivko
%Izhod:
% J baze po eni menjavi

V = size(A);
m = V(1);
n = V(2);
K = [1:n];
K(J) = 0;
K(K==0) = [];
x = zeros(length(c),1);

xj = A(:,J)\b;
vr = c(J)'*xj;
y = A(:,J)'\c(J);
neki = A(:,J)\A(:,K);
cpr = c(K)'-c(J)'*neki;

temp = cpr;
temp(temp == 0) = Inf;
[neki,indeks] = min(temp); %vzamemo najmanji indeks nebazne spremenljivke, ki nima koeficienta 0 ter ga damo v bazo
ro = K(indeks);

apr = A(:,J)\A(:,ro);

ro1 = izstopna; % za izstopno spremenljivko izberemo iztstopna

J(J == ro1) = [];
J = sort([J ro]);

xj1 = A(:,J)\b;
vr1 = c(J)'*xj;

%preverimo ce imamo se vedno dopustno resitev, ter da je vrednost
%kriterijske funkcije enaka kot pred menjavo baze
if prod(vr1 ~= vr) || prod(xj1 < 0)
    display('menjava baze ni uspela')
end

end

