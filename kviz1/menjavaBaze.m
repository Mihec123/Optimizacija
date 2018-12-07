function J = menjavaBaze(c,A,b,J,izstopna,opcija)
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
% opcija 
%       1 pomeni simpleksna metoda za vhodno spremenljivko izbere tisto z
%         najmanjso vrednostjo
%       2 pomeni simpleksna metoda za vhodno spremenljivko izbere tisto z
%         najmanjsim indeksom
%Izhod:
% J baze po eni menjavi

if nargin < 6
    opcija = 2;
end

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

if opcija == 1
    temp(temp == 0) = Inf;
    [neki,indeks] = min(temp); %vzamemo najmanjo vrednost nebazne spremenljivke, ki nima koeficienta 0 ter ga damo v bazo
    ro = K(indeks);
else
        %za vstopno spremenljivko vzamemo tisto z najmanjsim indeksom
    temp1 = temp>0; %temp1 vsebuje na i-tem mestu 1 ce je cpr(i) >0 in 0 sicer
    seznam = find(temp1 == 1); % seznam vseh mest enakih 1 (cpr(i) <0)
    indeks = seznam(1); %vzamemo prvega (najmanjsi indeks)
    ro = K(indeks);
    
end

apr = A(:,J)\A(:,ro);

ro1 = izstopna; % za izstopno spremenljivko izberemo iztstopna

J(J == ro1) = [];
J = sort([J ro]);

xj1 = A(:,J)\b;
vr1 = c(J)'*xj;

%preverimo ce imamo se vedno dopustno resitev, ter da je vrednost
%kriterijske funkcije enaka kot pred menjavo baze
if vr1 ~= vr || prod(xj1 < 0)
    display('menjava baze ni uspela')
end

end

