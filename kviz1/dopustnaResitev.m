function J = dopustnaResitev(c,A,b,maxstep,opcija)
%poisce dopustno resitev za linearni program ce ta obstaja. Izvedemo 1.fazo
%dvofazne simpleksne metode.
%Vhod:
% c nx1 matrika kriterijske funkcije
% A mxn matrika
% b mx1 matika
% maxstep najvecje stevilo korakov simpleksne metode, ce ga ne podamo
%   izvedemo najvec 100 korakov
% opcija 
%       1 pomeni simpleksna metoda za vhodno spremenljivko izbere tisto z
%         najmanjso vrednostjo
%       2 pomeni simpleksna metoda za vhodno spremenljivko izbere tisto z
%         najmanjsim indeksom
%Izhod:
% baza za simpleksno metodo originalnega problema

if nargin < 5
    opcija = 2;
end
if nargin < 4
    maxstep = 100;
end

Dim = size(A);
m = Dim(1);
n = Dim(2);

A1 = [A eye(m)];
J1 = [n+1:1:n+m];
c1 = [zeros(1,length(c)) ones(1,m)]';
[x,vr,y,st,J] = simpleksMetoda( c1,A1,b,J1,maxstep,opcija);
if vr ~=0
    display('ni dopustne resitve')
else
    display('baza pred menjavo')
    display(J)
    while ~isempty(intersect(J,J1))
        temp = intersect(J,J1); %preverimo ce je v trenutni bazi kaksna dopolnilna spremenljivka
        izstopna = temp(1);
        J = menjavaBaze(c1,A1,b,J,izstopna,opcija); %zamenjamo dopolnilno spremenljivko s prvo ki nam ne pridela singularne matrike
    end
    display('baza po menjavi')
    display(J)
    
end
end

