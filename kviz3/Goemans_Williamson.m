function [ opt] = Goemans_Williamson( e,X0,y0,iter )
% Opis:
%   metoda Goemans_Williamson izracuna optimalno vrednost lovaszeve theta
%   funkcije na grafu e, ter poisce priblizno vrednost najboljsega prereza
% Vhodni  podatki:
%   e matrika m+1 x 3 kjer e(1,1)=n pomeni stevilo vozlisc, e(1,2) = m pomeni
%       stevilo povezav, [e(i,1) e(i,2) e(i,3)] pomeni vozlisce
%       e(i,1) je povezano z vozliscem e(i,2) z utezjo e(i,3), i > 1
%   X0 zacetna matrika za ipm_sdp_predcor_maxCut
%   y0 zacetni vektor za ipm_sdp_predcor_maxCut
%
% Izhodni  podateki:
%   opt priblizno vrednost najboljsega prereza

if nargin < 4
    iter = 150;
end

[ b,C ] = convert_maxCut( e );

[vrednost,X] = ipm_sdp_predcor_maxCut( -C,b,X0,y0);
disp(vrednost);
[S,D] = eig(X);

V = sqrtm(D)*S';
opt = -Inf;

for j = 1:iter
    r = rand(length(V),1);
    %r = ones(length(V),1)/norm(ones(length(V),1));


    Vtemp = r'*V;
    V1 = Vtemp > 0;
    %V2 = Vtemp < 0; ne rabmo

    %vrednost maxcuta
    maxcut = 0;
    for i = 1:e(1,2)
        %pogledamo ce povezava v razlicni mnozici
        %prva tocka
        T1 = e(i+1,1);
        T2 = e(i+1,2);
        if (V1(T1)==0 && V1(T2) == 1) || (V1(T1)==1 && V1(T2) == 0)
            maxcut = maxcut + e(i+1,3);
        end
    end
    if opt < maxcut
        opt = maxcut;
    end
end
        


end

