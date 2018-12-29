function [ vrednost,X,y,Z, iter, napaka] = ipm_sdp_predcor( c,A,b,X0,y0,sigma,faktor,eps,maxit)
% Opis:
%   ipm_sdp_predcor metoda notranjih tock prediktor korektor
%   optimizacijskega problema
%       min <C,X>
%       p.p. A(X) = b
%            X >= 0 (X p.s.d.)
%
%       kjer je A(X) = [<A1,X>; ...;<An,X>] in A'(y) = sum (y_i.*Ai)
% Definicija:
%   [ vrednost,x,y,s, iter, napaka] = ipm_sdp_predcor( c,A,b,x0,y0,sigma,natancnost,maxit)
%
% Vhodni  podatki:
%   c minimizacijska funkcija (min c'*X(:)) predstavljena z matriko
%       velikosti m x 1, kjer je c = C(:) in C neka m1 x m1 matrika torej
%       za m velja m = m1*m1
%   A matrika velikosti n x m kjer vsaka vrstica predstavlja matriko Ai(:)'
%       in velja da je matrika Ai velikosti m1 x m1
%   b desna stran sistema A(X) = b predstavljena z matriko velikosti n x 1
%   X0 zacetni priblizek za resitev opt. problema velikosti m1 x m1 veljati
%       mora da je X0 strogo pozitivno definitna
%   y0 zacetni priblizek za resitev dualnega opt. problema velikosti n x 1
%       veljati mora da je y0 pozitiven
%   sigma vrednost, ki doloca premik proti srediscni poti default vrednost sigma
%       je sigma = 0.5
%   faktor je vrednost med 0<faktor<1, ki pove s kaksnim faktorjem manjsamo iskanje
%       alfe pri newtnovi metodi v smeri X in Z. Default value = 0.8
%   eps je vrednost, ki nam pove toleranco napake dobljene resitve
%       default vrednost natancnosti je eps = 1e-6
%   maxit pove maksimalno stevilo iteracij, ki jih izvedemo
%       default vrednost maxit je maxit = 100
%
% Izhodni  podateki:
%   vrednost nam pove vrednost optimizacijskega problema min <C,X>
%   X matrika m1 x m1 vrednosti X optimalnega problema 
%   y matrika n x 1 vrednosti y dualnega problema
%   Z matrika m1 x m1 vrednosti dopolnilnih spremenljivk dualnega problema
%   iter stevilo iteracij, ki jih je metoda izvedla
%   napaka velikost napake metode merjene kot X(:)'*Z(:)


%velikosti v funkciji se ne ujemajo z oznakami iz opisa
[n,~] = size(A);
[n1,m1] = size(X0);

if nargin < 6
    sigma = 0.5;
end
if nargin < 7
    faktor = 0.8;
end
if nargin < 8
    eps = 1e-6;
end
if nargin < 9
    maxit = 100;
end

X=X0;
y = y0;
At = transp_operator(A,y,n1);
C = reshape(c,[n1,m1]);
Z = C - At;

napaka = X(:)'*Z(:);

iter = 0;


while (napaka > eps) && (iter < maxit)

    iter = iter+1;
    %prvi korak predikotor

    rp = b-A*X(:); %tega ne rabmo razn za gledat razliko
    
    %dolocimo A'(y)
    At = transp_operator(A,y,n1);
    
    rd = C - At-Z;
    rc = -Z*X;
    
    
    M = zeros(n);
    Zinv = inv(Z);
    
    %dolocimo matriko M (mogoce ksn bols nacin to nardit)
    for i = 1:n
        for j=1:n
            M(i,j) = trace(reshape(A(i,:),[n1,m1])*Zinv*reshape(A(j,:),[n1,m1])*X);
        end
    end
    
    
    temp = Zinv*rd*X;
    dy1 = M\(b+A*temp(:));
    
    %dolocimo A'(dy1)
    At = transp_operator(A,dy1,n1);
    
    dZ1 = rd-At;
    dX1 = Z\(rc-dZ1*X);

    %drugi korak korektor

    tau = (X(:)'*Z(:))/(n1*m1); %to je vprasljivo ce prov
    mu = sigma*tau;
    rcc = mu*eye(n1)-Z*X-dZ1*dX1;
    
    temp1 =Zinv*dZ1*dX1;
    dy = M\(b+A*temp(:) - mu*(A*Zinv(:)) + A*temp1(:));

    %dolocimo A'(dy)
    At = transp_operator(A,dy,n1);
    dZ = rd - At;
    dX = Z \ (rcc-dZ*X);
    
    dX = (dX+dX')./2;

    %pogledamo kako velik korak lahko naredimo v smeri X da smo
    %se dopustni
    alfap = 1;
    napakax = true;
    while napakax
        %preverjamo psd matrike X tako da probamo nardit razcep choleskega
        try
            V = chol(X+alfap.*dX);
            napakax = false;
        catch
            alfap = faktor*alfap;
        end
    end

    %pogledamo kako velik korak lahko naredimo v smeri Z da smo
    %se dopustni
    alfad = 1;
    napakaz = true;
    while napakaz
        try
            V = chol(Z+alfad.*dZ);
            napakaz = false;
        catch
            alfad = faktor*alfad;
        end
     
    end


    X = X+alfap.*dX;
    y = y + alfad*dy;
    Z = Z+alfad.*dZ;
    
    napaka = X(:)'*Z(:);
    

end

vrednost = c'*X(:);

end



function At = transp_operator(A,vektor,velikost)

[v1,~] = size(A);
At = zeros(velikost);
for i= 1:v1
    At = At + vektor(i).*reshape(A(i,:),[velikost,velikost]);
end


end


