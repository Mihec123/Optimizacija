function [ vrednost,X,y,Z, iter, napaka] = ipm_sdp_predcor_maxCut( C,b,X0,y0,sigma,faktor,eps,maxit)
% Opis:
%   ipm_sdp_predcor_maxCut metoda notranjih tock prediktor korektor
%   prilagojena za problem maxcut
%   optimizacijskega problema
%       min 1/4<L,X>
%       p.p. diag(X) = e
%            X >= 0 (X p.s.d.)
%
% Definicija:
%   [ vrednost,x,y,s, iter, napaka] = ipm_sdp_predcor_maxCut( c,b,x0,y0,sigma,natancnost,maxit)
%
% Vhodni  podatki:
%   C minimizacijska funkcija ki mora predstavljati Laplaceovo matriko velikosti m1 x m1 
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
%n = length(b);
[n1,m1] = size(X0);
n = n1;
if nargin < 5
    sigma = 0.5;
end
if nargin < 6
    faktor = 0.8;
end
if nargin < 7
    eps = 1e-6;
end
if nargin < 8
    maxit = 50;
end

X=X0;
y = y0;
Z = C - diag(y);

napaka = Inf;
err_d = Inf;
err_p = Inf;

iter = 0;


while (iter < maxit) && ((napaka > eps) || err_p> eps || err_d>eps )

    iter = iter+1;
    %prvi korak predikotor

    rp = b-diag(X); %tega ne rabmo razn za gledat razliko
    
    %dolocimo A'(y)
%     At = transp_operator(A,y,n1);
    
%     rd1  =C - At-Z;   
    rd = C - diag(y)-Z;
    rc = -Z*X;

    M = zeros(n); %samo diagonala od M razlicna od nic si shranmo kot vektor
    Zinv = inv(Z);
    
    %dolocimo matriko M (mogoce ksn bols nacin to nardit)
    for i = 1:n
        for j=1:n
            M(i,j) = Zinv(i,j)*X(j,i);
        end
    end
    
%     M1 = zeros(n);
%     for i = 1:n
%         for j=1:n
%             M1(i,j) = trace(reshape(A(i,:),[n1,m1])*Zinv*reshape(A(j,:),[n1,m1])*X);
%         end
%     end
    
    
    temp = Zinv*rd*X;
    dy1 = M\(b+diag(temp));
    %dy1 = M\(rp+diag(temp)-diag(Zinv*rc));
    %dolocimo A'(dy1)
    %At = transp_operator(A,dy1,n1);
    
    dZ1 = rd-diag(dy1);
%     dZ11 = rd-At;
    dX1 = Z\(rc-dZ1*X);

    %drugi korak korektor

    tau = (X(:)'*Z(:))/(n1); %to je vprasljivo ce prov
    mu = sigma*tau;
    rcc = mu*eye(n1)-Z*X-dZ1*dX1;
    
    temp1 =Zinv*dZ1*dX1;
    dy = M\(b+diag(temp) - mu*(diag(Zinv)) + diag(temp1));

%     dyy = M1\(b+A*temp(:) - mu*(A*Zinv(:)) + A*temp1(:));
    %dolocimo A'(dy)
    %At = transp_operator(A,dy,n1);
    dZ = rd - diag(dy);
%     dZZ = rd - At;
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
    
    napaka = max(X(:)'*Z(:));
    err_d = norm(diag(y) + Z- C);
    err_p = norm(diag(X)-b);
    
     if mod(iter,2)
        vrednost = C(:)'*X(:);
        text = strcat('iter: ',int2str(iter));
        text = strcat(text,'  err: ');
        text = strcat(text,num2str(napaka));
        text = strcat(text,'  vrednost: ');
        text = strcat(text,num2str(vrednost));
        text = strcat(text,'  err_p: ');
        text = strcat(text,num2str(err_p));
        text = strcat(text,'  err_d: ');
        text = strcat(text,num2str(err_d));
        disp(text)
    end
    

end

vrednost = C(:)'*X(:);

end


function At = transp_operator(A,vektor,velikost)

[v1,~] = size(A);
At = zeros(velikost);
for i= 1:v1
    At = At + vektor(i).*reshape(A(i,:),[velikost,velikost]);
end


end


