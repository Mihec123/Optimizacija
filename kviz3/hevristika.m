function [ S,vrednost,Xk ] = hevristika( e )
% Opis:
%   hevristika izvaja hevristiko iz zadnji_kviz.pdf
% Vhodni  podatki:
%   e matrika m+1 x 3 kjer e(1,1)=n pomeni stevilo vozlisc, e(1,2) = m pomeni
%       stevilo povezav, [e(i,1) e(i,2) e(i,3)] pomeni vozlisce
%       e(i,1) je povezano z vozliscem e(i,2) z utezjo e(i,3), i > 1
%
% Izhodni  podateki:
%   S matrika ? x 1, vozlisc ki so v stabilni mnozici
%   vrednost resiteve z metodo Metoda_robnih_tock_lovasz
%   Xk reitev z metodo Metoda_robnih_tock_lovasz

[ A,b,C ] = convert_stable( e );
sigma = 0.05/e(1);
[ Xk,~,vrednost] = Metoda_robnih_tock_lovasz( A,C,b,sigma,5000,1e-4);
U = 1:e(1);
S = [];

while ~isempty(U)
   %poiscemo X_ii = max{X_jj,j in U}
%    xmax = -Inf;
%    indmax = 0;
%    for i = 1:length(U)
%        xtemp = Xk(i,i);
%        if xmax < xtemp
%            xmax = xtemp;
%            indmax = i;
%        end
%    end
   
   diagonala = diag(Xk);
   diagonala = diagonala(U);
   [xmax,indmax] = max(diagonala);
   indmax = U(indmax);
   S = [S, indmax];
   %odstranimo indmax in vse njegove sosede iz U
   odstrani = [];
   for i = 1:length(A);
       if A(i,1)==indmax || A(i,2)==indmax
           %je povezave med indmax in nekim vozliscem
           if A(i,1) == indmax
               odstrani = [odstrani,A(i,2)];
           else
               odstrani = [odstrani,A(i,1)];
           end
       end
   end
   
   odstrani = [odstrani,indmax];
   
   U = setdiff(U,odstrani);
    
    
end


end

