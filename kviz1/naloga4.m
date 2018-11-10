%naloga4
%graf incidencna matrika wikipedija

G = [1 0 0 0 1 0; 1 1 0 0 0 1; 0 1 1 0 0 0; 0 0 1 1 0 0; 0 0 0 1 1 1];
%zmisljene
utezi = [0.5, 0.8, 0.85, 0.8, 0.7, 0.65];

%[X,vr] = popolnoPrirejanjeCeloSt(utezi,G)

[X1,vr1] = popolnoPrirejanje(utezi,G)

%dvodeln graf


G1 = [1 1 0; 1 0 0; 0 0 1; 0 1 1];
utezi1 = [0.85,0.9,0.7];

%[X2,vr2] = popolnoPrirejanjeCeloSt(utezi1,G1)

[X3,vr3] = popolnoPrirejanje(utezi1,G1)
