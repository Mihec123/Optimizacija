% test QP solver

rng(1);

n = 2000;               % dimension of x
m = 1000;                % number of equality constraints
    
x0 = randn(n,1);        % create random solution vector

A = randn(m,n);         % create random  matrix A
b = A*x0;

% generate a well-conditioned positive definite matrix

P = rand(n);
P = P + P';
[V, ~] = eig(P);
P = V*diag(1+rand(n,1))*V';
P = 0.5*(P+P');

q = randn(n,1);


[ x,vrednost,err_p,err_d,iter] = okrepljena_lagrangeeva_metoda_kvadraticni( P,q,A,b,0.0005,1000,1e-4);