function [psi, x, y, iter, secs] = ipm_lp_longstep_hrga(c,A,b,x0,y0)
% IPM_LP_LONGSTEP (long-step primal-dual interior point method) for LP solves:
% 
% min c'x, subject to Ax = b, x >= 0
%
% and its dual
%
% max b'y, subject to A'y + s = c, s >= 0
%
% input:  c     ... cost vector
%         A, b  ... data determening feasible reagion
%         x0,y0 ... starting striclty feasible vectors, x0 > 0 and s0 =
%         c-A'y0 > 0
%
% output: psi  ... optimal value of LP
%         x    ... optimal primal vector
%         y    ... optimal dual vector
%         iter ... number of iterations
%         secs ... running time in seconds
%
% call:   [psi, x, y, iter, secs] = ipm_lp_longstep(c,A,b,x0,y0)
%
% Note: This is long-step method:
% - mu is reduced more aggressively
% - update step: damped Newton step

tic;

digits = 6;                   % significant digits of psi
n = length(c);                % size of problem


% initial strictly feasible vectors
x = x0;
y = y0;
s = c - A'*y0;

mu = s'*x/n;                  % initial complementarity

iter = 0;                     %iteration count

phi = b'*y;                   % objective values
psi = c'*x;

% while duality gap too large
while (psi - phi) > max(1,abs(phi))*10^(-digits) 
    
    iter = iter + 1;                            % start new iteration
    mu = 0.5*mu;                                % decrease mu
    
    X = diag(x);
    Si = diag(s.^(-1));                         % explicitly compute inv(S)
    ASi = A*Si;
    
    dy = ( ASi*X*A' ) \ ( b-mu*ASi*ones(n,1) ); % solve for dy
    ds = - A'*dy;                               % back substitute for ds 
    dx = mu*Si*ones(n,1) - x - Si*X*ds;         % back substitute for dx 

     % line search on primal: x  = x + alpha_p * dx  positive
    alpha_p = 1;
    positive = all((x + alpha_p*dx)>0);  % test if positive definite
    
    while positive == 0
       alpha_p = 0.8 * alpha_p;
       positive = all((x + alpha_p*dx)>0);
    end
    
 
    % line search on dual
    alpha_d = 1;
    positive = all((s + alpha_d*ds)>0);
    
    while positive == 0
       alpha_d = 0.8 * alpha_d;
       positive = all((s + alpha_d*ds)>0);
    end
    
    % update
    x = x + alpha_p * dx;
    y = y + alpha_d * dy;
    s = s + alpha_d * ds;
    mu = s'*x/n;
    
    % objective values
    phi = b'*y;
    psi = c'*x;
    
end

secs = toc;


end


