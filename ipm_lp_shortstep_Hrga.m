function [psi,x, y, iter, secs] = ipm_lp_shortstep_Hrga(c,A,b,x0,y0)
% IPM_LP_SHORTSTEP (short-step primal-dual interior point method) for LP solves:
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
% call:   [psi, x, y, iter, secs] = ipm_lp_shortstep(c,A,b,x0,y0)
%
% Note: This is short-step method:
% - mu is reduced more slowly
% - update step: full Newton step

tic;

digits = 6;                   % significant digits of psi
n = length(c);                % size of problem

% initial strictly feasible vectors
x = x0;
y = y0;
s = c - A'*y0;

tau = 1-0.4/(sqrt(n));        % centering parameter

mu = s'*x/n;                  % initial complementarity

iter = 0;                     % iteration count

phi = b'*y;                   % objective values
psi = c'*x;

% while duality gap too large
while (psi - phi) > 10^(-digits) 
    
    iter = iter + 1;                            % start new iteration
    mu = tau*mu;                                % decrease mu
    
    X = diag(x);
    Si = diag(s.^(-1));                         % explicitly compute inv(S)
    ASi = A*Si;
    
    dy = ( ASi*X*A' ) \ ( b-mu*ASi*ones(n,1) ); % solve for dy
    ds = - A'*dy;                               % back substitute for ds 
    dx = mu*Si*ones(n,1) - x - Si*X*ds;         % back substitute for dx 
    
    % update: full Newton step
    x = x + dx;
    y = y + dy;
    s = s + ds;
    
    mu = s'*x/n;
    
    % objective values
    phi = b'*y;
    psi = c'*x;
    
end

secs = toc;


end