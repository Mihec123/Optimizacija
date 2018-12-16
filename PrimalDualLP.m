function [x,y,s,f,iter] = PrimalDualLP(A,b,c)
% Experimental version of Mehrotra?s primal-dual interior-point method
% for linear programming
%
% primal problem: min c?x s.t. Ax = b, x >= 0
% dual problem: max b?y s.t. A?y + s = c, s >= 0
Maxiter = 100; Tol = 1.e-8;
UpperB = 1.e10*max([norm(A),norm(b),norm(c)]);
[m,n] = size(A);
% starting point
e = ones(n,1); M = A*A'; R = chol(M);
x = A'*(R \ (R' \ b)); y = R \ (R' \ (A*c)); s = c - A'*y;
delta_x = max(-1.5*min(x),0)+1/n; x = x+delta_x*e;
delta_s = max(-1.5*min(s),0)+1/n; s = s+delta_s*e;
for iter = 0 : Maxiter
    
f = c'*x;
% residuals
rp = A*x-b; % primal residual

rd = A'*y+s-c; % dual residual
rc = x.*s; % complementarity
tau = x'*s/n; % duality measure
residue = norm(rc,1)/(1+abs(b'*y));
primalR = norm(rp)/(1+norm(x)); dualR = norm(rd)/(1+norm(y));
STR1 = 'iter %2i: f = %14.5e, residue = %10.2e';
STR2 = 'primalR = %10.2e, dualR = %10.2e\n';
fprintf([STR1 STR2], iter, f, residue, primalR, dualR);
if (norm(x)+norm(s) >= UpperB)
error('Problem possibly infeasible!'); end
if (max([residue; primalR; dualR]) <= Tol) break; end
% coefficient matrix of the linear systems
M = A*diag(x./s)*A'; R = chol(M); % M = R?R
% predictor step with maximal Newton step size
rhs = rp - A*((rc-x.*rd) ./ s);
dy = R \ (R' \ rhs); % dy = M \ rhs
ds = rd-A'*dy; dx = (rc-x.*ds) ./ s;
alpha_p = 1/max([1; dx./x]); alpha_d = 1/max([1; ds./s]);
tau_N = ((x-alpha_p*dx)'*(s-alpha_d*ds))/n;
% corrector step: correct towards center path
sigma = (tau_N/tau)^3; % Mehrotra
rc = rc - sigma*tau + dx.*ds;
rhs = rp - A*((rc-x.*rd) ./ s);
dy = R \ (R' \ rhs); % dy = M \ rhs
ds = rd-A'*dy; dx = (rc-x.*ds) ./ s;
eta = max(0.99,1-tau);
alpha_p = eta/max([eta; dx./x]); alpha_d = eta/max([eta; ds./s]);
x = x-alpha_p*dx; y = y-alpha_d*dy; s = s-alpha_d*ds;
end % for-loop

