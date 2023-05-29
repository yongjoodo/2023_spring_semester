function [x, r, iter] = gmresb(A,b,x0,m,tol,maxiter,reorthog)
%gmresb Solve Ax = b using the GMRES method.
%
%   [x r iter] = gmresb(A,b,x0,m,tol,maxiter,reorthog).
%   x0 is the initial approximation to the solution of the sparse system
%   Ax = b, and integer m < n is the dimension of the Krylov subspace used.
%   tol is the error tolerance, maxiter is the maximum number of
%   iterations to perform, and if reorthog = 1, reorthogonalization
%   is performed. r is the final residual norm, and iter is the number
%   of iterations required. iter = -1 if the error tolerance was not met.
%   If this method fails, try pregmres or the normal equations using
%   either cg or precg.

if nargin == 6
   reorthog = 0;
end

iter = 1;
mm = m;
while iter <= maxiter
   r = b-A*x0;
   [Q, H, m] = arnoldi(A,r,m,reorthog);
   beta = norm(r);
   e1 = zeros(m+1,1);
   e1(1) = 1;
   y = hesslqsolve(H,beta*e1);
   x = x0 + Q(:,1:m)*y;
   r = norm(b - A*x);
   if r < tol
      return;
   end
   x0 = x;
   iter = iter + 1;
   m = mm;
end

iter = -1;
