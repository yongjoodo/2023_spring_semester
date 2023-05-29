function [x, residual, iter] = minresb(A,b,x0,m,tol,maxiter)
%MINRESB Solve the symmetric system Ax = b using the MINRES method.
%
%   [x r iter] = minresb(A,b,x0,m,tol,maxiter).
%   x0 is the initial approximation to the solution of the sparse system,
%   and integer m < n is the dimension of the Krylov subspace used.
%   tol is the error tolerance, and maxiter is the maximum number of
%   iterations to perform. r is the final residual, and iter is the number
%   of iterations required. iter = -1 if the error tolerance was not met.
%   If this method fails, try gmres or the normal equations with
%   either cg or precg.

iter = 1;
while iter <= maxiter
   r = b-A*x0;
   % do reorthogonalization
   [Q, T] = lanczos(A,r,m,1);
   beta = norm(r);
   e1 = zeros(m+1,1);
   e1(1) = 1;
   y = hesslqsolve(T,beta*e1);
   x = x0 + Q(:,1:m)*y;
   residual = norm(b - A*x);
   if residual < tol
      return;
   end
   x0 = x;
   iter = iter + 1;
end

iter = -1;
