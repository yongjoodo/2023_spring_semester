function [x, residual, iter] = cg(A,b,x0,tol,maxiter)
%cg  Conjugate gradient method
%
%   [x,residual,iter] = cg(A,b,x0,tol,numitr) computes the solution x 
%   of a sparse symmetric positive definite system Ax = b using
%   the Conjugate Gradient method. x0 is the initial approximation, tol is
%   the error tolerance, and maxiter is the maximum number of iterations
%   to execute. If the conjugate gradient method converges, iter contains
%   the number of iterations required to converge; otherwise, iter = -1.
%   For consistent failure of convergence, try precg.

r = b - A*x0;
p = r;
x = x0;
tol = tol^2;

for i = 1:maxiter
   normrsqr = r'*r;
   w = A*p;
   alpha = normrsqr/(p'*w);
   x = x + alpha*p;
   r = r - alpha*w;
   normrnewsqr = r'*r;
   if normrnewsqr < tol
      iter = i;
      residual = sqrt(normrnewsqr);
      return;
   end
   beta = normrnewsqr/normrsqr;
   p = r + beta*p;
end

iter = -1;
residual = sqrt(normrnewsqr);
