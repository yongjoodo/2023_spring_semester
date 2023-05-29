function [x, residual, iter] = precg(A,b,x0,tol,maxiter,method,droptol)
%precg  Preconditioned conjugate gradient method using either
%incomplete Cholesky or SSOR
%
%   [x,iter] = precg(A,x0,b,tol,numitr) computes the solution x 
%   of a symmetric positive definite system Ax = b using
%   the preditioned Conjugate Gradient method. The preconditioning is
%   performed using the incomplete Cholesky factorization or
%   SSOR.
%Input: A - the sparse coefficient matrix.
%       b - right-hand side of the system Ax = b.
%       x0 - the initial approximation to x.
%       tol - the error tolerance required
%       maxiter - the maximum number of iterations to execute.
%       method - the string 'incomplete Cholesky' or 'SSOR'.
%       droptol - the drop tolerance if the preconditioning
%                 method is incomplete Cholesky.
%If method is omitted, incomplete Cholesky is assumed with zero-fill. If
%method is specified and droptol is not, droptol defaults to 1.0e-4.
%Output:
%       x - approximate solution.
%       residual - the value of norm(b - A*x).
%       iter = the number of iterations required if the method attains
%       the error tolerance and -1 otherwise.

if nargin < 5
    error('precg requires at least 5 arguments.');
end

% prepare for possible SSOR preconditioner
L = tril(A,-1);
D = diag(diag(A));
U = triu(A,1);

incompleteCholesky = false;
r = b - A*x0;

if nargin == 5
    method = 'incomplete Cholesky';
    opt.shape = 'upper';
    opt.type = 'nofill';
    R = ichol(A,opt);
    z = R\(R'\r);
	 incompleteCholesky = true;
elseif strcmp(method,'incomplete Cholesky')
    if nargin == 7
        opt.type = 'ict';
        opt.droptol = droptol;
    else
        opt.type = 'ict';
        opt.droptol = 1.0e-4;
	 end
	 incompleteCholesky = true;
    opt.shape = 'upper';
    R = ichol(A,opt);
    z = R\(R'\r);
elseif strcmp(method,'SSOR');
    z = (U+D)\(D*((L+D)\r));
else
	 error('The method must be ''incomplete Cholesky'' or ''SSOR''');
end
p = z;
x = x0;

for i = 1:maxiter
   alpha = (r'*z)/(p'*(A*p));
   x = x + alpha*p;
   prevr = r;
   prevz = z;
   r = r - alpha*(A*p);
   residual = norm(r);
   if residual < tol
      iter = i;
      return;
   end
   if incompleteCholesky
       z = R\(R'\r);
   else
       z = (U+D)\(D*((L+D)\r));
   end
   beta = (r'*z)/(prevr'*prevz);
   p = z + beta*p;     
end

iter = -1;
