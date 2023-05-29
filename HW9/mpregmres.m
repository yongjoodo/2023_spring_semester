function [x, r, iter] = mpregmres(A,b,x0,m,tol,maxiter,droptol,reorthog)
%mpregmres Solve the n x n linear system Ax = b using incomplete
%LU factorization as a preconditioner for GMRES.
%
%   [x r iter] = mpregmres(A,b,x0,m,tol,maxiter) preconditions the
%   sparse system Ax = b using incomplete LU factorization and attempts
%   to find a solution to Ax = b. x0 is the initial approximation,
%   m < n is the dimension of the Krylov subspace used, tol is
%   the requested error tolerance, and maxiter is the maximum
%   number of iterations to perform. r is the residual obtained,and
%   iter = -1 if the iteration fails.

[p, q] = size(A);
if p ~= q
   error('Matrix must be square.');
end

setup.type = 'ilutp';
if nargin == 6
   reorthog = 0;
   droptol = 1.0e-6;
elseif nargin == 7
   reorthog = 0;
end
setup.droptol = droptol;
% setup.milu = 'row';
setup.udiag = 1;
[L, U, P] = ilu(A,setup);

if nargin == 6
   reorthog = 0;
end

b_bar = (P'*L)\b;
B = (P'*L)\(A/U);
[x_bar , ~, iter] = gmresb(B,b_bar,x0,m,tol,maxiter,reorthog);
x = U\x_bar;
r = norm(A*x-b);
