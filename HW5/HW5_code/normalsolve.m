function [x, residual] = normalsolve(A, b)
%normalsolve Solve the least-squares problem Ax = b using the
%normal equations.
%
%   x = normalsolve(A,b) solves the least-squares problem Ax = b,
%   where A is an m x n full rank matrix, m >= n.
%   This not normally a satisfactory way to solve a least-squares
%   problem.

R = chol(A'*A);
x = cholsolve(R, A'*b);
residual = norm(b-A*x);
