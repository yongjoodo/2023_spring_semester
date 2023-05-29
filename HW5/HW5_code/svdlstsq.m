function [x, residual] = svdlstsq(A,b)
% svdsolve Use the reduced SVD to solve the least-squares
% problem Ax = b.
%
%   x = svdsolve(A,b) uses the SVD to solve the least-squares
%   problem Ax = b, where A is an m x n matrix, m >= n.
%   If A has full rank, the solution is unique.

n = size(A,1);

[U, S, V] = svd(A,0);
c = U'*b;
y = c./diag(S);
x = V*y;
residual = norm(b-A*x);
