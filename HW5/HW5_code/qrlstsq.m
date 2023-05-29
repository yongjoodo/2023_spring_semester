function [x, residual] = qrlstsq(A, b)
% qrsolve Solve the least-squares problem Ax = b using the
% QR factorization.
%
%   x = qrsolve(A,b) solves the m x n least squares problem, m >= n,
%   using the thin QR factorization of A.

[Q, R] = qr(A,0);
c = Q'*b;
x = backsolve(R,c);
residual = norm(b-A*x);
