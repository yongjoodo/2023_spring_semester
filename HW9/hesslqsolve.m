function x = hesslqsolve(A, b)
%HESSLQSOLVE solves the least-squares problem Ax = b using the
%QR factorization.
%
%   x = hesslqsolve(A,b) returns the unique least squares solution to
%   Ax = b, where the m x n matrix A has full rank, and m >= n.

[Q, R] = givenshess(A);
c = Q'*b;
x = backsolve(R,c);
