function A = pentd(a,b,c,d,e)
%pentd Construct an n x n pentadiagonal matrix
%
%   A = pentd(a,b,c,d,e) builds a pentadiagonal matrix, where a [(n-2) x 1]
%   and b [(n-1) x 1] are the subdiagonals, c [n x 1] is the primary
%   diagonal, and d [(n-1) x 1], e [(n-2) x 1] are the superdiagonals.

A = diag(a,-2) + diag(b,-1) + diag(c) + diag(d,1) + diag(e,2);
% A = spdiags([a b c d e], -2:2, n, n);

% A = diag(a,-1) + diag(b) + diag(c,1);