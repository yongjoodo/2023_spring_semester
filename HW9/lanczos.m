function [Q, T] = lanczos(A,v0,m,reorthog)
% LANCZOS  Lanczos decomposition of the symmetric matrix A
%
%  Input:   A -- an n by n symmetric matrix
%           v0 -- an initial approximation for beginning the
%                 Lanczos steps
%           m -- dimension of the Krylov subspace, m << n.
%           reorthog -- 1 if reorthogonalization and 0 if not
%
%   Output: Q -- an n by (m+1) orthogonal matrix
%           T -- an (m+1) x m matrix.
%                T(1:m,1:m) is a symmetric tridiagonal matrix
%
%           AQ(:,1:m) = QT

[p n] = size(A);
if p ~= n
   error('The matrix is not square.');
end

if m > n
   error('m must be <= n');
end

if nargin == 3
   reorthog = 0;
end

Q = zeros(n,m+1);
alpha = zeros(m,1);
beta = zeros(m,1);
q = v0/norm(v0);
Q(:,1) = q;

for k = 1:m 
   w = A*q;
   alpha(k) = q'*w;
   if k == 1
      w = w - alpha(1)*q;
   else
      w = w - beta(k-1)*Q(:,k-1) - alpha(k)*q;
   end
	if reorthog == 1
      % full reorthogonalization. slow but effective.
      for i = 1:k-1
         h = Q(:,i)'*w;
         w = w - Q(:,i)*h;
      end
   end
   beta(k) = norm(w);
   if beta(k) < eps
      return;
   end
   q = w/beta(k);
   Q(:,k+1) = q;
end

% construct the (m+1) x m matrix. T(1:m,1:m) is symmetric tridiagonal
T = zeros(m+1,m);
T(1:m,1:m) = diag(beta(1:m-1),-1) + diag(alpha) + diag(beta(1:m-1),1);
T(m+1,m) = beta(m);
