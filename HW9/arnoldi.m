function [Q, H, m] = arnoldi(A,v0,m,reorthog)
% ARNOLDI  computes the Arnoldi decomposition of matrix A
%
%
%  Input:   A -- an n by n matrix
%           v0 -- an initial approximation for beginning the
%                 Arnoldi steps
%           m --  dimension of the Krylov subspace, m << n.
%           reorthog -- 1 if full reorthogonalization and 0 if not
%
%   Output: Q -- an n by (m+1) orthogonal matrix
%           H -- an (m+1) x m upper Hessenberg matrix
%           m -- if H(k,k+1) is very small, we have found
%                 an invariant subspace. return the new m.
%
%           AQ(:,1:m) = Q*H

[p, n] = size(A);
if p ~= n
   error('The matrix is not square.');
end

if m > n
   error('m must be <= n');
end

if nargin == 3
	reorthog = 0;
end

H = zeros(m+1,m);
Q = zeros(n,m+1);
tol = n*eps;

Q(:,1) = v0/norm(v0);
for k = 1:m
   w = A*Q(:,k);
	ow = norm(w);
   for j = 1:k
      H(j,k) = (Q(:,j))'*w;
      w = w - H(j,k)*Q(:,j);
	end
	if reorthog == 1
      for j=1:k
         tmp = (Q(:,j))'*w;
         w = w - tmp*Q(:,j);
         H(j,k)=H(j,k)+tmp;
      end
   end

	H(k+1,k) = norm(w);
   if H(k+1,k) <= tol*ow
      m = k;
      H = H(1:m+1,1:m);
      Q = Q(1:n,1:m+1);
      return;
	end
   Q(:,k+1) = w/H(k+1,k);
   
end
