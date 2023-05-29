
function [V,T,f] = lanczosf(A,V,T,f,k,v,m,reorthog)
%LANCZOSF computes the Lanczos decomposition of matrix A
%
%  Input:   A -- an n by n symmetric matrix
%           V -- an n x m matrix of eigenvectors, zero before first call
%           T -- an m x m tridiagonal matrix, zero before first call
%           f -- a column vector of dimension n
%           k -- start of the Lanczos steps
%           v -- on first call, a random vector or a known approximation
%                to an eigenvector of A
%           m -- dimension of the Krylov subspace, m << n.
%           reorthog -- 1 if full reorthogonalization and 0 if not
%
%   Output: V -- an n by k orthogonal matrix
%           T -- an m x m symmetric tridiagonal matrix
%           f -- an n vector
%
%           AV = VT + f*e_m', where e_m is the standard basis
%           vector [0 0 ...0 1]'

%   lanczos can be coded more efficiently by representing T as two
%   vectors, sub and super diagonal beta and diagonal alpha. speed
%   improvement can be done using selective reorthogonalization

%   Original author: D.C. Sorensen
%   modified by William Ford
%   15 May 2013

n = length(v);
v1 = v/norm(v);
if nargin == 7
   reorthog = 1;
end

if k == 1
   f = A*v1;
   alpha = v1'*f;
   f = f - v1*alpha;

   V(:,1) = v1;
   T(1,1) = alpha;
   k = 2;
end

for j = k:m,
   beta = norm(f); 
   if beta < eps
      beta = eps;
   end
   v0 = v1;
   v1 = f/beta;

   % full reorthogonalization. slow but effective.
   if reorthog == 1
      for i = 1:j-1
         h = V(:,i)'*v1;
         v1 = v1 - V(:,i)*h;
      end
   end

   f = A*v1 - v0*beta;

   alpha = v1'*f;
   f = f - v1*alpha;

   T(j,j-1) = beta;
   T(j-1,j) = beta;
   T(j,j) = alpha;
   V(:,j) = v1;
end
