function X = forsolve(L,B)
% forsolve Perform forward substitution.
%
%   x = forsolve(L,B) uses forward substitution to
%   compute the solution to LX = B, where U is a square
%   lower triangular matrix, and B is a matrix of right-hand sides.

[m, n] = size(L);
if m ~= n
   disp('The matrix is not square');
   return;
end

%preallocate an n x k matrix to contain the solutions
[n, k] = size(B);
X = zeros(n,k);
for i = 1:k
%  Solve Lx = b by forward substitution
	x = X(:,i);
	b = B(:,i);
	x(1) = b(1)/L(1,1);
	for j = 2:n
	   x(j) = (b(j) - L(j,1:j-1)*x(1:j-1))/L(j,j);
	end
	X(:,i) = x;
end
