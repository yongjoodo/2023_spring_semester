function X = backsolve(U,B)
% backsolve Perform back substitution.
%
%   x = backsolve(U,b) finds the solution to UX = B,
%   where U is a square upper triangular matrix, and
%	 B is a matrix of right-hand sides.

[m, n] = size(U);
if m ~= n
   error('The matrix is not square');
end

[n, k] = size(B);
%preallocate an n x k matrix to contain the solutions
X = zeros(n,k);
for i = 1:k
%  Solve Ux = b by back substitution
	x = X(:,i);
	b = B(:,i);
	x(n) = b(n)/U(n,n);
	for j = n-1:-1:1
	   x(j) = (b(j) - U(j, j+1:n)*x(j+1:n))/U(j,j);
% 		sum = 0.0;
% 		for j = i+1:n
% 			sum = sum + U(i,j)*x(j);
% 		end
% 		x(i) = (b(i)-sum)/U(i,i);
	end
	X(:,i) = x;
end
