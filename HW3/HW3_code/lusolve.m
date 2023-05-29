function X = lusolve(L, U, P, B)
%lusolve Solve multiple equations Ax = b using
%the result of the LU factorization.
%
%   X = lusolve(L,U,P,B), where B is an n x k matrix containing
%   k right-hand sides for which a solution to the linear system
%   Ax = b is required. L, U, P are the result of the LU factorization
%   P*A = L*U. The solutions are in the k columns of X.

% find the dimension of L, U, P, and the number of
% columns of B
[n, k] = size(B);
%preallocate an n x k matrix to contain the solutions
X = zeros(n,k);   
%each column of pb will be a right-hand side during forward
%substitution
pb = P*B;
%solve k equations
for i=1:k  
	%forward substitution using column i of pb to get y
	y = forsolve(L, pb(:,i));

	%back substitution using y to get solution vector
	%X(:,i)
	X(:,i) = backsolve(U, y);
end
   