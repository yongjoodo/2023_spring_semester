function [Q,R] = clqrgrsch(A)
%clgrsch Classical Gram-Schmidt for QR factorization.
%
%   [Q,R] = clgrsch(A) computes the reduced QR factorization
%   of an m x n  matrix A using the classical Gram-Schmidt method.
%   A = QR, where Q is m x n and has orthonormal columns, and
%   R is n x n upper triangular matrix.
 
[m, n] = size(A);
% allocate the output matrices Q and R
Q = zeros(m, n);
R = zeros(n, n);

% begin the computation
for i=1:n
	%compute entries of the upper triangular matrix in column i.
	%we must form the inner products of the vectors in columns
	%1 to i-1 of Q and the m elements in column i of A
	R(1:i-1,i) = Q(:,1:i-1)'* A(:,i);
	sum_proj = zeros(m,1);
	% sum the projections of Q(:,j) onto A(:,i), 1 <= j <= i-1
	% with the sum assigned to sum_proj
	for j=1:i-1
		sum_proj = sum_proj + Q(:,j)'*A(:,i)*Q(:,j);
	end
	% remove the projections from A(:,i)
	Q(:,i) = A(:,i) - sum_proj;
	% complete column i of R by computing the diagonal element
	R(i,i) = norm(Q(:,i));
	% normalize Q(:,i)
	Q(:,i) = Q(:,i)/R(i,i);
end
