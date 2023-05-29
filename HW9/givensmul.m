function A = givensmul(A,i,j,c,s)
%givensmul Computes JA, where J is the
%Givens rotation with parameters c and s.
%
%   A = givensmul(A,i,j,c,s) computes the product J(i,j,c,s)*A,
%   where J is a Givens rotation that zeros-out A(j,i)

a = A(i, :);
b = A(j, :);
A(i,:) = c*a + s*b;
A(j,:) = -s*a + c*b;

