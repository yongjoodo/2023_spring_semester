function [Q,R] = modqrgrsch(A)
%modqrgrsh Modified Gram-Schmidt method for QR factorization
%   [Q,R] = modgrsch computes the reduced QR factorization
%   of an m x n matrix A using the modified Gram-Schmidt method.
%   A = QR, where R is n x n upper triangular, and Q is m x n
%   with orthonormal columns

   [m, n] = size(A);
   % allocate the output matrices Q and R
   Q = zeros(m, n);
   R = zeros(n, n);

   % begin the computation
   for i=1:n
      Q(:,i) = A(:,i);
      for j=1:i-1
         R(j,i) = Q(:,j)'* Q(:,i);
         Q(:,i) = Q(:,i) - R(j,i)*Q(:,j);
      end
      % complete column i of R by computing the diagonal element
      R(i,i) = norm(Q(:,i));
      % normalize Q(:,i)
      Q(:,i) = Q(:,i)/R(i,i);
   end