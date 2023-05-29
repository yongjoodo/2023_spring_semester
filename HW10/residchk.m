function average = residchk(A,V,D,avgonly)
%RESIDCHK View the residuals after computation of
%the eigenvalues/eigenvectors of a large sparse matrix
%  Input:
%     A is the sparse matrix.
%     V, D are the result of [V,D] = eigsb(...) or
%     [V, D] = eigssymb(...)
%  Output:
%     Listing of residual for each of the eigenpairs computed
%     D(i,i)/V(:,i).
%     average, the average of the residuals.

average = 0.0;
n = size(D,1);

if nargin == 3
   avgonly = false;
end

for i = 1:n
   residual = norm(A*V(:,i) - D(i,i)*V(:,i));
   average = average + residual;
   if avgonly == false
      fprintf('Eigenpair %d residual = %g\n',i,residual);
   end
end

average = average/n;