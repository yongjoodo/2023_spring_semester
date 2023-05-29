function x = cholsolve(R,b)
% cholsolve  Solve a system after cholesky factorization.
%
%   x = cholsolve(R,b) solves the system Ax = b, given
%   upper triangular matrix R obtained from the Cholesky
%   factorization of positive definite matrix A.

[m, n] = size(R);
if m ~= n
	disp('The system is not square.');
	return;
end

% Solve the lower triangular system R'y = b
y = forsolve(R', b);
% Solve the upper triangular system Rx = y
% by back substitution
x = backsolve(R, y);
