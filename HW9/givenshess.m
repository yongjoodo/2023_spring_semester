function [Q, R] = givenshess(H)
%givenshess reduced QR factorization of a (k+1) x k
%upper Hessenberg matrix by Givens rotations.
%
%   [Q,R] = givenshess(H) returns a reduced QR factorization
%   of the (k+1) x k upper Hessenberg matrix H. Q is a
%   (k+1) x k orthogonal matrix, and R is a k x k upper
%   triangular matrix such that H = QR.

[kp1,k] = size(H);
Q = eye(kp1,kp1);

for i= 1:k
   [c,s] = givensparms(H(i,i),H(i+1,i));
   H = givensmul(H,i,i+1,c,s);
   Q = givensmul(Q,i,i+1,c,s);
end

R = H(1:k,1:k);
Q = Q';
Q = Q(:,1:k);
