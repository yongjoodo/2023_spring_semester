function A = givensmult(A,i,j,c,s)
% GIVENSMULT Form A*J(i,j,c,s)

a = A(:,i);
b = A(:,j);
A(:,i) = c*a + s*b;
A(:,j) = -s*a + c*b;
