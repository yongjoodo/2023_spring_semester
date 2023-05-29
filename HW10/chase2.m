function varargout = chase2(T,sigma)
% CHASE Bulge chase in the symmetric Francis algorithm.
%
%   [P T1] = chase(T) returns an orthogonal matrix P such that P'*T*P = T1,
%   where T and T1 are symmetric tridiagonal matrices, and P is
%   an orthogonal matrix.
%
%   T1 = chase(T) returns the result of a bulge chase down the symmetric
%   tridiagonal matrix T.
%

n = size(T,1);

% handle (1,2) and (2,1) with shift sigma
[c,s] = givensparms(T(1,1)-sigma,T(2,1));
T = givensmul(T,1,2,c,s);
T = givensmult(T,1,2,c,s);
if nargout == 2
	P = eye(n);
	P = givensmul(P,1,2,c,s);
end

% chase the bulge
for i = 1:n-2
   [c,s] = givensparms(T(i+1,i),T(i+2,i));
   T = givensmul(T,i+1,i+2,c,s);
   T = givensmult(T,i+1,i+2,c,s);
	T(i+2,i) = 0;
	if nargout == 2
		P = givensmul(P,i+1,i+2,c,s);
	end
end

if nargout == 2
	varargout{1} = P';
	varargout{2} = T;
else
	varargout{1} = T;
end
	
