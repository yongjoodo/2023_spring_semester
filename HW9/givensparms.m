function [c, s] = givensparms(xi,xj)
%givensparm Computes the parameters required to
%build a Givens rotation matrix.
%
%   [c s] = givensparam(xi,xj) computes the Givens parameters c and s
%   to use in building a Givens rotation using matrix values
%   xi and xj.

if xj == 0
   c = 1;
   s = 0;
elseif abs(xj) >= abs(xi)
   t = xi/xj;
   s = 1.0/sqrt(1 + t^2);
   c = s*t;
else
   t = xj/xi;
   c = 1.0/sqrt(1+t^2);
   s = c*t;
end
