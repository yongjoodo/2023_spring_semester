function [x,iter] = iterimp(A,L,U,P,b,x0,tol,maxiter)
%iterref Iterative refinement of the solution to Ax = b
%obtained using LU factorization.
%  
%    [x,iter] = iterref(A,L,U,P,b,x0,tol,numitr)
%    iteratively computes succesive refinements
%    of an initial solution x0 of the system Ax = b obtained using the
%    LU factorization PA = LU. tol is the tolerance, and maxiter
%    is the user supplied number of iterations.
%    If the algorithm converges to the desired accuracy, iter contains
%    the number of iterations needed to converge. 
%    If the algorithm does not converge, iter = -1.
%    

	[m,n] = size(A);
	if m~=n
		disp('matrix A	is not square')  ;
		return;
	end
	
	xk = x0;
	for k = 1:maxiter	
		iter = k;	  %   always store iteration count in iter
		%  compute the residual vector
		r = b - A*xk;
		%  compute the correction
		c = lusolve(L, U, P, r);
		%  add the correction to form a new approximate solution
		x_kp1 = xk + c;
		%  see if we have attained the desired tolerance
		if norm(x_kp1 - xk) / norm(xk) < tol
			%  we have attained the tolerance
			x = x_kp1;
			return;
		end
		%  xk for the next iteration is x_kp1
		xk = x_kp1;
	end
	%  if we get here, the tolerance was not attained
	iter = -1;
	x = x_kp1;

