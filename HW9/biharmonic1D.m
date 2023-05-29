function [x, y, residual, iter] = biharmonic1D(n, f, tol, maxiter) 
% BIHARMONIC1D Approximate and plot the solution to the biharmonic % equation y’’’’ = f(x), 0 <= x <= 1, y(0) = y(1) = 0,
% y’’(0) = y’’(1) = 0. Solve the system of equations using the
% preconditioned CG method.
h = 1/n;
a = ones(n-3,1);
b = -4*ones(n-2,1);
c = 6*ones(n-1,1);
d = b;
e = a;
B = diag(a,-2) + diag(b,-1) + diag(c) + diag(d,1) + diag(e,2); B(1,1) = 5;
B(n-1,n-1) = 5;
B = sparse(B);
x = (0.0:h:1.0)';
y = zeros(n+1,1);
b = h^4*feval(f,x(2:n));
x0 = zeros(n-1,1);
[y(2:n), residual, iter] = precg(B,b,x0,tol,maxiter);
plot(x,y,'color','k','LineWidth',1.5);
xlabel('x');
ylabel('y(x)');
title({'$$\frac{d^{4}u}{dx^{4}}=x,\,0\leq x\leq1\\$$'; '$u\left(0\right)=u\left(1\right)=\frac{d^{2}u}{dx^{2}} \left(0\right)=\frac{d^{2}u}{dx^{2}}=0,\,0\leq x\leq1$'},...
'fontsize', 14, 'FontWeight','bold',... 
'interpreter','latex');
