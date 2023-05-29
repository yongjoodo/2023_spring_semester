function varargout = eigssymb1(A, nev, m, tol, maxiter,reortho)
%EIGSSYMB executes the Lanczos process with implicit restarts to
%compute a specified number of the largest eigenvalues of the sparse
%symmetric matrix A.
%
%   [V D] = eigssymb(A,v0,nev,m,tol,maxiter) returns a matrix V of
%   eigenvectors and diagonal matrix D of eigenvalues for the sparse
%   symmetric A. m is the size of the Krylov subspace to be
%   used and should be as small as possible. nev is the number of
%   eigenvalues desired and should also be small. The error tolerance
%   tol defaults to 1.0e-6, maxiter defaults to 100.
%
%   E = eigssymb(A, nev, m, tol, maxiter) - assigns E a column
%   vector containing the eigenvalues. The default values of tol and
%   maxiter are the same as the previous calling sequence.
%   reortho chooses whether reorthogonalization at the Lanczosf.m

% initialization
n = size(A, 1);

if nargin < 3
   error('eigssymb requires at least 3 arguments.');
end

% take care of defaults
if nargin == 3
	tol = 1.0e-6;
	maxiter = 100;
elseif nargin == 4
	maxiter = 100;
end

% we use deflation as k increases
k = 1;
I = eye(n);
em = zeros(m,1);
em(m) = 1;
ek = zeros(nev,1);
ek(nev) = 1;

V = zeros(n, m);
T = zeros(m, m);
% start with random v
v = randn(n,1);
f = zeros(n,1);
eigenvalues = zeros(nev,1);
%peform an m-step Lanczos decomposition
[V,T,f] = lanczosf(A,V,T,f,k,v,m,reortho);
% we only allow maxiter iterations
iter = 0;                                                               

while true
   iter = iter + 1;

% find eigenvalues and eigenvectors for T.
   [UT, DT] = eig(T);
	sigma = diag(DT);

% reorder the eigenvalues and the corresponding eigenvectors.
	[~,ind] = sort(abs(sigma),'descend');
	sigma = sigma(ind);

   VTtmp = UT;
   for j =1:nev
      UT(:,j) = VTtmp(:, ind(j));
   end

   Q = eye(m);
   % use sigma(nev+1), ..., sigma(m) as the shifts for
   % the QR iteration
   j = m;
   while j >= nev+1
      lambda = sigma(j);
      % lambda is real. use an implicit single shift
      
      % NOTE: If you don't want to use chase modified for this
      % iteration, the same thing can be accomplished with the
      % statements
      
%       [Qj,~] = qr(T - lambda*eye(m));
%       j = j - 1;
%       T = Qj'*T*Qj;
      
      [Qj,T] = chase2(T,lambda);
      j = j-1;
      Q = Q*Qj;
    end
    
    % compute the residual norm for the kth eigenpair
    u = UT(:,k);
    resid_norm = norm(f)*abs(u(m));
    % lock v_k if the tolerance is obtained
    if(resid_norm < tol)
      if(k <= nev)
         eigenvalues(k) = sigma(k);
         if k ~= nev
            k = k+1;
         else
            break;
         end
		end
    end
    
    % build an m-step factorization from the nev step one
    betak = T(nev+1,nev);
    sigmak = Q(m,nev);
    fk = V(:,nev+1)*betak + f*sigmak;
    V(:,1:nev) = V(:,1:m)*Q(:,1:nev);
    [V, T, f] = lanczosf(A,V(:,1:nev),T(1:nev,1:nev),fk,nev+1,V(:,nev),m);
   
	if(iter >= maxiter)
      error('Cannot compute the requested eigenvalues in the specfied number of iterations');
	end
end

if nargout == 2
   V(:,1:nev) = V(:,1:m)*UT(:,1:nev);
	varargout{1} =  V(:,1:nev);
	varargout{2} =  diag(eigenvalues);
else
   varargout{1} = eigenvalues;
end