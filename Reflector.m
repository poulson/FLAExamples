%
% Copyright Jack Poulson, 2014
%
% Translation of include/elemental/lapack-like/util/Reflector.hpp
% from Elemental, which is loosely based upon LAPACK's zlarfg.f.
% The main difference is that this routine explicitly returns v with a '1'
% in the j'th entry and optionally also returns the nonzero entry of the result,
% 'beta'. After completion,
%
%  (eye(n,n)-conj(tau)*v*v')x = [zeros(j-1,1);beta;zeros(n-j,1)]
%
% which can, of course, be more efficiently applied as
%
%  x - conj(tau)*v*(v'*x)
%
function [tau,x,beta] = Reflector(x,j,safeInv)
if nargin < 1,
  error('Reflector requires at least one argument');
end
if nargin < 2,
   j = 1; 
end
if nargin < 3,
   % Ideally this should be safeMin/epsilon
   safeInv = 2e-292;
end

alpha=x(j);
nrmExl=norm(x([1:j-1 j+1:end]));
if nrmExl == 0 && imag(alpha) == 0,
  beta=-x(j);
  x(j)=1;
  tau=2;
  return;
end
if real(alpha) <= 0,
  beta = norm([alpha;nrmExl]);
else 
  beta = -norm([alpha;nrmExl]);
end

% Rescale if necessary
count=0;
if abs(beta) < safeInv,
  invOfSafeInv = 1/safeInv;
  while abs(beta) < safeInv,
    count=count+1;
    x([1:j-1 j+1:end])=x([1:j-1 j+1:end])*invOfSafeInv;
    alpha = alpha*invOfSafeInv;
    beta = beta*invOfSafeInv;
  end  
  nrmExl = norm(x([1:j-1 j+1:end]));
  if real(alpha) <= 0,
    beta = norm([alpha;nrmExl]);
  else
    beta = -norm([alpha;nrmExl]);
  end
end

tau = (beta-alpha)/beta;
x(j)=1;
x([1:j-1 j+1:end])=x([1:j-1 j+1:end])/(alpha-beta);

% Undo the scaling
for k=1:count,
  beta = beta*safeInv;
end
