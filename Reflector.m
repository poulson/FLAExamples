%
% Copyright Jack Poulson, 2014
%
% Translation of include/elemental/lapack-like/util/Reflector.hpp
% from Elemental, which is loosely based upon LAPACK's zlarfg.f.
% The main differences are: 
%   1. Returns v with a '1' in the j'th entry,
%   2. Legal (complex) Householder reflectors are always used, e.g.,
%      if x is already proportional to ej, then the result is -x.
%   3. The reflector is defined as (I-tau v v') rather than (I-conj(tau) v v')
%
% More specifically, this routine returns [tau,v,beta] such that
%
%  (eye(n,n)-tau*v*v')x = [zeros(j-1,1);beta;zeros(n-j,1)]
%
% which can, of course, be more efficiently applied as
%
%  x - tau*v*(v'*x)
%
function [tau,v,beta] = Reflector(x,j,safeInv)
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

v=x;
if nrmExl == 0 && imag(alpha) == 0,
  beta=-x(j);
  v(j)=1;
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
    v([1:j-1 j+1:end])=x([1:j-1 j+1:end])*invOfSafeInv;
    alpha = alpha*invOfSafeInv;
    beta = beta*invOfSafeInv;
  end  
  nrmExl = norm(v([1:j-1 j+1:end]));
  if real(alpha) <= 0,
    beta = norm([alpha;nrmExl]);
  else
    beta = -norm([alpha;nrmExl]);
  end
end

tau = (beta-conj(alpha))/beta;
v(j)=1;
v([1:j-1 j+1:end])=v([1:j-1 j+1:end])/(alpha-beta);

% Undo the scaling
for k=1:count,
  beta = beta*safeInv;
end
