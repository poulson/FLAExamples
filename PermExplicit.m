%
% Copyright Jack Poulson, 2014
%
% Returns the explicit n x n permutation matrix corresponding to
% the sequence of column swap specified by the pivot vector p.
%
function P = PermExplicit(p,n)

k=size(p);
if nargin < 2,
  n=k;
end

P=eye(n,n);
for j=1:k
  pTmp=P(:,p(j));
  P(:,p(j))=P(:,j);
  P(:,j)=pTmp;
end
