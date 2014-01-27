%
% Copyright Jack Poulson, 2014
%
% Returns the explicit n x n permutation matrix corresponding to
% the sequence of column swaps specified by the pivot vector p.
% If row swaps were desired, then the result should be transposed.
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
