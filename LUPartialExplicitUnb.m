%
% Copyright Jack Poulson, 2014
%
% Returns [L,U,P] such that P A = L U
%
function [L,U,P] = LUPartialExplicitUnb(A)
[m,n]=size(A);
k=min(m,n);
[A,p]=LUPartialUnb(A);
L=tril(A(:,1:k),-1)+diag(ones(k,1),m,k);
U=triu(A(1:k,:));
P=PermExplicit(p,m)';
