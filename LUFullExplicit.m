%
% Copyright Jack Poulson, 2014
%
% Returns [L,U,P,Q] such that P A Q = L U
%
function [L,U,P,Q] = LUFullExplicit(A)
[m,n]=size(A);
k=min(m,n);
[A,p,q]=LUFull(A);
L=tril(A(:,1:k),-1)+diag(ones(k,1),m,k);
U=triu(A(1:k,:));
P=PermExplicit(p,m)';
Q=PermExplicit(q,n);
