%
% Copyright Jack Poulson, 2014
%
% Returns [Q,R] such that A = Q R, where Q is unitary and R is upper-triangular
%
function [Q,R] = QRExplicit(A)
[m,n]=size(A);
k=min(m,n);

[A,t]=QRUnb(A);
R=triu(A(1:k,:));
Q=QRExpandQUnb(A,t);
