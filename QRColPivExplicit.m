%
% Copyright Jack Poulson, 2014
%
function [Q,R,P] = QRColPivExplicit(A,numSteps)
[m,n]=size(A);
k=min(m,n);
if nargin < 2,
  numSteps = k;
end

[A,t,p]=QRColPiv(A,numSteps);
R=triu(A(1:numSteps,:));
A=A(:,1:numSteps);
Q=QRExpandQUnb(A,t);
P=PermExplicit(p,n);
