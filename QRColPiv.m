%
% Copyright Jack Poulson, 2014
% 
% Returns the in-place QR factorization, with R stored in the upper-triangle
% and the Householder vectors (with their first entry implicitly equal to one)
% stored in the strictly lower triangle. The Householder scalars are stored in 
% t such that the j'th Householder reflector is (I - t(j) v v'), where
% v = [1;A(j+1:end,j)]. The column swaps are stored in the pivot vector p.
%
function [A,t,p] = QRColPiv(A,numSteps)
[m,n]=size(A);

k=min(m,n);
if nargin < 2,
  numSteps = k;
end

t=zeros(numSteps,1);
p=zeros(numSteps,1);

% Initialize the column norms
colNorms=zeros(n,1);
for j=1:n,
  colNorms(j)=norm(A(:,j));  
end

% Perform numSteps iterations of Businger-Golub pivoted QR
for j=1:numSteps,
  % Swap the maximum column norm into the first position
  % (and then swap the column norms)
  [val,idxOff]=max(colNorms(j:end));
  idx=idxOff+j-1;
  p(j)=idx;
  aTmp=A(:,idx);        A(:,idx)=A(:,j);           A(:,j)=aTmp;
  nrmTmp=colNorms(idx); colNorms(idx)=colNorms(j); colNorms(j)=nrmTmp; 

  % Perform a step of pivoted QR
  [tau,v,beta]=Reflector(A(j:end,j));
  A(j,      j) = beta;
  A(j+1:end,j) = v(2:end);
  t(j)=tau;
  z = v'*A(j:end,j+1:end);
  A(j:end,j+1:end) = A(j:end,j+1:end) - conj(tau)*v*z;

  % Recompute the column norms (TODO: Use a safe updating formula)
  for l=j+1:n,
    colNorms(l)=norm(A(j+1:end,l));
  end
end
