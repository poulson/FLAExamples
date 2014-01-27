%
% Copyright Jack Poulson, 2014
%
function [A,p,q] = LUFull(A)
[m,n]=size(A);
k=min(m,n);
p=zeros(k,1);
q=zeros(k,1);
for j=1:k,
  % Find the maximum absolute value in A(j:end,j:end) and
  % then swap with it. The reason for swapping the first j-1 columns 
  % of A has to do with the fact that the row permutation commutes with the 
  % last Gauss transform (you should prove this!).
  [val,idxOff]=max(abs(A(j:end,j:end))(:));  
  iIdx = j + mod(idxOff-1,m-j+1);
  jIdx = j + floor((idxOff-1)/(m-j+1));
  p(j)=iIdx;
  q(j)=jIdx;
  % Perform the row swap
  aTmp=A(iIdx,:);A(iIdx,:)=A(j,:);A(j,:)=aTmp;
  % Perform the column swap
  aTmp=A(:,jIdx);A(:,jIdx)=A(:,j);A(:,j)=aTmp;

  % Now perform a step of Gaussian elimination
  i1=j; j1=j; I2B=j+1; I2E=m; J2B=j+1; J2E=n;
  A(I2B:I2E,j1)=A(I2B:I2E,j1)/A(i1,j1);
  A(I2B:I2E,J2B:J2E)=A(I2B:I2E,J2B:J2E)-A(I2B:I2E,j1)*A(i1,J2B:J2E);
end
