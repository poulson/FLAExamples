%
% Copyright Jack Poulson, 2014
%
function [A,p] = LUPartialUnb(A)
[m,n]=size(A);
k=min(m,n);
p=zeros(k,1);
for j=1:k,
  % Find the maximum absolute value in this column of A(j:end,j:end) and
  % then swap with it. The reason for swapping the first j-1 columns of A
  % have to do with the fact that the permutation commutes with the last Gauss 
  % transform (you should prove this!).
  [val,idxOff]=max(abs(A(j:end,j)));  
  idx=idxOff+j-1;
  p(j)=idx;
  aTmp=A(idx,:);
  A(idx,:)=A(j,:);
  A(j,:)=aTmp;

  % Now perform a step of Gaussian elimination
  i1=j; j1=j; I2B=j+1; I2E=m; J2B=j+1; J2E=n;
  A(I2B:I2E,j1)=A(I2B:I2E,j1)/A(i1,j1);
  A(I2B:I2E,J2B:J2E)=A(I2B:I2E,J2B:J2E)-A(I2B:I2E,j1)*A(i1,J2B:J2E);
end
