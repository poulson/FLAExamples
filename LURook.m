%
% Copyright Jack Poulson, 2014
%
function [A,p,q] = LURook(A)
[m,n]=size(A);
k=min(m,n);
p=zeros(k,1);
q=zeros(k,1);
count=0;
for j=1:k,
  % Perform rook pivoting 
  done=false; currCol=j;
  while ~done,
    count=count+1;
    [val,iOff]=max(abs(A(j:end,currCol)));
    iIdx=iOff+j-1;
    [val,jOff]=max(abs(A(iIdx,j:end)));
    jIdx=jOff+j-1;
    if jIdx == currCol,
      done=true;
    else
      currCol=jIdx;
    end
  end
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
avgSteps=count/k;
fprintf('Average number of rook pivoting steps: %f\n',avgSteps);
