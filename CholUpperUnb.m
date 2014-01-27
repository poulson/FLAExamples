%
% Copyright Jack Poulson, 2014
%
function A = CholUpperUnb(A)
n=size(A,1);
for j=1:n,
  j1=j; J2B=j+1; J2E=n;
  A(j1,j1)=sqrt(A(j1,j1));  
  A(j1,J2B:J2E)=A(j1,J2B:J2E)/A(j1,j1);
  A(J2B:J2E,J2B:J2E)=A(J2B:J2E,J2B:J2E)-A(j1,J2B:J2E)'*A(j1,J2B:J2E);
end
A=triu(A);
