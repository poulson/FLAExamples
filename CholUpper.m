%
% Copyright Jack Poulson, 2014
%
function A = CholUpper(A,nb)
n=size(A,1);
for j=1:nb:n,
  b=min(n-j+1,nb); 
  J1B=j;   J1E=j+b-1;
  J2B=j+b; J2E=n;
  A(J1B:J1E,J1B:J1E)=UpperCholUnb(A(J1B:J1E,J1B:J1E));
  A(J1B:J1E,J2B:J2E)=A(J1B:J1E,J1B:J1E)'\A(J1B:J1E,J2B:J2E);
  A(J2B:J2E,J2B:J2E)=A(J2B:J2E,J2B:J2E)-A(J1B:J1E,J2B:J2E)'*A(J1B:J1E,J2B:J2E);
end
A=triu(A);
