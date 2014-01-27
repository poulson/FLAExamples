%
% Copyright Jack Poulson, 2014
%

% Test for m > n
m=500;n=400;
fprintf('m=%d,n=%d\n',m,n);
A=randn(m,n)+i*randn(m,n); 
[L,U,P,Q]=LUFullExplicit(A);
err=norm(P*A*Q-L*U,'fro');
froA=norm(A,'fro');
relErr=err/froA;
fprintf('||A||_F              = %e\n',froA);
fprintf('||PAQ-LU||_F         = %e\n',err);
fprintf('||PAQ-LU||_F/||A||_F = %e\n',relErr);
fprintf('\n');

% Test for m < n
m=400;n=500;
fprintf('m=%d,n=%d\n',m,n);
A=randn(m,n)+i*randn(m,n); 
[L,U,P,Q]=LUFullExplicit(A);
err=norm(P*A*Q-L*U,'fro');
froA=norm(A,'fro');
relErr=err/froA;
fprintf('||A||_F              = %e\n',froA);
fprintf('||PAQ-LU||_F         = %e\n',err);
fprintf('||PAQ-LU||_F/||A||_F = %e\n',relErr);
fprintf('\n');
