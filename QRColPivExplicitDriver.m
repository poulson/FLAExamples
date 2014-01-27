%
% Copyright Jack Poulson, 2014
%

% Test for m > n
m=500;n=400;r=10;
fprintf('m=%d,n=%d,r=%d\n',m,n,r);
U=randn(m,r)+i*randn(m,r); V=randn(n,r)+i*randn(n,r); A=U*V';
[Q,R,P]=QRColPivExplicit(A,r);
err=norm(A*P-Q*R,'fro');
froA=norm(A,'fro');
relErr=err/froA;
fprintf('||A||_F             = %e\n',froA);
fprintf('||AP-QR||_F         = %e\n',err);
fprintf('||AP-QR||_F/||A||_F = %e\n',relErr);
fprintf('\n');

% Test for m < n
m=400;n=500;r=10;
fprintf('m=%d,n=%d,r=%d\n',m,n,r);
U=randn(m,r)+i*randn(m,r); V=randn(n,r)+i*randn(n,r); A=U*V';
[Q,R,P]=QRColPivExplicit(A,r);
err=norm(A*P-Q*R,'fro');
froA=norm(A,'fro');
relErr=err/froA;
fprintf('||A||_F             = %e\n',froA);
fprintf('||AP-QR||_F         = %e\n',err);
fprintf('||AP-QR||_F/||A||_F = %e\n',relErr);
fprintf('\n');
