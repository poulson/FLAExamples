%
% Copyright Jack Poulson, 2014
%
% Returns y := inv(U) y, where U is upper-triangular
%
function y = TriangSolveUpper(U,y)
n=size(y,1);
for j=n:-1:1,
  y(j)=y(j)/U(j,j);
  y(1:j-1)=y(1:j-1)-U(1:j-1,j)*y(j);
end
