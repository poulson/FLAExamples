%
% Copyright Jack Poulson, 2014
%
% Returns y := inv(L) y, where L is lower-triangular
%
function y = TriangSolveLower(L,y)
n=size(y,1);
for j=1:n,
  y(j)=y(j)/L(j,j);
  y(j+1:end)=y(j+1:end)-L(j+1:end,j)*y(j);
end
