%
% Copyright Jack Poulson, 2014
%
% Returns the in-place QR factorization of a matrix A, where R is stored in
% the upper-triangle and the Householder vectors (with their first entry 
% implicitly set to one) stored in the strictly lower triangle.
% The Householder scalars are stored in t such that the j'th Householder 
% reflector is (I - t(j) v v'), where v = [1;A(j+1:end,j)].
%
function [A,t] = QRUnb(A)
[m,n]=size(A);
k=min(m,n);
t=zeros(k,1);
for j=1:k,
  [tau,v,beta]=Reflector(A(j:end,j));
  % Overwrite with (I - conj(tau) v v')A(j:end,j:end)
  %               = A(j:end,j:end) - conj(tau) v (v' A(j:end,j:end))
  % and then, since v is equal to 1 in its first entry, store v(2:end) 
  % below the diagonal of R
  A(j,      j) = beta;
  A(j+1:end,j) = v(2:end);
  t(j)=tau;

  z = v'*A(j:end,j+1:end);
  A(j:end,j+1:end) = A(j:end,j+1:end) - conj(tau)*v*z;
end
