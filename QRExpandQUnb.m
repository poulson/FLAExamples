%
% Copyright Jack Poulson, 2014
%
function Q = QRExpandQUnb(A,t)
[m,n]=size(A);
k=min(m,n);
Q=eye(m,k);
for j=k:-1:1,
  tau=t(j);
  v=[1;A(j+1:end,j)];
  % conj(H) = I - conj(tau)*v*v'
  % We only need to apply to the j:end columns of Q
  z = v'*Q(j:end,j:end);
  Q(j:end,j:end)=Q(j:end,j:end)-conj(tau)*v*z;
end
