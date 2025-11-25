function [M] = Lorentz(k,j)
  J = zeros(4,4);
  K = zeros(4,4);
  for n=1:3
    n1 = mod(n-2,3)+1;
    n0 = mod(n,3)+1;
    K(1,n+1) = k(n);
    J(n1+1,n0+1) = -j(n);
    K(n+1,1) = k(n);
    J(n0+1,n1+1) = j(n);
  endfor
  M = expm(K+J);
endfunction
