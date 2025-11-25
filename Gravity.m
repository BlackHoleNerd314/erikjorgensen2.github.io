function [data] = Gravity(m,a,X,V,n)

  for s=1:n
    x = X(2);
    y = X(3);
    z = X(4);
    r = sqrt(x^2+y^2+z^2);
    t = X(1);
    V(1) = 1;
    A = zeros(1,4);
    A(2) = -m*x/r^3;
    A(3) = -m*y/r^3;
    A(4) = -m*z/r^3;
    V = V+A;
    X = X+V;
    data(:,s)=X(:);
  endfor
end
