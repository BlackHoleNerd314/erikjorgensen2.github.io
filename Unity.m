function [data] = Unity(m,a,X,V,n)
for s=1:n
    x = X(2);
    y = X(3);
    z = X(4);
    r = sqrt(x^2+y^2+z^2);
    t = X(1);
    V(1) = 1;
    A = zeros(1,4);
    V = V+A;
    X = X+V;
    data(:,s)=X(:);
  endfor
end
