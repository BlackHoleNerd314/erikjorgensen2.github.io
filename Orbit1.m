function [data] = Orbit1(m,a,X,V,J,n)
  for n0 = 1:n
    data = Orbit(m,a,X,V,J,1);
    X0 = data;
    V0 = X0-X;
    X = X0;
    V = V0;
  end
end
