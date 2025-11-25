function [data] = Relativity(m,a,X,V,n)
  for s=1:n
    x = X(2);
    y = X(3);
    z = X(4);
    G = Metric(0,0,3,4,12);
    m2 = V*G*V';
    V = V/sqrt(-m2);
    A = zeros(1,4);
    V = V+A;
    X = X+V;
    data(:,s)=X(:);
  endfor
end
