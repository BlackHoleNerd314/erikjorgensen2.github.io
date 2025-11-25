function [data] = Geodesic0(m,a,X,V,n)


  for s=1:n
    x = X(2);
    y = X(3);
    z = X(4);
    G = KerrMetric(m,a,x,y,z);
    L = Connection(m,a,x,y,z);
    m2 = V*G*V';
    V = V/sqrt(-m2);
    A = zeros(1,4);
    for u = 1:4
      for v = 1:4
        A(:) = A(:) +L(:,u,v)*V(u)*V(v);
      endfor
    endfor
    V = V+A;
    X = X+V;
    data(:,s)=X(:);
  endfor

endfunction
