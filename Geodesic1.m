function [data] = Geodesic1(m,a,X,V,J,n)


  for s=1:n
    x = X(2);
    y = X(3);
    z = X(4);
    G = Metric(m,a,x,y,z);
    L = Connection0(m,a,x,y,z);
    m2 = V*G*V';
    V = V/sqrt(-m2);
    j2 = J*G*J';
    J = J/sqrt(j2);
    A = zeros(1,4);
    dJ = zeros(1,4);
    for u = 1:4
      for v = 1:4
        A(:) = A(:) +L(:,u,v)*V(u)*V(v);
        dJ(:) = dJ(:) +L(:,u,v)*V(u)*J(v);
      endfor
    endfor
    J = J+dJ;
    V = V+A;
    X = X+V;
    data(:,s)=[X(:);V(:);J(:)];
  endfor

endfunction
