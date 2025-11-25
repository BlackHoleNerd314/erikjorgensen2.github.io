function [data] = Orbit(m,a,X,V,J,n)

  constG0 = 2*m/((V(2)^2+V(3)^2+V(4)^2)*sqrt(X(2)^2+X(3)^2+X(4)^2));
  constc0 = 1/sqrt(V(2)^2+V(3)^2+V(4)^2);


  c = 299792458;
  h = 6.62607015/10^34;
  K = (c^2)/(10^7);
  e = 1.602176634/10^19;



  acc0 = (K*e^2)/((h/(2*pi))*c);

  if (constG0>acc0)&&(constc0<(1/acc0))
    data = Geodesic1(m,a,X,V,J,n);
  end
  if (constG0<acc0)&&(constc0<(1/acc0))
    data = Relativity(m,a,X,V,n);
  end
  if (constG0>acc0)&&(constc0>(1/acc0))
    data = Gravity(m,a,X,V,n);
  end
  if (constG0<acc0)&&(constc0>(1/acc0))
    data = Unity(m,a,X,V,n);
  end

  data = data';
endfunction
