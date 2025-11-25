function [data0] = LightCone(m,X,V,J)
  ma = PropBH(m);
  m = ma(1);
  a = ma(2);
% Extract Coordinates
  t = X(1);
  x = X(2);
  y = X(3);
  z = X(4);
% Define Radius
  R2 = x^2 + y^2 + z^2;
  Ra = R2 - a^2;
  AZ = a^2*z^2;
  b3 = sqrt(Ra^2 + 4*AZ);
  r2 = (Ra + b3)/2;
  r = sqrt(r2);
% Initialize
  n = 0;
  while t<r
    n = n + 1;
    % Extract Coordinates
    t = X(1);
    x = X(2);
    y = X(3);
    z = X(4);
    % Define Radius
    R2 = x^2 + y^2 + z^2;
    Ra = R2 - a^2;
    AZ = a^2*z^2;
    b3 = sqrt(Ra^2 + 4*AZ);
    r2 = (Ra + b3)/2;
    r = sqrt(r2);
    % Next Iteration
    data = Orbit(m,a,X,V,J,1);
    X0 = data;
    V0 = X0-X;
    X = X0;
    V = V0;
    % data tracking
    data0(n,:) = data;
  end

end
