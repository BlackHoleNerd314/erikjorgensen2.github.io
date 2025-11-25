function [position] = CGH(M,t,x,y,z,vx,vy,vz,s0)

for s = 1:s0;
  v = sqrt(vx^2+vy^2+vz^2);
  dr = sqrt(1/M);
  x = x + sinh(v)*vx/v + randn(1)*dr;
  y = y + sinh(v)*vy/v + randn(1)*dr;
  z = z + sinh(v)*vz/v + randn(1)*dr;
  r = sqrt(x^2+y^2+z^2);
  vx = vx - M*x/r^3;
  vy = vy - M*y/r^3;
  vz = vz - M*z/r^3;
  t = t + cosh(v);
  t0 = t + r;
  position(s,:) = [t0,x,y,z];
end
%


end
