function [position] = CGH4(M,t,x,y,z,vx,vy,vz,n,c,G,h,dt)


for s = 1:n;
  r = sqrt(x^2+y^2+z^2);
  t0 = t - r/c;
  v = sqrt(vx^2+vy^2+vz^2)/c;
  dr = sqrt(h*dt/M);
  x = x + sinh(v)*vx*dt/v + randn(1)*dr;
  y = y + sinh(v)*vy*dt/v + randn(1)*dr;
  z = z + sinh(v)*vz*dt/v + randn(1)*dr;
  r1 = sqrt(x^2+y^2+z^2);
  vx = vx - G*M*x*dt/r1^3;
  vy = vy - G*M*y*dt/r1^3;
  vz = vz - G*M*z*dt/r1^3;
  t0 = t0 + cosh(v)*dt;
  t = t0 + r1/c;
  position(s,:) = [t,x,y,z];
end
















end
