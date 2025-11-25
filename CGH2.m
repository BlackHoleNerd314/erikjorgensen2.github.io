function [position] = CGH2(M,t,x,y,z,vx,vy,vz,s0)

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
%Please remember these points about my code: 1). the diffusion constant hbar/(2*m) has dimensions of meters^2 per second. so to get a spacial displacement one must take its square root. sqrt(meters^2) = sqrt(1/m) => meters = deltaX. The test particle is quantum only 2). Lorentz boosting a test particle in special relativity does indeed use rapidity as its hyperbolic angle parameter. Slopes are not invariant under rotation, but angles are. The test particle is relativistic only 3). The test particle is gravitating only, but if we frame hop into the observers frame, the scene is symmetric between frames, and in this rest frame it seems as if the gravitational mass is produced by the observer, which is why we set M=m and the equivalence principle holds. 4). In the rest frame of the source, the relativistic effects due to retardation vanish, its just a geodesic, and the painleve gullstrand tetrad is reproduced, then due to the escape velocity being galilean boost of a lorentz frame, there is vanishing temporal acceleration, and the acceleration reduced back to newtonian. 5). Entanglement between dimensions in the schrodinger/heat equations works exactly the same as the probability density function for a 3d random walk. in 3d, the central limit theorem applies such that the exponential for the probability is quadratic and as such the entaglement holds in each dimension independently due to coorelations in all three dimensions producing a sum over squares in the probability phase/damping, which works due to the 3d euclidean metric 6). For similar reasons, the  retarded time in the observer's frame which is a rest makes the retarded time simple. This is due to the fact that no lorentz boosting between events on the obervers worldline is needed, as it and the test particle are following geodesics. Then the retarded time simplifies drastically, because the spacial distances are not length contracted and the retarded time is simply the light travel distance, of which the time is also simulantious due to trivial lorentz transforms. 

















end
