function [data0] = ReferenceBody(M,X,V,A,x,v,j)
  % relative spacetime position and initial normalized velocity:
  x = x-X;
  X0 = X;
  G = Metric(0,0,3,4,12);
  m2 = v*G*v';
  v = v/sqrt(-m2);
%% Lorentz Boost scene such that big M is stationary:
  G = Metric(0,0,3,4,12);
  M2 = V*G*V';
  V = V/sqrt(-M2);
  K0 = [V(2),V(3),V(4)];
  if norm(K0)==0
    boostT = eye(4);
  else
    K = (acosh(V(1)))*K0/norm(K0);
    boostT = Lorentz(-K,[0,0,0]);
  end
  A0 = boostT*A';
  j0 = boostT*j';
  x = boostT*x';
  x = x';
  v = boostT*v';
  v = v';
  A = [A(2),A(3),A(4)];
%% Rotate scene such that big M has spin axis align with Z axis:
  if A(1)^2+A(2)^2==0
    rotZ = eye(4);
  else
    A = A/norm(A);
    J = cross(A,[0,0,1]);
    if A(3)<0
      J = (pi-asin(norm(J)))*J/norm(J);
    else
      J = asin(norm(J))*J/norm(J);
    end
    rotZ = Lorentz([0,0,0],-J);
  end
  A = [0,A/norm(A)];
  A0 = rotZ*A';
  j1 = rotZ*j0';
  x = rotZ*x';
  x = x';
  v = rotZ*v';
  v = v';
%% Calculate trajectory of little m around big M:
  data = LightCone(M,x,v,j1);
%% Undo rotation:
  for s = 1:length(data)
    data(s,:) = (inv(rotZ)*data(s,:)')';
  end
%% Undo boost:
  for s = 1:length(data)
    data(s,:) = (inv(boostT)*data(s,:)')';
  end
%% Undo translation:
  for s = 1:length(data)
    data(s,:) = data(s,:)+X0;
  end
%% Collect data:
  data0 = data;

end
