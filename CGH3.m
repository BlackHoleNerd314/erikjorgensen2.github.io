function [position] = CGH3(M,t,x,y,z,vx,vy,vz,s0)
% CGH2: Simulates a quantum-relativistic test particle in a gravitational field.
% Assumptions incorporated:
% 1) Diffusion constant D = hbar/(2*m) → spatial step sqrt(1/M)
% 2) Relativistic motion via rapidity (sinh/cosh)
% 3) Equivalence principle: M=m in observer's rest frame
% 4) In source rest frame, relativistic retardation vanishes
% 5) 3D entanglement modeled via Gaussian diffusion (central limit theorem)
% 6) Retarded time simplified in observer's rest frame

for s = 1:s0
    % Compute relativistic speed magnitude (rapidity parameter)
    v = sqrt(vx^2 + vy^2 + vz^2);  % Test particle is relativistic only
    
    % Spatial diffusion step (quantum effect)
    dr = sqrt(1/M);  % sqrt(D*dt) dimensionally gives meters
    
    % Relativistic position update + stochastic diffusion (field collapse)
    x = x + sinh(v)*vx/v + randn(1)*dr;
    y = y + sinh(v)*vy/v + randn(1)*dr;
    z = z + sinh(v)*vz/v + randn(1)*dr;
    
    % Radial distance from gravitational source
    r = sqrt(x^2 + y^2 + z^2);
    
    % Gravitational acceleration (Newtonian form; valid in source rest frame)
    vx = vx - M*x/r^3;
    vy = vy - M*y/r^3;
    vz = vz - M*z/r^3;
    
    % Relativistic time step (proper time increment)
    t = t + cosh(v);
    
    % Retarded time for light travel from particle to observer
    t0 = t + r;  % simplified because observer is in source rest frame
    
    % Store position and retarded time
    position(s,:) = [t0, x, y, z];
end
%
end
%


































