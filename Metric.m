function [G]=Metric(m,a,x,y,z)
    % Standard Minkowski Spacetime
    N = diag([-1,ones(1,3)]);
    % Define Radius
    R2 = x^2 + y^2 + z^2;
    Ra = R2 - a^2;
    AZ = a^2*z^2;
    b3 = sqrt(Ra^2 + 4*AZ);
    r2 = (Ra + b3)/2;
    r = sqrt(r2);
    % Define Scalar Perturbation
    H1 = 2*m*r^3;
    H2 = r^4 + AZ;
    H = H1/H2;
    % Define Vector Perturbation
    Lt = 1;
    Lx = (r*x + a*y)/(r^2+a^2);
    Ly = (r*y - a*x)/(r^2+a^2);
    Lz = z/r;
    L = zeros(1,4);
    L(2) = Lx;
    L(3) = Ly;
    L(4) = Lz;
    L(1) = Lt;
    % Construct Kerr Spacetime
    G = N + H*L'*L;

end
