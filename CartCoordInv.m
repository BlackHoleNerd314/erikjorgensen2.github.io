function [data] = CartCoordInv(a,x,y,z)

    % Define Radius
    R2 = x^2 + y^2 + z^2;
    Ra = R2 - a^2;
    AZ = a^2*z^2;
    b3 = sqrt(Ra^2 + 4*AZ);
    r2 = (Ra + b3)/2;
    r = sqrt(r2);
    % Angles
    theta = acos(z/r);
    phi = atan2(y,x) + atan2(a,r);
    %
    data = [r,theta,phi];


end
