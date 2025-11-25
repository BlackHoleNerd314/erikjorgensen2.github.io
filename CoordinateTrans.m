function [data] = CoordinateTrans(a,r,theta,phi)

    n = (r-i*a)*exp(i*phi)*sin(theta);
    x = real(n);
    y = imag(n);
    z = r*cos(theta);

    data = [x,y,z];


end
