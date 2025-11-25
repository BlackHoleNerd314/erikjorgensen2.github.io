function [e0] = Frame(m,a,x,y,z)

    a = -a;

    data0 = CartCoordInv(a,x,y,z);


    r = data0(1);
    theta = data0(2);
    phi = data0(3);






    L = [0,-1,0,0];
    N = [1,(1-(2*m*r)/(r^2+a^2*cos(theta)^2))/2,0,0];
    M = [i*a*sin(theta),i*a*sin(theta),1,i/(sin(theta))]/(sqrt(2)*(r-i*a*cos(theta)));
    Mbar = conj(M);

    e0(1,:) = L;
    e0(2,:) = N;
    e0(3,:) = Mbar;
    e0(4,:) = M;

    nullretard = [1,0,0,0;-1,1,0,0;0,0,1,0;0,0,0,1];

    e0 = e0*nullretard;


    d = 0.0000001*sqrt(x^2+y^2+z^2);;%d = 1/2^32;

    xyz = CoordinateTrans(-a,r,theta,phi);

    xyz_r = CoordinateTrans(-a,r+d,theta,phi);
    xyz_theta = CoordinateTrans(-a,r,theta+d,phi);
    xyz_phi = CoordinateTrans(-a,r,theta,phi+d);

    d_dr = (xyz_r-xyz)/d;
    d_dtheta = (xyz_theta-xyz)/d;
    d_dphi = (xyz_phi-xyz)/d;

    jacobian = [1,0,0,0;0,d_dr;0,d_dtheta;0,d_dphi];


    e0 = e0*jacobian;




    %n = (r-i*a)*exp(i*phi)*sin(theta);


    %g0 = (Metric(m,-a,x,y,z));

    %data = g;







end
