function [] = TempTest1(z)
u = z;
du = 0.01;
eps1 = 10000;
eps0 = 0.01;
gpoints = 1;
E1 = 0;
N1 = 0;
e1 = 0;
n1 = 0;
zeta = 0;
for iter0 = 1:1000000
  zeta = zeta+1/iter0^3;
end
for iter0 = 1:eps1;
    u = u + du;
    if u > eps1;
        E0 = 0;
    else;
        E0 = (15 * gpoints / pi ^ 4) * (u ^ 2 * sqrt(u ^ 2 - z ^ 2)) / (exp(u) - 1);
        e0 = (15 * gpoints / pi ^ 4) * (u ^ 2 * sqrt(u ^ 2 - z ^ 2)) / (exp(u) + 1);
    end;
    N0 = (E0 / u);
    n0 = (e0 / u);
    E1 = E1 + E0 * du;
    N1 = N1 + N0 * du;
    e1 = e1 + e0 * du;
    n1 = n1 + n0 * du;
end;
if z<eps0
  E1 = gpoints;
  N1 = gpoints*((2*zeta)/(pi^4/15));
  e1 = gpoints*(7/8);
  n1 = gpoints*(3/4)*2*zeta/(pi^4/15);
end
disp([E1,N1,e1,n1])
  

end;
%




























































































