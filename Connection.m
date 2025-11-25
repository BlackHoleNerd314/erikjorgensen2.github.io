function [L12]=Connection(m,a,x,y,z)
    d = 0.0001*sqrt(x^2+y^2+z^2);%d = 1/2^8*(x^2+y^2+z^2)/m^2;
    dt = d;
    dx = d;
    dy = d;
    dz = d;
    G = KerrMetric(m,a,x,y,z);
    g = inv(G);
    %Gt = G;KerrMetric(m,a,x,y,z);
    Gx = KerrMetric(m,a,x+dx,y,z);
    Gy = KerrMetric(m,a,x,y+dy,z);
    Gz = KerrMetric(m,a,x,y,z+dz);
    dG = zeros(4,4,4);
    %dG(:,:,1) = 0;
    dG(:,:,2) = (Gx-G)/dx;
    dG(:,:,3) = (Gy-G)/dy;
    dG(:,:,4) = (Gz-G)/dz;
    dG = - dG;
    for u = 1:4
        for v = 1:4
            for p = 1:4
                B1(u,v,p)=dG(v,p,u);
                B2(u,v,p)=dG(p,u,v);
                B3(u,v,p)=-dG(u,v,p);
            end
        end
    end
    B = B1+B2+B3;
    L12 = zeros(4,4,4);
    for o = 1:4
        for p = 1:4
            for u = 1:4
                for v = 1:4
                    L12(o,u,v)=L12(o,u,v)+(g(o,p)*B(u,v,p)/2);
                end
            end
        end
    end

end
