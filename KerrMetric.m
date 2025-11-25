function [g] = KerrMetric(m,a,x,y,z)

    eta = NPmetric();

    e0 = Frame(m,a,x,y,z);

    g = zeros(4,4);
    for u = 1:4
        for v = 1:4
            for a0 = 1:4
                for b = 1:4
                    g(u,v) = g(u,v) - e0(a0,u)*e0(b,v)*eta(a0,b);
                end
            end
        end
    end


    g = inv(g);


end
