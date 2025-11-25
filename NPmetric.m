function [g] = NPmetric()

    %g = diag([ones(1,3),-1]);
    g = zeros(4,4);
    g(1,2) = 1;
    g(2,1) = 1;
    g(3,4) = -1;
    g(4,3) = -1;


end
