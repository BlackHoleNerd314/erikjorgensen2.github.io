function [data] = Nbody1(M,A,X,V,m,a,x,v,n0)
  N0 = 1;
  for N = 1:n0
    % Orbit m around M:
    data1 = ReferenceBody(M,X,V,A,x,v,a);
    n1 = length(data1);
    x = data1(n1-1,:);
    v = data1(n1,:)-data1(n1-2,:);
    % collect data:
    for N1 = 1:n1
      data0(N1+N0,:) = data1(N1,:);
    end
    N0 = N0 + n1;
    % Orbit M around m:
    data2 = ReferenceBody(m,x,v,a,X,V,A);
    n2 = length(data2);
    X = data2(n2-1,:);
    V = data2(n2,:)-data2(n2-2,:);
    % collect data:
    for N2 = 1:n2
      data0(N2+N0,:) = data2(N2,:);
    end
    N0 = N0 + n2;
  end
  data = data0;
end
