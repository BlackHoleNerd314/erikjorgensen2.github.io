function [data] = MajoranaGamma(m,Psi,n)

  O2 = zeros(2,2);
  i = sqrt(-1);

  pauliX = [0,1;1,0];
  pauliY = [0,-i;i,0];
  pauliZ = [1,0;0,-1];

  Gamma0 = [O2,pauliY;pauliY,O2];
  Gamma1 = [i*pauliZ,O2;O2,i*pauliZ];
  Gamma2 = [O2,-pauliY;pauliY,O2];
  Gamma3 = [-i*pauliX,O2;O2,-i*pauliX];
  Gamma5 = -i*Gamma0*Gamma1*Gamma2*Gamma3;

  T = Gamma0*i;
  X = Gamma1*i;
  Y = Gamma2*i;
  Z = Gamma3*i;
  M = Gamma5*i;


  freq_mat = inv(T)*i;



  for spinor = 1:4
    Psi0(:,:,:,spinor) = fftn(Psi(:,:,:,spinor));
  end

  for px = (0:(n-1))/n;
    for py = (0:(n-1))/n;
      for pz = (0:(n-1))/n;
        px0 = px*n+1;
        py0 = py*n+1;
        pz0 = pz*n+1;
        Psi0(px0,py0,pz0,:) = Psi0(px0,py0,pz0,:) * expm(2*pi*(-i*inv(T)*M*m + px*inv(T)*X + py*inv(T)*Y + pz*inv(T)*Z));
      end
    end
  end
  %


  for spinor = 1:4
    Psi(:,:,:,spinor) = ifftn(Psi0(:,:,:,spinor));
  end

  data = Psi;


end
