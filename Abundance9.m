function [] = Abundance9(a)
% Temperature in MeV:
  a0 = a;
  eV = (2.72548/11604.51812)/a;
  MeV = eV/10^6;
% Nuclei Masses:
  Melectron0 = 5.48579909065 / 10^4;
  %Malpha = 7294.29954142;
  %Mproton = 1836.15267343;
  Mproton0 = 1.007276466621;
  Malpha0 = 4.001506179127;
  Mdeuteron0 = 2.013553212745;
  Mhelion0 = 3.014932247175;
  Mtriton0 = 3.01550071621;
  Mneutron0 = 1.00866491595;
% Particle Masses: 
  Mtop = 173210;
  Mhiggs = 125250;
  Mzero = 91187.6;
  Mweak = 80377;
  Mbottom = 4180;
  Mtauon = 1776.86;
  Mcharm = 1275;
  Mmuon = 105.6583755;
  Mstrange = 93.4;
  Mdown = 4.79;
  Mup = 2.01;
  Melectron = 0.51099895;
% Heavy Particle Masses:
  Mparticle1 = [Mhiggs];
  Mparticle2Lepton = [Mmuon,Mtauon];
  Mparticle2Quark = [Mup,Mcharm,Mtop,Mdown,Mstrange,Mbottom];
  Mparticle3 = [Mweak,Mzero,Mweak];
% Hadron Masses:
  Mhadron1 = [134.9766,139.57018,139.57018,493.677,493.677,497.614,497.614,497.614,497.614,547.862,957.78];
  Mhadron2 = [938.272,938.272,939.565,939.565,1115.683,1115.683,1189.37,1189.37,1192.642,1192.642,1197.449,1197.449,1314.86,1314.86,1321.71,1321.71];
  Mhadron3 = [775.11,775.11,775.26,782.65,891.66,891.66,895.81,895.81,1019.461];
  Mhadron4 = [1232,1232,1232,1232,1232,1232,1232,1232,1382.8,1382.8,1383.7,1383.7,1387.2,1387.2,1531.80,1531.80,1535.0,1535.0,1672.45,1672.45];
%%% Calculate Thermal History:
eps0 = 0.01;
eps1 = 10000;
eps2 = 0.01;
zeta = 0;
for iter0 = 1:1000000
  zeta = zeta+1/iter0^3;
end
% Fundamental Particles:
    Equark0 = 16;
    Ngluon0 = Equark0*1;
    Nparticle2Quark = 0*Mparticle2Quark;
    for n2 = 1:length(Mparticle2Quark);
        gpoints = 12;
        Mpoint = Mparticle2Quark(n2);
        z = Mpoint / MeV;
        du = eps2;
        u = z;
        E1 = 0;
        N1 = 0;
        for iter0 = 1:eps1;
            u = u + du;
            if u > eps1;
                E0 = 0;
            else;
                E0 = (15 * gpoints / pi ^ 4) * (u ^ 2 * sqrt(u ^ 2 - z ^ 2)) / (exp(u) + 1);
            end;
            N0 = (E0 / u);
            E1 = E1 + E0 * du;
            N1 = N1 + N0 * du;
        end;
        if z<eps0
          E1 = gpoints*7/8;
          N1 = gpoints*3/4;
        end
        Nparticle2Quark(n2) = N1;
        Equark0 = Equark0 + E1;
    end;
    Epoint = 0;
    Nparticle1 = 0*Mparticle1;
    for n1 = 1:length(Mparticle1);
        gpoints = 1;
        Mpoint = Mparticle1(n1);
        z = Mpoint / MeV;
        du = eps2;
        u = z;
        E1 = 0;
        N1 = 0;
        for iter0 = 1:eps1;
            u = u + du;
            if u > eps1;
                E0 = 0;
            else;
                E0 = (15 * gpoints / pi ^ 4) * (u ^ 2 * sqrt(u ^ 2 - z ^ 2)) / (exp(u) - 1);
            end;
            N0 = (E0 / u);
            E1 = E1 + E0 * du;
            N1 = N1 + N0 * du;
        end;
        if z<eps0
          E1 = gpoints*1;
          N1 = gpoints*1;
        end
        Nparticle1(n1) = N1;
        Epoint = Epoint + E1;
    end;
    Nparticle2Lepton = 0*Mparticle2Lepton;
    for n2 = 1:length(Mparticle2Lepton);
        gpoints = 4;
        Mpoint = Mparticle2Lepton(n2);
        z = Mpoint / MeV;
        du = eps2;
        u = z;
        E1 = 0;
        N1 = 0;
        for iter0 = 1:eps1;
            u = u + du;
            if u > eps1;
                E0 = 0;
            else;
                E0 = (15 * gpoints / pi ^ 4) * (u ^ 2 * sqrt(u ^ 2 - z ^ 2)) / (exp(u) + 1);
            end;
            N0 = (E0 / u);
            E1 = E1 + E0 * du;
            N1 = N1 + N0 * du;
        end;
        if z<eps0
          E1 = gpoints*7/8;
          N1 = gpoints*3/4;
        end
        Nparticle2Lepton(n2) = N1;
        Epoint = Epoint + E1;
    end;
    Nparticle3 = 0*Mparticle3;
    for n3 = 1:length(Mparticle3);
        gpoints = 3;
        Mpoint = Mparticle3(n3);
        z = Mpoint / MeV;
        du = eps2;
        u = z;
        E1 = 0;
        N1 = 0;
        for iter0 = 1:eps1;
            u = u + du;
            if u > eps1;
                E0 = 0;
            else;
                E0 = (15 * gpoints / pi ^ 4) * (u ^ 2 * sqrt(u ^ 2 - z ^ 2)) / (exp(u) - 1);
            end;
            N0 = (E0 / u);
            E1 = E1 + E0 * du;
            N1 = N1 + N0 * du;
        end;
        if z<eps0
          E1 = gpoints*1;
          N1 = gpoints*1;
        end
        Nparticle3(n3) = N1;
        Epoint = Epoint + E1;
    end;
    QuarkDensity0 = Equark0;
    ParticleDensity0 = Epoint;
% Composite Particles:
    Epoint = 0;
    N_1 = 0*Mhadron1;
    for n1 = 1:length(Mhadron1);
        gpoints = 1;
        Mpoint = Mhadron1(n1);
        z = Mpoint / MeV;
        du = eps2;
        u = z;
        E1 = 0;
        N1 = 0;
        for iter0 = 1:eps1;
            u = u + du;
            if u > eps1;
                E0 = 0;
            else;
                E0 = (15 * gpoints / pi ^ 4) * (u ^ 2 * sqrt(u ^ 2 - z ^ 2)) / (exp(u) - 1);
            end;
            N0 = (E0 / u);
            E1 = E1 + E0 * du;
            N1 = N1 + N0 * du;
        end;
        if z<eps0
          E1 = gpoints*1;
          N1 = gpoints*1;
        end
        N_1(n1) = N1;
        Epoint = Epoint + E1;
    end;
    N2 = 0*Mhadron2;
    for n2 = 1:length(Mhadron2);
        gpoints = 2;
        Mpoint = Mhadron2(n2);
        z = Mpoint / MeV;
        du = eps2;
        u = z;
        E1 = 0;
        N1 = 0;
        for iter0 = 1:eps1;
            u = u + du;
            if u > eps1;
                E0 = 0;
            else;
                E0 = (15 * gpoints / pi ^ 4) * (u ^ 2 * sqrt(u ^ 2 - z ^ 2)) / (exp(u) + 1);
            end;
            N0 = (E0 / u);
            E1 = E1 + E0 * du;
            N1 = N1 + N0 * du;
        end;
        if z<eps0
          E1 = gpoints*7/8;
          N1 = gpoints*3/4;
        end
        N2(n2) = N1;
        Epoint = Epoint + E1;
    end;
    N3 = 0*Mhadron3;
    for n3 = 1:length(Mhadron3);
        gpoints = 3;
        Mpoint = Mhadron3(n3);
        z = Mpoint / MeV;
        du = eps2;
        u = z;
        E1 = 0;
        N1 = 0;
        for iter0 = 1:eps1;
            u = u + du;
            if u > eps1;
                E0 = 0;
            else;
                E0 = (15 * gpoints / pi ^ 4) * (u ^ 2 * sqrt(u ^ 2 - z ^ 2)) / (exp(u) - 1);
            end;
            N0 = (E0 / u);
            E1 = E1 + E0 * du;
            N1 = N1 + N0 * du;
        end;
        if z<eps0
          E1 = gpoints*1;
          N1 = gpoints*1;
        end
        N3(n3) = N1;
        Epoint = Epoint + E1;
    end;
    N4 = 0*Mhadron4;
    for n4 = 1:length(Mhadron4);
        gpoints = 4;
        Mpoint = Mhadron4(n4);
        z = Mpoint / MeV;
        du = eps2;
        u = z;
        E1 = 0;
        N1 = 0;
        for iter0 = 1:eps1;
            u = u + du;
            if u > eps1;
                E0 = 0;
            else;
                E0 = (15 * gpoints / pi ^ 4) * (u ^ 2 * sqrt(u ^ 2 - z ^ 2)) / (exp(u) + 1);
            end;
            N0 = (E0 / u);
            E1 = E1 + E0 * du;
            N1 = N1 + N0 * du;
        end;
        if z<eps0
          E1 = gpoints*7/8;
          N1 = gpoints*3/4;
        end
        N4(n4) = N1;
        Epoint = Epoint + E1;
    end;
    N1 = N_1;
    HadronDensity0 = Epoint;
% Quark Gluon Plasma to Baryon Meson Gas Transition:
    if QuarkDensity0>HadronDensity0
        Density0 = ParticleDensity0 + HadronDensity0;
        NWeakons = [Nparticle3',Mparticle3'];
        NGluons = [0*Ngluon0,0];
        NQuarks = [0*Nparticle2Quark',Mparticle2Quark'];
        NLeptons = [Nparticle2Lepton',Mparticle2Lepton'];
        NHiggs = [Nparticle1',Mparticle1'];
        NHadrons4 = [N4',Mhadron4'];
        NHadrons3 = [N3',Mhadron3'];
        NHadrons2 = [N2',Mhadron2'];
        NHadrons1 = [N1',Mhadron1'];
        iscomposite = 1;
    else
        Density0 = ParticleDensity0 + QuarkDensity0;
        NHadrons4 = [0*N4',Mhadron4'];
        NHadrons3 = [0*N3',Mhadron3'];
        NHadrons2 = [0*N2',Mhadron2'];
        NHadrons1 = [0*N1',Mhadron1'];
        NWeakons = [Nparticle3',Mparticle3'];
        NGluons = [Ngluon0,0];
        NQuarks = [Nparticle2Quark',Mparticle2Quark'];
        NLeptons = [Nparticle2Lepton',Mparticle2Lepton'];
        NHiggs = [Nparticle1',Mparticle1'];
        iscomposite = 0;
    end
%% Input parameters of modern universe
  H0 = 0.9568;
  number = 6.1273089426/10^10; 
  photon = 1;
  neutrino = (5.25/2)*(4/11)^(4/3);   
  Vac = 0.6847;
  Mat = 0.2589+0.0486;
  RadPhoton = 0.0486*(2.72548/11604.51812)/(938.272041*10^6*number);
%% Energy Densities of Thermal Radiation
  gElectrons = 4;
  gPhotons = 2;
  z = Melectron/MeV;
  Nelectrons = 0;
  Eelectrons = 0;
  Pelectrons = 0;
  du = 0.01;
  u = z;
  for iter0 = 1:10000
    u = u + du;
    E0 = (15*gElectrons/pi^4)*(u^2*sqrt(u^2-z^2))/(exp(u)+1);
    N0 = (E0/u);
    Eelectrons = Eelectrons + E0*du;
    Nelectrons = Nelectrons + N0*du;
  end
  if z<eps0
    Eelectrons = 3.5;
    Nelectrons = 3;
  end
  Ephotons = gPhotons;
  T_neutrinos = MeV*((2+Eelectrons)/(2+3.5))^(1/3);
  Eneutrinos = 6*(7/8)*(T_neutrinos/MeV)^4;
  g0 = Eelectrons + Ephotons + Eneutrinos + Density0;
  Nneutrinos = 6*(3/4)*(T_neutrinos/MeV)^3;
  Nphotons = gPhotons;
%% Space Expansion and Time Evolution

  radiation0 = (g0/Ephotons)*RadPhoton;
  matter0 = Mat;
  vacuum0 = Vac;
  
  arm = radiation0/matter0;
  a2_term1 = (1-(1-(a/(2*arm)))*sqrt(1+(a/arm)));
  a2_term2 = (4*arm^2)/(3*sqrt(radiation0));
  t1 = (1/H0)*a2_term1*a2_term2;

  aml = (matter0/vacuum0)^(1/3);
  a3_term1 = log(sqrt(a/aml)^3+sqrt(1+(a/aml)^3));
  a3_term2 = 2/(3*sqrt(1-matter0));
  t2 = (1/H0)*a3_term1*a3_term2;

  if t1>t2
    t = t2;
  else
    t = t1;
  end
  
  arm = radiation0/matter0;
  a2_term1 = a^2/(2);
  a2_term2 = 1/(sqrt(radiation0));
  t3 = (1/H0)*a2_term1*a2_term2;

  if t<(1/10^15)
    t = t3;
  end

  s = t*(4.35456*10^17);

  s
  eV
  

  
  %%% Nuclear Fusion Abundances
  MeV = eV/10^6;
  secs = s;
  days = secs/86400;
  yrs = days/365.2425;
  keV = MeV*1000;
  n1_293=(-Melectron*(Mproton0-Mneutron0)/Melectron0);
  s0 = (1.1^2)*(n1_293-Melectron)^(-2);
  s1 = (1.47^2)*(78/1000)^(-2);
  Alpha1 = (0.248568317/4);
  Neutron1 = exp(-(s-s0)/880.2)*exp(-n1_293/sqrt(MeV^2+(n1_293-Melectron)^2))/(1+exp(-n1_293/sqrt(MeV^2+(n1_293-Melectron)^2)));
  n1 = (erf(8.5*log2(keV/78))+1)/2;
  Neutron2 = Neutron1*(n1);
  Alpha2 = Alpha1*(1-n1);
  Proton2 = 1-Alpha2*4-Neutron2;
  data=[Neutron2,Proton2,Alpha2];
  n2 = (erf(8.5*log2(keV/70))+1)/2;
  Deuteron_end = (1-n1)*26.4/10^6;
  Helion_end = (1-n1)*10.201/10^6;
  Triton_end = (1-n1)*0.99/10^7;
  Neutron_end = (1-n1)*1/10^6;
  Be7_end = (1-n1)*3.8/10^10;
  Li7_end = (1-n1)*2.9/10^11;
  Li6_end = (1-n1)*1.1/10^14;
  n3 = (-1.5*(n1-n2));
  Deuteron_mid = Deuteron_end + n3*4.5/10^3;
  Triton_mid = Triton_end + n3*0.83/10^4;
  Helion_mid = Helion_end - n3/10^6;
  Be7_mid = Be7_end * (1-n2);
  Li7_mid = Li7_end + n3*3.5/10^9;
  Li6_mid = Li6_end + 0.95*n3/10^12;
  n4 = 1 - (erf(8.5*log2(keV/280))+1)/2;
  n5 = n1*n4;
  n6 = n5/MeV^7;
  n7 = (n1)/MeV^10;
  Alpha_bump = 5*n6/10^13;
  Helion_bump = n6/10^16;
  Triton_bump = 0.7*(n6/MeV^2)/10^15;
  Deuteron_bump = (1.15627803/10^12) + n7/10^14.5;
  Neutron_sum = exp(-(s-s0)/880.2)*Neutron_end + Neutron2;
  Proton_sum = Proton2;
  Deuteron_sum = Deuteron_bump + Deuteron_mid;
  Triton_sum1 = Triton_bump + Triton_mid + (1.0005/3)/10^24;
  Helion_sum2 = Helion_bump + Helion_mid + (1/3)/10^24;
  Alpha_sum = Alpha_bump + Alpha2;
  Li6_sum = Li6_mid;
  Li7_sum2 = Li7_mid;
  Be7_sum1 = Be7_mid;
  Triton_decay = 0.5^(yrs/12.32);
  Be7_decay = 0.5^(days/53.217592);
  Triton_sum = Triton_sum1*Triton_decay;
  Helion_sum = Helion_sum2 + Triton_sum1*(1-Triton_decay);
  Be7_sum = Be7_sum1*Be7_decay;
  Li7_sum = Li7_sum2 + Be7_sum1*(1-Be7_decay);
  Nukes = zeros(5,5);
  Sum = 1*Neutron_sum+1*Proton_sum+2*Deuteron_sum+3*Triton_sum+3*Helion_sum+4*Alpha_sum+6*Li6_sum+7*Li7_sum+7*Be7_sum;
  Nukes(0+1,1+1) = Neutron_sum/Sum;
  Nukes(1+1,0+1) = Proton_sum/Sum;
  Nukes(1+1,1+1) = Deuteron_sum/Sum;
  Nukes(1+1,2+1) = Triton_sum/Sum;
  Nukes(2+1,1+1) = Helion_sum/Sum;
  Nukes(2+1,2+1) = Alpha_sum/Sum;
  Nukes(3+1,3+1) = Li6_sum/Sum;
  Nukes(3+1,4+1) = Li7_sum/Sum;
  Nukes(4+1,3+1) = Be7_sum/Sum;
  %disp(Nukes);
  Nukes0 = zeros(5,5);
  Li6_mass = 6;
  Li7_mass = 7;
  Be7_mass = 7;
  Nukes0(0+1,1+1) = Mneutron0;
  Nukes0(1+1,0+1) = Mproton0;
  Nukes0(1+1,1+1) = Mdeuteron0;
  Nukes0(1+1,2+1) = Mtriton0;
  Nukes0(2+1,1+1) = Mhelion0;
  Nukes0(2+1,2+1) = Malpha0;
  Nukes0(3+1,3+1) = Li6_mass;
  Nukes0(3+1,4+1) = Li7_mass;
  Nukes0(4+1,3+1) = Be7_mass;
  
  %%% Recombination Abundances
  
  z = eV*11604.51812/2.72548;
  
  z = z + 8.8;
  %%z = z0;
  
  
  p = Nukes(1+1,0+1);
  pn = Nukes(1+1,1+1);
  pnn = Nukes(1+1,2+1);
  ppn = Nukes(2+1,1+1);
  ppnn = Nukes(2+1,2+1);
  
  Hi = 1425;
  Hei = 2604;
  Heii = 6101;
  
  recomb = 1090.1;
  
  rate0 = (erf(z/10^4)/(z/10^4))*((z/recomb)*1.5)*(1090-1425)/log(578.6834137264191/10^7);
  
  ion0=(tanh((z-1425)/(rate0))+1)/2;
  ion1=(tanh((z-2604)/(rate0))+1)/2;
  ion2=(tanh((z-6101)/(rate0))+1)/2;
  p = p*ion0 + (578.6834137264191/10^7)*(1-ion0);
  
  pn_e = pn*(1-ion0);
  pn = pn*ion0;
  
  pnn_e = pnn*(1-ion0);
  pnn = pnn*ion0;
  
  ppn_e = ppn*(1-ion2);
  ppn = ppn*ion2;
  ppn_ee = ppn_e*(1-ion1);
  ppn_e = ppn_e*ion1;
  
  ppnn_e = ppnn*(1-ion2);
  ppnn = ppnn*ion2;
  ppnn_ee = ppnn_e*(1-ion1);
  ppnn_e = ppnn_e*ion1;
  
  
  
  n = Nukes(0+1,1+1);
  
  pppnnn = Nukes(3+1,3+1);
  pppnnnn = Nukes(3+1,4+1);
  ppppnnn = Nukes(4+1,3+1);
  
  ion3=(tanh((z-520.2)/(rate0))+1)/2;
  ion4=(tanh((z-7298.1)/(rate0))+1)/2;
  ion5=(tanh((z-11815.0)/(rate0))+1)/2;
  
  
  Li = 1;
  Li_e = Li*(1-ion5);
  Li = Li*ion5;
  Li_ee = Li_e*(1-ion4);
  Li_e = Li_e*ion4;
  Li_eee = Li_ee*(1-ion3)/2;
  Li_ee = Li_ee*(1-(1-ion3)/2);
  
  pppnnnn_e = pppnnnn*Li_e;
  pppnnn_e = pppnnn*Li_e;
  pppnnnn_ee = pppnnnn*Li_ee;
  pppnnn_ee = pppnnn*Li_ee;
  pppnnnn_eee = pppnnnn*Li_eee;
  pppnnn_eee = pppnnn*Li_eee;
  
  
  p_e = 1-(n+p+2*pn+2*pn_e+3*pnn+3*pnn_e+3*ppn+3*ppn_e+3*ppn_ee+4*ppnn+4*ppnn_e+4*ppnn_ee+7*pppnnnn+6*pppnnn+7*ppppnnn);
  
  pppnnnn = pppnnnn*Li;
  pppnnn = pppnnn*Li;
  
  
  
  %molecular abundances
  
  H2_0 = (15/10^7)*(1-erf((z-108)/50))/2;
  H2_1 = (2/10^7)*(1-erf((z-425)/20))/2;
  H2_2 = (1.3/10^12)*(1-erf((z-770)/30))/2;
  
  p_p_ee = H2_2+H2_1+H2_0;
  p_ee = H2_0/(2*10^4)*(z/100)^2;
  p_p_e = H2_1/(4*10^5)*(z/400)^3 + H2_0/(9*10^7)/(z/100);
  
  
  ppnn_p_ee = ((1-erf((z-320)/45))/2)/(z*(10^12));
  
  pn0 = pn + ((1-erf((z-1200)/275))/2)*(z/10^11)*((1+erf((z-120)/35))/2);
  
  pn_p_ee = ((1-erf((z-1000)/10))/2)*4/((3.7/10^10)/((3/10^12)/(z/300)^4)+1)/10^10;
  
  p_e = p_e-(2*p_p_ee+p_ee+2*p_p_e+ppnn_p_ee+pn_p_ee);
  
  pn_e = pn_e-(pn0-pn);
  pn = pn0;
  
  ppnn_ee = ppnn_ee-(ppnn_p_ee);
  
  pn_e = pn_e-(pn_p_ee);

  % molecular abundances rare
  recomb1 = (1-erf((z-1425)/194))/2;
  p_p_p_ee = (recomb1/10^20)+(p_p_ee*(z/10^11));
  bump12 = exp(-((z-1425)/194)^2/2);
  p_p_e = (p_p_e + bump12/10^18)*(p_e/0.75);
  p_ee = (p_ee + bump12/10^18)*(p_e/0.75);
  

  % isotopes abundances rare

  p_p_pn_ee = p_p_p_ee/(z*10^5);
  p_pn_e = p_p_e/10^6;
  pn_ee = p_ee*((4/3)*6*2.64)/10^15;
  
  % isotopes of common molecules
  
  recomb0 = (1-erf((z-recomb)/194))/2;





ppnn_pn_ee = recomb0*ppnn_p_ee*Nukes(1+1,1+1)/Nukes(1+1,0+1);
ppn_p_ee = recomb0*ppnn_p_ee*Nukes(2+1,1+1)/Nukes(2+1,2+1);





  
  % error fixing

  p_e = p_e * (1-erf((z-5000)/1000))/2;
  
  
  %%% normalizing nuclide abundances at chemical level



  p_sum0 = Nukes(1+1,0+1);
  pn_sum0 = Nukes(1+1,1+1);
  ppn_sum0 = Nukes(2+1,1+1);
  ppnn_sum0 = Nukes(2+1,2+1);

  p_sum1 = 2*p_p_ee+2*p_p_e+1*p_ee+3*p_p_p_ee+1*p_pn_e+2*p_p_pn_ee+1*pn_p_ee+1*ppnn_p_ee+1*ppn_p_ee;
  pn_sum1 = p_pn_e+p_p_pn_ee+pn_ee+1*pn_p_ee+1*ppnn_pn_ee;
  ppn_sum1 = ppn_p_ee;
  ppnn_sum1 = ppnn_p_ee+ppnn_pn_ee;

  p_sum2 = p+p_e;
  pn_sum2 = pn+pn_e;
  ppn_sum2 = ppn+ppn_e+ppn_ee;
  ppnn_sum2 = ppnn+ppnn_e+ppnn_ee;


  % balancing chemical reactions

  p_e = p_sum0-(p_sum1+p);
  pn_e = pn_sum0-(pn_sum1+pn);
  ppn_ee = ppn_sum0-(ppn_sum1+ppn+ppn_e);
  ppnn_ee = ppnn_sum0-(ppnn_sum1+ppnn+ppnn_e);
  
%% error fixing
a = 1/z;
if a<(1/10^9)
  p_e = 0;
  ppn_ee = 0;
  ppn_p_ee = 0;
end

  %%% collecting data


  
  
  
  data = iscomposite * [n,p,p_e,pn,pn_e,pnn,pnn_e,ppn,ppn_e,ppn_ee,ppnn,ppnn_e,ppnn_ee,...
pppnnn,pppnnn_e,pppnnn_ee,pppnnn_eee,pppnnnn,pppnnnn_e,pppnnnn_ee,pppnnnn_eee,ppppnnn,...
p_p_ee,p_p_e,p_ee,p_p_p_ee,p_pn_e,p_p_pn_ee,pn_ee,ppnn_p_ee,pn_p_ee,ppn_p_ee,ppnn_pn_ee];
  
  
  Q = [0,1,0,1,0,1,0,2,1,0,2,1,0,3,2,1,0,3,2,1,0,4,...
0,1,-1,1,1,1,-1,1,0,1,1];
  
  e0 = sum(data.*Q);
  

e1 = Nelectrons/(2*number) + e0;
E1 = Nelectrons/(2*number);
y1 = Nphotons/number;
v1 = Nneutrinos/number;

  data1 = [data,e1,E1,y1,v1];

  Q1 = [Q,-1,1,0,0];
  
  %% particle mass precision data

HI = 13.598434599702;
HeI = 79.005154539;
HeII = 54.4177655282;
LiI = 203.4861711;
LiII = 198.0944562;
LiIII = 122.45435913;
BeI = 399.14864;
BeII = 389.82595;
BeIII = 371.614789;
BeIV = 217.71858459;

LiBind0 = (LiI / (Melectron * 10^6)) * Melectron0;
BeBind0 = (BeI / (Melectron * 10^6)) * Melectron0;

Li7mass0 = 7.016003434 - 3*Melectron0 + LiBind0;
Li6mass0 = 6.0151228874 - 3*Melectron0 + LiBind0;
Be7mass0 = 7.01692871 - 4*Melectron0 + BeBind0;


  # Particle Masses
nM = (Mneutron0 * Melectron / Melectron0);
pM = (Mproton0 * Melectron / Melectron0);
p_eM = ((Mproton0 + Melectron0) * Melectron / Melectron0) - HI * (1 / 10 ^ 6);
pnM = (Mdeuteron0 * Melectron / Melectron0);
pn_eM = ((Mdeuteron0 + Melectron0) * Melectron / Melectron0) - HI * (1 / 10 ^ 6);
pnnM = (Mtriton0 * Melectron / Melectron0);
pnn_eM = ((Mtriton0 + Melectron0) * Melectron / Melectron0) - HI * (1 / 10 ^ 6);
ppnM = (Mhelion0 * Melectron / Melectron0);
ppn_eM = ((Mhelion0 + Melectron0) * Melectron / Melectron0) - HeII * (1 / 10 ^ 6);
ppn_eeM = ((Mhelion0 + 2 * Melectron0) * Melectron / Melectron0) - HeI * (1 / 10 ^ 6);
ppnnM = (Malpha0 * Melectron / Melectron0);
ppnn_eM = ((Malpha0 + Melectron0) * Melectron / Melectron0) - HeII * (1 / 10 ^ 6);
ppnn_eeM = ((Malpha0 + 2 * Melectron0) * Melectron / Melectron0) - HeI * (1 / 10 ^ 6);
pppnnnM = ((Li6mass0 + 0 * Melectron0) * Melectron / Melectron0) - 0 * (1 / 10 ^ 6);
pppnnn_eM = ((Li6mass0 + 1 * Melectron0) * Melectron / Melectron0) - LiIII * (1 / 10 ^ 6);
pppnnn_eeM = ((Li6mass0 + 2 * Melectron0) * Melectron / Melectron0) - LiII * (1 / 10 ^ 6);
pppnnn_eeeM = ((Li6mass0 + 3 * Melectron0) * Melectron / Melectron0) - LiI * (1 / 10 ^ 6);
pppnnnnM = ((Li7mass0 + 0 * Melectron0) * Melectron / Melectron0) - 0 * (1 / 10 ^ 6);
pppnnnn_eM = ((Li7mass0 + 1 * Melectron0) * Melectron / Melectron0) - LiIII * (1 / 10 ^ 6);
pppnnnn_eeM = ((Li7mass0 + 2 * Melectron0) * Melectron / Melectron0) - LiII * (1 / 10 ^ 6);
pppnnnn_eeeM = ((Li7mass0 + 3 * Melectron0) * Melectron / Melectron0) - LiI * (1 / 10 ^ 6);
ppppnnnM = ((Be7mass0) * Melectron / Melectron0);
p_p_eeM = (2 * p_eM) - (4.52 / 10 ^ 6);
p_p_eM = (-0.597 * 2 * HI / 10 ^ 6) + (2 * pM + Melectron);
p_eeM = (p_eM + Melectron) - (0.754195 / 10 ^ 6);
p_p_p_eeM = (3 * p_eM - Melectron) - (55.78 / 10 ^ 9);
p_pn_eM = p_p_eM - pM + pnM;
p_p_pn_eeM = p_p_p_eeM - pM + pnM;
pn_eeM = p_eeM - pM + pnM;
pn_p_eeM = p_p_eeM - pM + pnM;
ppnn_p_eeM = (ppnn_eeM + pM) - ((4.52 * 25.1) / (436 * 10 ^ 6));
ppn_p_eeM = ppnn_p_eeM - ppnnM + ppnM;
ppnn_pn_eeM = ppnn_p_eeM - pM + pnM;
eM = Melectron;

## Putting it all Together
dataM = [nM, pM, p_eM, pnM, pn_eM, pnnM, pnn_eM, ppnM, ppn_eM, ppn_eeM, ppnnM, ppnn_eM, ppnn_eeM, pppnnnM, pppnnn_eM, pppnnn_eeM,...
            pppnnn_eeeM, pppnnnnM, pppnnnn_eM, pppnnnn_eeM, pppnnnn_eeeM, ppppnnnM, p_p_eeM, p_p_eM, p_eeM, p_p_p_eeM, p_pn_eM,...
            p_p_pn_eeM, pn_eeM, ppnn_p_eeM, pn_p_eeM, ppn_p_eeM, ppnn_pn_eeM, eM, eM, 0 , 0];

  
matter1 = (data1.*dataM)*0.0486/938.272041;


names0 = ["n";"p";"p_e";"pn";"pn_e";"pnn";"pnn_e";"ppn";"ppn_e";"ppn_ee";"ppnn";"ppnn_e";"ppnn_ee";...
"pppnnn";"pppnnn_e";"pppnnn_ee";"pppnnn_eee";"pppnnnn";"pppnnnn_e";"pppnnnn_ee";"pppnnnn_eee";"ppppnnn";...
"p_p_ee";"p_p_e";"p_ee";"p_p_p_ee";"p_pn_e";"p_p_pn_ee";"pn_ee";"ppnn_p_ee";"pn_p_ee";"ppn_p_ee";"ppnn_pn_ee";...
%];%"E";"e";"y";"v"];
"e";"E";"y";"v"];

radiation1 = radiation0/a0;
Matter = matter1;#*radiation1;
PBHs = 0.2589;#*radiation1;
LambdaVac = Vac*a0^3;#*radiation1;

Total0 = radiation1+sum(Matter)+PBHs+LambdaVac;
Total = Total0/10^15;

Total1 = 0.2589/(41*42*43/6);

OtherFrac = int2str(round([LambdaVac,PBHs,radiation1]/Total1)');
names1 = ["LambdaVac";"PBHs";"RadTot"];

MatterFrac = int2str(round(Matter/Total1)');
#disp([names0,MatterFrac]);

MeV0 = (2.72548/11604.51812)/(10^6);

T_y = MeV/MeV0;
T_CMB = T_y;
T_v = T_neutrinos/MeV0;
T_qqqM = T_CMB * ((194/T_CMB)^(1/0.3)+1)^(-0.3);

o_v = (erf(log2(s/s0))+1)/2;
o_y = (erf(10*log2(a0*1090.1))+1)/2;


#disp([names1,OtherFrac]);

ans=[int2str([data1*10^15]'),names0,int2str([dataM/MeV0]')]
for ind3 = 1:max(size(NHadrons1))
  ans1 = [int2str([(NHadrons1(ind3,1)/number)*(10^15)]),"qq",int2str([Mhadron1(ind3)/MeV0]')]
end
for ind3 = 1:max(size(NHadrons2))
  ans2 = [int2str([(NHadrons2(ind3,1)/number)*(10^15)]),"qqq",int2str([Mhadron2(ind3)/MeV0]')]
end
for ind3 = 1:max(size(NHadrons3))
  ans3 = [int2str([(NHadrons3(ind3,1)/number)*(10^15)]),"qq",int2str([Mhadron3(ind3)/MeV0]')]
end
for ind3 = 1:max(size(NHadrons4))
  ans4 = [int2str([(NHadrons4(ind3,1)/number)*(10^15)]),"qqq",int2str([Mhadron4(ind3)/MeV0]')]
end
%
  name1 = ["higgs"];
  name2Lepton = ["muon";"tauon"];
  name2Quark = ["up";"charm";"top";"down";"strange";"bottom"];
  name3 = ["weak";"zero";"weak"];
ansHiggs = [int2str([(NHiggs(:,1)/number)*(10^15)]),name1,int2str([NHiggs(:,2)/MeV0])]
for ind4 = 1:2
  ansLepton = [int2str([(NLeptons(ind4,1)/number)*(10^15)]),name2Lepton(ind4,:),int2str([NLeptons(ind4,2)/MeV0])]
end
for ind4 = 1:6
  ansQuark = [int2str([(NQuarks(ind4,1)/number)*(10^15)]),name2Quark(ind4,:),int2str([NQuarks(ind4,2)/MeV0])]
end
%
ansGluon = [int2str([(NGluons(:,1)/number)*(10^15)]),["gluons"],int2str([NGluons(:,2)/MeV0])]
for ind4 = 1:3
  ansWeak = [int2str([(NWeakons(ind4,1)/number)*(10^15)]),name3(ind4,:),int2str([NWeakons(ind4,2)/MeV0])]
end
%
o_y
o_v
T_qqqM
T_y
T_v




###########################
IGM = [0.59,4,5,-6,-4];
WHIM = [0.24,5,7,-6,-4];
ICM = [0.04,7,8,-3,-2];
CGM = [0.05,4,6,-4,-2];
ISM_MM = [0.0013,1,1.7,2,4];
ISM_CNM = [0.0030,1.7,2,1.3,1.7];
ISM_WNM = [0.0038,3.8,4,-0.7,-0.3];
ISM_WIM = [0.0014,4.9,4.9,-1,-0.5];
ISM_HIM = [0.0004,6,7,-2.5,-2];
Stars = [0.07,6,8,24,26];
WhiteDwarf = [0.05,0.6];
NeutronStar = [0.005,1.4];
StellarBlackHole = [0.003,20];
GalacticBlackHole = [0.0001,4*10^6];
MainSequenceStars = [40,30,15,5,2;0.08,0.5,1.5,8,25];


nuke_pn = zeros(137,137+256);
nuke_II = zeros(137,137+256);
nuke_1a = zeros(137,137+256);

nuke_pn(6,12) = 69;
nuke_pn(7,14) = 13;
nuke_1a(26,56) = 3;
nuke_1a(14,28) = 2;
nuke_1a(16,32) = 1;
nuke_II(8,16) = 156;
nuke_II(6,12) = 23;
nuke_II(7,14) = 4;
nuke_II(10,20) = 17;
nuke_II(12,24) = 6;
nuke_II(14,28) = 4;
nuke_II(16,32) = 2;
nuke_II(26,56) = 2;





end















































































































































































































