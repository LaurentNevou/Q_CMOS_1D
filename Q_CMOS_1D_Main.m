%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% last update 17January2021, lne %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code solves the Poisson equation in 1D CMOS transistor.
% As a results, it gives the band bending profile for any applied voltage.
% The program computes also the typical CV curve of a MOS transistor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If the code doesn t converge:
% -> increase the damping, tau0
% -> increase the amount of loops, Nloops. Nloops should be more than 3 times higher than tau0
% -> The Newton-Raphson algo should NOT start too early because the "guess" won t be good enough
% -> increase the amount of points
% -> increase the resolution dE
% -> increase the temperature (T=0K is very bad while T=10K is already much better)
% -> decrease the doping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ref: Newton-Raphson algorithm
%
% "Newton-Raphson solution of Poisson's equation in a pn diode"
% R. A. Jabr, M. Hamad and Y. M. Mohanna
% International Journal of Electrical Engineering Education
%  Volume: 44 issue: 1, page(s): 23-33, Issue published: January 1, 2007
% https://doi.org/10.7227/IJEEE.44.1.3
% https://journals.sagepub.com/doi/10.7227/IJEEE.44.1.3
%
% Sun Hee Lee
% "DEVELOPMENT OF A MASSIVELY PARALLEL NANOELECTRONIC MODELING TOOL AND ITS APPLICATION TO QUANTUM COMPUTING DEVICES"
% Purdue University, 2011
% APPENDICE-B. ITERATIVE METHOD FOR SOLVING POISSON EQUATION, page 104
% https://engineering.purdue.edu/gekcogrp/publications/theses/PhD_11_2011_Sunhee_Lee_PhD_Thesis_main.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h     = 6.62606896E-34;           %% Planck constant J.s
hbar  = h/(2*pi);
e     = 1.602176487E-19;          %% charge de l electron Coulomb
m0    = 9.10938188E-31;           %% electron mass kg
c     = 2.99792458e8;             %% speed of light (m/s)
Epsi0 = 8.854187817620E-12;       %% constant dielectric du vide F/m
mu0   = 1/(Epsi0*c^2);            %% permeabiliy du vide
kB    = 1.3806488E-23;            %% Boltzmann's constant (J/K)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Convergence parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NNewton = 250;      % number of the loop at which starts the Newton Raphson algorithm
Nloops  = 1500;     % number of loops at which it stops
tau0    = 800;      % Damping value for the guessed solution

T=300;              % Temperature in Kelvin (never put zero)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Turm on Graph and Saving %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0 for turn off
% 1 for turn on
CBand1D=1;    Xfig=10;Yfig=100;Wfig=1000;Hfig=800;
CBand3D=1;

Convergence=1;
Video_convergence=0;
Video_Voltage=0;
TotCharge_density=0;
Charge_density=1;
band3D=0;

Capa_graph=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Voltage sweep %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In order to compute the Capacitance C=dQ/dV, we need dV
% To make it easier, dV must be constant (homogeneous voltage grid)
dV=0.02;
%Voltage=2;
Voltage=-5:dV:5;
%Voltage=5:-dV:-5;
%Voltage=-2:dV:2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Library;                  % load material parameter DB from "materialDB_ZB.csv"
ExtractParameters;        % extract parameter from the Library

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% import the layer structure file %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

input_CMOS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Energy grid definition: the grid is moving respect to the bending %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Electron Energy grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

En1 = linspace( 0 , 0.35, 51 );
En2 = linspace( En1(end)+0.01 , 1, 10 );
En  = [En1 En2]; En = sort(En);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Holes Energy grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ep1 = linspace( 0 , -0.15, 20 );
Ep2 = linspace( Ep1(end)-0.01 , -1, 10 );
Ep=[Ep1 Ep2]; Ep = sort(Ep);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% NOTHING TO CHANGE ANYMORE !!! %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Grabbing the parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zt   = M(:,end-3)*1e-9;    % conversion of the length from Angstrom to meter

%%%%%%%%%%%%%%%%%%%%%%% Eg = Eg0 - (a*T.^2)./(T + b) %%%%%%%%%%%%%%%%%%%%%%%%%%%
EgG  = M(:,idx_Eg6c) - (M(:,idx_alphaG)*T^2) ./ (T+M(:,idx_betaG));   % Bandgap at Gamma point
EgX  = M(:,idx_EgX)  - (M(:,idx_alphaX)*T^2) ./ (T+M(:,idx_betaX));   % Bandgap at X point
EgL  = M(:,idx_EgL)  - (M(:,idx_alphaL)*T^2) ./ (T+M(:,idx_betaL));   % Bandgap at L point

Egt=min([EgG EgX EgL],[],2);

VBOt = M(:,idx_VBO);
CBOt = Egt+VBOt;         % CBO form band gap difference and temperature
Epsit= M(:,idx_Epsi);   %(used for Poisson solver only)

Doptn=M(:,end-2)*1e18*1e6;  % n doping conversion from cm-3 to m-3
Doptp=M(:,end-1)*1e18*1e6;  % p doping conversion from cm-3 to m-3
Dopt=Doptn-Doptp;

Masstn = M(:,idx_me);
Masstp = M(:,idx_mhh);
%Masstp = ( M(:,idx_mhh).^(3/2) + M(:,idx_mlh).^(3/2) ).^(2/3) ;

Pt=M(:,end);

Ntott=Dopt.*zt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Discretisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here, I descretize the grid z and a lot of other parameters

z(1)=0;  V0(1)=CBOt(1); Dop(1)=Dopt(1); Mass_n(1)=Masstn(1); Mass_p(1)= Masstp(1);Eg(1)=Egt(1);Epsi(1)=Epsit(1);
dzz=1E-12;

for i=1:length(zt)
    t=zt(i);
    zv     = linspace( z(end)+dzz , z(end) + t , Pt(i) );
    z      = [ z                        zv        ];
    V0     = [ V0      ones(size(zv)) * CBOt(i)   ];
    Dop    = [ Dop     ones(size(zv)) * Dopt(i)   ];
    Mass_n = [ Mass_n  ones(size(zv)) * Masstn(i) ];
    Mass_p = [ Mass_p  ones(size(zv)) * Masstp(i) ];
    Eg     = [ Eg      ones(size(zv)) * Egt(i)    ];
    Epsi   = [ Epsi    ones(size(zv)) * Epsit(i)  ];
end

[ZZn,EEn]=meshgrid(z,En);
[ZZp,EEp]=meshgrid(z,Ep);
V0=V0-V0(1); % It seems without this offset, the NewtonRaphson algo does not work...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Finding the boundary conditions %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EfL = Fermi3D_np( Eg(1),T,Masstn(1),Masstp(1),Dop(1) );          % Fermi level on the left side
EfR = Fermi3D_np( Eg(end),T,Masstn(end),Masstp(end),Dop(end) );  % Fermi level on the right side

EfL = EfL+V0(1);
EfR = EfR+V0(end);

distance_L=0;
for i=1:Fermi_layerbreak_L
    distance_L=distance_L + zt(i);
end

distance_R=0;
for i=1:Fermi_layerbreak_R
    distance_R=distance_R + zt(i);
end

idxg = find( abs( z-distance_L ) <1e-10 , 1); % Here the index in the vector z where the Fermi level will be broken
idxd = find( abs( z-distance_R ) <1e-10 , 1); % Here the index in the vector z where the Fermi level will be broken

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% loop over the Voltage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Video_Voltage==1
    figure('position',[Xfig Yfig Wfig Hfig],'color','w')
    subplot(1,1,1,'fontsize',20)
end

ErrVec=[];

for k=1:length(Voltage)     % big loop over the applied voltage
  
  display(sprintf('Simulation at : %.2f Volt' ,Voltage(k)))
  
  %%%%%%%%%%%%%%%%%%% building of the quasi Fermi level %%%%%%%%%%%%%%%%%%%%%%%%
    
  for i=1:idxg
      EfXX(i)=EfL;
  end
  for i=idxg+1:idxd-1
      EfXX(i)= Voltage(k)*z(i)/(z(idxd)-z(idxg)) + EfL - Voltage(k)*z(idxg)/(z(idxd)-z(idxg));
  end
  for i=idxd:length(z)
      EfXX(i)=EfL+Voltage(k);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  Structure=[z' V0' Eg' Dop' Epsi' EfXX' Mass_n' Mass_p'];

  if k==1
    Guess=z*0;
  else
    Guess=Vbitot;
    Video_convergence=0;
  end

  [OUTPUT,ro3DEfn,ro3DEfp,Err] = Poisson_f(Structure,En,Ep,T,EfL,EfR,Nloops,tau0,NNewton,k,Guess,Video_convergence);
  
  Vbitot = OUTPUT(:,1)';
  F      = OUTPUT(:,2)';
  NtotX  = OUTPUT(:,3)';
  PtotX  = OUTPUT(:,4)';
  ntot   = OUTPUT(:,5)';
  ErrVec = [ErrVec Err];
    
  if Video_Voltage==1
    cla
    hold on
    plot(z*1e9,Vbitot,'b-')
    plot(z*1e9,Vbitot-Eg,'b-')
    plot(z*1e9,EfXX,'m.-','linewidth',1)
    hold off
    xlim([0 z(end)*1e9])
    ylim([min(Vbitot-Eg )-0.5 max(Vbitot)+0.5])
    title(strcat('Voltage=',num2str(Voltage(k)),'V'))
    %xlabel('z (nm)')
    %ylabel('Energy (eV)')
    pause(0.01)
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % storing the total number of charge in a matrix in order to plot the capacitance

  Qtot_n(k,:) = NtotX ;
  Qtot_p(k,:) = PtotX ;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

if CBand1D==1 || CBand3D==1;
  figure('position',[Xfig Yfig Wfig Hfig],'color','w')

  xscale =[z(1) z(end)]*1e9;
  yscale=[min(Vbitot-Eg)-0.5 max(Vbitot)+0.5];

  subplot(1,1,1,'fontsize',20)
  hold on; grid on;box on;
  colormap(jet)
  col=colormap;

  if CBand3D==1
      
    [Vbitot_Mn]=meshgrid(Vbitot,En);
    [Vbitot_Mp]=meshgrid(Vbitot,Ep);
    [Eg_M]=meshgrid(Eg,Ep);
    grid off
    RRRn=log10(ro3DEfn*1e-6); RRRn(RRRn==inf)=0; RRRn(RRRn==-inf)=0;
    RRRp=log10(ro3DEfp*1e-6); RRRp(RRRp==inf)=0; RRRp(RRRp==-inf)=0;
    
    pcolor(ZZn*1e9,EEn+Vbitot_Mn,RRRn)
    pcolor(ZZp*1e9,EEp+Vbitot_Mp-Eg_M,RRRp)
    plot(z*1e9,Vbitot,'w-','linewidth',2)
    plot(z*1e9,Vbitot-Eg,'w-','linewidth',2)
    %plot(z*1e9,V0,'b--','linewidth',1)
    
    set(gca,'color',col(1,:))
    shading flat
    caxis([16 20])
    hcb=colorbar;
    title(hcb,'\fontsize{8}log10[cm-3]')

  else
    plot(z*1e9,Vbitot   ,'b-' ,'linewidth',2)
    plot(z*1e9,Vbitot-Eg,'b-' ,'linewidth',2)
    plot(z*1e9,V0,'k--','linewidth',1)
    plot(z*1e9,V0-Eg,'k--','linewidth',1)
  end

  plot(z*1e9,EfXX,'g-','linewidth',1)
  text(z(1)*1e9,EfXX(1)-0.15,'\color{green}Fermi')
  text(z(end)*1e9*0.95,EfXX(end)+0.15,'\color{green}Fermi')

  xlim(xscale)
  ylim(yscale)
  xlabel('z (nm)')
  ylabel('Energy (eV)')
  title(strcat('\fontsize{15}T=',num2str(T),'K; Voltage=',num2str(Voltage(end)),'V'))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Convergence==1;
    
    figure
    semilogy(ErrVec,'b.-')  % convergence graph
    grid on;
    xlim([0 length(ErrVec)])
    xlabel('Convergence cycle');
    ylabel('Oscillation of the Build-in potential (%)');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if TotCharge_density==1;
    
    figure('position',[410 50 400 400]);

    subplot(3,1,1)
    plotyy(z*1e9,ntot*1e-6,z*1e9,Vbitot)
    xlabel('z (nm)');
    ylabel('Charge density (cm-3)');
    legend('ntot','Vbitot')
    xlim([z(1) z(end)]*1e9)
    
    subplot(3,1,2)
    plot(z*1e9,(NtotX-PtotX-Dop)*1e-6,'r')
    hold on
    plot(z*1e9,ntot*1e-6,'b')
    plot(z*1e9,ntot*1e-6 + Dop*1e-6,'g')
    xlabel('z (nm)');
    ylabel('Charge density (cm-3)');
    legend('NPtotX','ntot','ntot+Dop')
    xlim([z(1) z(end)]*1e9)
    
    subplot(3,1,3)
    hold on
    plot(z*1e9,F,'g')
    xlabel('z (nm)');
    ylabel('Electric field (V/m)');
    legend('F')
    xlim([z(1) z(end)]*1e9)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Charge_density==1;

    figure('position',[520 50 400 400]);

    subplot(2,1,1)
    yscale=[1e0 1e20];
    ylim(yscale)
    %hold on
    semilogy(z*1e9,NtotX*1e-6,'r')
    hold on;grid on;box on;
    semilogy(z*1e9,PtotX*1e-6,'g')
    xlabel('z (nm)');
    ylabel('Charge density (cm-3)');
    xlim([z(1) z(end)]*1e9)
    ylim(yscale)
    
    subplot(2,1,2)
    plot(z*1e9,NtotX*1e-6,'r')
    hold on;grid on;box on;
    plot(z*1e9,PtotX*1e-6,'g')
    xlabel('z (nm)');
    ylabel('Charge density (cm-3)');
    xlim([z(1) z(end)]*1e9)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if band3D==1

    [Vbitot_Mn]=meshgrid(Vbitot,En);
    [Vbitot_Mp]=meshgrid(Vbitot,Ep);
    [Eg_M]=meshgrid(Eg,Ep);
    
    figure('position',[1230 50 600 500]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    subplot(2,1,1)
    hold on
    pcolor(ZZn*1e9,EEn+Vbitot_Mn,(ro3DEfn))
    plot(z*1e9,EfXX,'g.-','linewidth',1)
    plot(z*1e9,Vbitot,'r.-')
    shading flat
    colormap(jet)
    xlabel('z (nm)')
    ylabel('Energy (eV)')
    %zlabel('density of states * Fermi distribution')
    xlim([z(1) z(end)]*1e9)
    
    subplot(2,1,2)
    hold on
    pcolor(ZZp*1e9,EEp+Vbitot_Mp-Eg_M,(ro3DEfp))
    plot(z*1e9,EfXX,'g.-','linewidth',1)
    plot(z*1e9,Vbitot-Eg,'r.-')
    shading flat
    colormap(jet)
    xlabel('z (nm)')
    ylabel('Energy (eV)')
    %zlabel('density of states * Fermi distribution')
    xlim([z(1) z(end)]*1e9)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Computes the capa from book formula %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% which book?
Vbi = EfL-EfR;
Cj0 = sqrt(e*Epsi0*mean(Epsit)/(2*abs(Vbi)) * abs(  Dopt(1)*Dopt(end) / ( abs(Dopt(1))+abs(Dopt(end)) ) ) ); %% (F/m2)
Cj  = Cj0./sqrt(1-(+Voltage/abs(Vbi)) );  %% (F/m2) for NMOS
%Cj  = Cj0./sqrt(1-(-Voltage/abs(Vbi)) );  %% (F/m2) for PMOS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( Capa_graph==1 ) && ( length(Qtot_n(:,1)) > 2 )

    figure('position',[830 50 400 400]);

    subplot(3,1,1)
    hold on
    pcolor(z*1e9,Voltage,Qtot_n*1e-20)
    shading flat
    colormap(jet)
    xlabel('z (nm)');
    ylabel('Voltage (Volt)');
    title('electrons charge')

    subplot(3,1,2)
    hold on
    pcolor(z*1e9,Voltage,Qtot_p*1e-20)
    shading flat
    colormap(jet)
    xlabel('z (nm)');
    ylabel('Voltage (Volt)');
    title('holes charge')

    subplot(3,1,3)
    hold on
    pcolor( z*1e9 , Voltage , (Qtot_n-Qtot_p)*1e-20 )
    shading flat
    colormap(jet)
    xlabel('z (nm)');
    ylabel('Voltage (Volt)');
    title('Total charge')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    capa_dist_L=0;
    for i=1:capa_layerbreak_L
        capa_dist_L=capa_dist_L + zt(i);
    end
    capa_dist_R=0;
    for i=1:capa_layerbreak_R
        capa_dist_R=capa_dist_R + zt(i);
    end
    capa_L = find( abs(z-capa_dist_L)<1e-10 , 1 );
    capa_R = find( abs( z - capa_dist_R ) <1e-10  , 1 );
 
    % here, we cut the map of the charge because the charges must be
    % conserved, therefore, we should keep either the gate charges, or the
    % substrate charges
    
    Qtot_n=trapz( z(capa_L:capa_R) , Qtot_n(:,capa_L:capa_R) , 2 );
    Qtot_p=trapz( z(capa_L:capa_R) , Qtot_p(:,capa_L:capa_R) , 2 );
    
    Cap_n=e*diff( Qtot_n , 1 , 1 ) / dV;
    Cap_p=e*diff( Qtot_p , 1 , 1 ) / dV;
    Cap=Cap_n-Cap_p;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure('position',[1200 50 700 600],'color','w');
    
    subplot(2,1,1,'fontsize',15)
    hold on; grid on; box on;
    plot(Voltage,Qtot_n,'r.-') 
    plot(Voltage,Qtot_p,'g.-') 
    plot(Voltage,Qtot_n-Qtot_p,'b.-') 
   
    xlabel('Voltage (Volt)');
    ylabel('Charge (m-2)');
    legend('N charges','N charges','N-P charges')
    
    subplot(2,1,2,'fontsize',15)
    hold on; grid on; box on;
    plot(Voltage(2:end),abs(Cap),'b-','linewidth',2)
    plot(-Voltage,real(Cj),'m--')
    plot(Voltage(2:end),abs(Cap_n),'r-')
    plot(Voltage(2:end),abs(Cap_p),'g-')
    
    xlabel('Voltage (Volt)');
    ylabel('Capacitance (F/m2 or pF/um2)');
    legend('Poisson','Formula')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%