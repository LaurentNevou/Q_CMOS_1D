function[OUTPUT,ro3DEfn,ro3DEfp,ErrVec] = Poisson_f(Structure,En,Ep,T,EfL,EfR,Nloops,tau0,NNewton,k,Guess,Video_convergence)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h     = 6.62606896E-34;           %% Planck constant J.s
hbar  = h/(2*pi);
e     = 1.602176487E-19;          %% charge de l electron Coulomb
m0    = 9.10938188E-31;           %% electron mass kg
%c     = 2.99792458e8;             %% speed of light (m/s)
Epsi0 = 8.854187817620E-12;       %% constant dielectric du vide F/m
%mu0   = 1/(Epsi0*c^2);            %% permeabiliy du vide
kB    = 1.3806488E-23;            %% Boltzmann's constant (J/K)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z     = Structure(:,1)';
V0    = Structure(:,2)';
Eg    = Structure(:,3)';
Dop   = Structure(:,4)';
Epsi  = Structure(:,5)';
EfXX  = Structure(:,6)';
Mass_n= Structure(:,7)';
Mass_p= Structure(:,8)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Meshgrid of density matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ro3Dn_const = (1/(2*pi^2)) * ( (2*e*Mass_n*m0/(hbar^2)).^(3/2) );
ro3Dp_const = (1/(2*pi^2)) * ( (2*e*Mass_p*m0/(hbar^2)).^(3/2) );

[ro3Dn_const_M,EEn]=meshgrid(ro3Dn_const,En); % put the vector Mass_n in a matrix En-long
[ro3Dp_const_M,EEp]=meshgrid(ro3Dp_const,Ep); % put the vector Mass_p in a matrix Ep-long

ro3Dn = ro3Dn_const_M .* sqrt(  EEn );
ro3Dp = ro3Dp_const_M .* sqrt( -EEp );

[Eg_M]=meshgrid(Eg,Ep);       % put the vector Gap in a matrix E-long

[EfXX_Mn,EEn]=meshgrid(EfXX,En);
[EfXX_Mp,EEp]=meshgrid(EfXX,Ep);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Building operator matrix for Newton-Raphson Algorithm %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HH=zeros(length(z),length(z));

for ii=2:length(z)-1,
    
    dxb = z(ii)-z(ii-1);                % backward difference   (i)-(i-i)
    dxf = z(ii+1)-z(ii);                % forward difference    (i+1)-(i)
    
    HH(ii,ii-1) =  ( Epsi(ii) + Epsi(ii-1) ) / ( dxb*(dxb+dxf) );
    HH(ii,ii+1) =  ( Epsi(ii) + Epsi(ii+1) ) / ( dxf*(dxb+dxf) );
    HH(ii,ii)   = -(HH(ii,ii-1)+HH(ii,ii+1)); %% good and working
end

% Boundary conditions
HH(1,1)      =  HH(2,2);
HH(1,2)      =  HH(2,1);
HH(end,end)  =  HH(end-1,end-1);
HH(end,end-1)=  HH(end-1,end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Voltage=EfXX(end)-EfXX(1);
Vbi=(EfL-EfR);
Fbi=(EfL-EfR)/(z(end)-z(1));
Vs1 = -(EfR-EfL-Voltage)/(z(end)-z(1))*z;
Vs2 = Guess-V0;

ntot=0;
dVV=1e-5;
minErr = 1e-10;     % minimum error on the potential at which the program stop
ErrVec=[];
sumVbitotVec=1;
nloop=1;

if Video_convergence==1
    figure('position',[100 100 1000 800]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Start of the Poisson s loop %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while (nloop<Nloops)
         
  nloop;

  Vbending=Vs1; 
  if (nloop>NNewton) || k>1
     Vbending=Vs2; 
  end

  Vbitot=V0+Vbending;
  tau = tau0*(1 + 2^((nloop - Nloops*0.8 )/10)); %the tau will increase at each loop

  %%%%%%%%%%%%%%%%%%%%%%%%%%% matrix density calcul %%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Vbitot_Mn=meshgrid(Vbitot,En);
  Vbitot_Mp=meshgrid(Vbitot,Ep);

  %%%%%%%%%%%%%%%%%%% calcul of the electrons density %%%%%%%%%%%%%%%%%%%%%%%%%%

  FEc = 1./(1+exp((EEn +Vbitot_Mn -EfXX_Mn)/(kB*T/e))) ; 
  ro3DEfn = ro3Dn .* FEc  ;
  NtotX = trapz(En,ro3DEfn);

  %%%%%%%%%%%%%%%%%%%%%% calcul of the holes density %%%%%%%%%%%%%%%%%%%%%%%%%%%

  FEv = 1./(1+exp(-( EEp-Eg_M +Vbitot_Mp -EfXX_Mp )/(kB*T/e))) ; 
  ro3DEfp=ro3Dp.*FEv;   
  PtotX=trapz(Ep,ro3DEfp);
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  NPtotX=NtotX-PtotX-Dop;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
  if (nloop>NNewton-1) || k>1 
    % this is needed only for the Newton-Raphson loop
    dVbitot_Mn=meshgrid(Vbitot+dVV,En);
    dVbitot_Mp=meshgrid(Vbitot+dVV,Ep);

    dFEc = 1./(1+exp((EEn +dVbitot_Mn -EfXX_Mn)/(kB*T/e))) ; 
    dro3DEfn = ro3Dn .* dFEc  ;
    dNtotX = trapz(En,dro3DEfn);

    dFEv = 1./(1+exp(-( EEp-Eg_M +dVbitot_Mp -EfXX_Mp )/(kB*T/e))) ; 
    dro3DEfp=ro3Dp.*dFEv;   
    dPtotX=trapz(Ep,dro3DEfp);

    dNPtotX=(dNtotX-dPtotX-Dop)-NPtotX;
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
  if (nloop<NNewton) && k==1    % => Damping injection method
       
    ntot = ntot + (NPtotX-ntot)/tau;    %% It add slowly the total number of electrons in order to converge

    %%%%%%%%%%%%%%%%%%%%% Electrical Field calculation %%%%%%%%%%%%%%%%%%%%%%%%%

    F   = e*cumtrapz(z,ntot)./(Epsi0*Epsi);     % integal on a nonlinear grid z
    MF  = trapz(z,F)/(z(end)-z(1));  % MF is the mean(F) function on a nonlinear grid z
    F   = F - MF - Fbi  - Voltage/(z(end)-z(1)) ;
          
    %%%%%%%%%%%%%%%%%%%%% New potentiel calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Vsold=  Vs1;                % storing the old value
    Vs1  = -cumtrapz(z,F);      % integal on a nonlinear grid z
    Vs   =  Vs1;                % storing in order to compute the error at the next iteration
       
  elseif (nloop>NNewton-1) || k>1   % => Newton Raphson algorithm
    
    ntot=NPtotX;
    Vbegin=V0(1);
    Vend=V0(end)+Vbi+Voltage;
    
    NPtotX(end) = NPtotX(end) +  Vend  * Epsi(end) * Epsi0/(e* (z(end) - z(end-1))^2) ;
    NPtotX(1)   = NPtotX(1)   +  Vbegin * Epsi(1)   * Epsi0/(e* (z(2)   - z(1)    )^2) ;
    %NPtotX(end) =   Vend   * Epsi(end) * Epsi0/(e* (z(end) - z(end-1))^2) ;
    %NPtotX(1)   =   Vbegin * Epsi(1)   * Epsi0/(e* (z(2)   - z(1)    )^2) ;
    
    FFF = -(HH*Vbitot')*Epsi0/e - NPtotX';  %% FFF is the fonctionnel that should converge to zero
    Jacobian = -HH*Epsi0/e - diag(dNPtotX/dVV);

    %Vs2 = Vbitot' - inv(Jacobian)*FFF;
    Vs2 = Vbitot' - Jacobian\FFF;
    Vs2=Vs2';
    Vs=Vs2;            % storing in order to compute the error at the next iteration

    F=-e*cumtrapz(z,ntot)./Epsi/Epsi0; %% Calcul du champs just for plotting
     
  end
   
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
  if Video_convergence==1
    cla
    hold on
    plot(z*1e9,Vbitot   ,'b-')
    plot(z*1e9,Vbitot-Eg,'b-')
    plot(z*1e9,EfXX     ,'m-')
    xlim([0 z(end)*1e9])
    ylim([min([Vbitot-Eg Voltage])-0.5 max([Vbitot Voltage])+0.5])
    title(strcat('nloop=',num2str(nloop)))
    %xlabel('z (nm)')
    %ylabel('Energy (eV)')
    pause(0.01)
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%% for the plotting of the graph of the convergence %%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Err=abs(  1 - sumVbitotVec(end)/sum(Vs)  );
  sumVbitotVec(k) = sum(Vs);
  ErrVec = [ErrVec Err];
  nloop=nloop+1;

  if Err<minErr
     Err;
     break 
  end

  OUTPUT = [Vbitot' F' NtotX' PtotX' ntot'];

end
