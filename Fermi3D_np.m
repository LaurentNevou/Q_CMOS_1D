function[Ef]=Fermi3D_np(Eg,T,meffn,meffp,Dop3D)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h    = 6.62606896E-34;              %% Planck constant [J.s]
hbar = h/(2*pi);
e    = 1.602176487E-19;             %% electron charge [C]
m0   = 9.10938188E-31;              %% electron mass [kg]
kB   = 1.3806488E-23;               %% Boltzmann's constant [J/K]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Dop3D==0
   Dop3D=1;
   Ef=-Eg/2;
   return
end

%%%%%%%%%%%%%%%%%%%%%%%%%% Fermi level at T=0K %%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Dop3D>0
    Ef0 =  hbar^2 / (2*e*meffn*m0)   *  (3*pi^2* Dop3D )^(2/3);
elseif Dop3D<0
    Ef0= -Eg - hbar^2 / (2*e*meffp*m0)   *   (-3*(pi^2)* Dop3D )^(2/3);
elseif Dop3D==0
    Ef0=0;
end

if T==0
    Ef=Ef0;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%% 3D density of states %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Here, I try to optimize the meshing
if Dop3D>0
    Emax=Ef0+0.3;
    Emin=-Eg;
else
    Emin=Ef0-0.3;
    Emax=0;
end

En1=linspace( 0 , Emax ,1e4);
En2=linspace( Emax , 3 ,50);
En=sort([En1 En2]);

Ep1=linspace( -Eg-3   , Emin , 50 );
Ep2=linspace( Emin , -Eg , 1e4 );
Ep=sort([Ep1 Ep2]);

ro3Dn = (1/(2*pi^2))*( (2*e*meffn*m0/(hbar^2))^(3/2) ) * sqrt(En) ;
ro3Dp = (1/(2*pi^2))*( (2*e*meffp*m0/(hbar^2))^(3/2) ) * sqrt( -(Ep+Eg) ) ;

%%%%%%%%%%%%%%%%%%%%% Fermi level at any temperature %%%%%%%%%%%%%%%%%%%%%%

Ef=Ef0;
FEc = 1./(1+exp(+(En-Ef)/(kB*T/e))) ;
FEv = 1./(1+exp(-(Ep-Ef)/(kB*T/e))) ; 

ro3DEfn=ro3Dn.*FEc;
NtotX=trapz(En,ro3DEfn);

ro3DEfp=ro3Dp.*FEv;
PtotX=trapz(Ep,ro3DEfp);

NPtotX=NtotX-PtotX;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the Fermi level is obiously at T=0K the max possible for n doped
% so now, it scans down with big step (ddE) 
% to find the Fermi level when T is not zero
    
ddE=0.005;   % step of the scan in eV to pass below the Fermi level
 
    
    
while ( ((NPtotX - Dop3D)) > 0)
    
    Ef  = Ef - ddE;    
    FEc = 1./(1+exp(+(En-Ef)/(kB*T/e))) ;
    FEv = 1./(1+exp(-(Ep-Ef)/(kB*T/e))) ; 

    ro3DEfn = ro3Dn.*FEc;
    NtotX   = trapz(En,ro3DEfn);
    ro3DEfp = ro3Dp.*FEv;
    PtotX   = trapz(Ep,ro3DEfp);

    NPtotX  = NtotX-PtotX;
                            
end         

 % the Fermi level is obiously at T=0K the min possible for p doped
 % so now, it scans up with big step (ddE) 
 % to found the Fermi level when T is not zero
 % step of the scan in eV to pass above the Fermi level
 
while ( ((NPtotX - Dop3D)) < 0)
        
    Ef  = Ef + ddE;    
    FEc = 1./(1+exp(+(En-Ef)/(kB*T/e))) ;
    FEv = 1./(1+exp(-(Ep-Ef)/(kB*T/e))) ; 

    ro3DEfn = ro3Dn.*FEc;
    NtotX   = trapz(En,ro3DEfn);
    ro3DEfp = ro3Dp.*FEv;
    PtotX   = trapz(Ep,ro3DEfp);

    NPtotX  = NtotX-PtotX;
            
end        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now, it will try to get as close as posible to the real Ef with an
% error of 1% by dichotomy

Ef1=Ef;
Ef2=Ef1+ddE;
Epsilon = 1e-10;
        
while  abs((NPtotX - Dop3D)/Dop3D) > Epsilon  % find the Fermi level at any temperature
   
    if NPtotX > Dop3D
        Ef= Ef - abs(Ef1-Ef2)/2 ;  
        Ef1=Ef ;
        FEc= 1./(1+exp((En-Ef)/(kB*T/e))) ;
        FEv = 1./(1+exp(-(Ep-Ef)/(kB*T/e))) ; 
    else
        Ef= Ef + abs(Ef1-Ef2)/2 ;
        Ef2=Ef ;
        FEc= 1./(1+exp((En-Ef)/(kB*T/e))) ;
        FEv = 1./(1+exp(-(Ep-Ef)/(kB*T/e))) ; 
    end

    ro3DEfn = ro3Dn.*FEc;
    NtotX   = trapz(En,ro3DEfn);
    ro3DEfp = ro3Dp.*FEv;
    PtotX   = trapz(Ep,ro3DEfp);

    NPtotX  = NtotX-PtotX;
       
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure('position',[150 30 1000 650]);
% plot(ro3Dn,En,'b')
% hold on;grid on;box on;
% plot(ro3Dp,Ep,'b')
% 
% plot(ro3DEfn,En,'r.')
% plot(ro3DEfp,Ep,'g.')
% 
% plot([0 max(max(ro3Dp),max(ro3Dn))],Ef*[1 1],'m')
% 
% xlabel('Densité d etat Bulk (eV-1.m-3)');
% ylabel('Energy (eV)');

end
