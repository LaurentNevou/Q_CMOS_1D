%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first column is the material used from the "library"
% second column is the length of the layer in nm
% third column is the n doping volumique of that layer in 1E18 cm-3 
% fourth column is the p doping volumique of that layer in 1E18 cm-3 
% fifth column is the number of points (meshing) of that layer 
% You have to put a resonable amount of doping! Otherwise, it will diverge!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 5V NMOS with 12nm oxide gate

M = [
Poly    3    50       0      3
Poly    2    50       0     15
Oxide  12     0       0      2
Si     10     0      0.2    40
Si    100     0      0.2    30
];

Fermi_layerbreak_L = 2;   %% it chooses the layer number at which the Fermi level breaks on the left
Fermi_layerbreak_R = 3;   %% it chooses the layer number at which the Fermi level breaks on the right

capa_layerbreak_L  = 3;   %% it chooses the layer number at which the charge starts to be counted (0 is OK)
capa_layerbreak_R  = 5;   %% it chooses the layer number at which the charge stops to be counted

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 5V PMOS with 12nm oxide gate

% M=[
% Poly    4    0       50      5
% Poly    3    0       50     10
% Oxide  12    0        0      3
% Si     10    0.2      0     50
% Si     95    0.2      0     30
% ];
% 
% Fermi_layerbreak_L = 2;   %% it chooses the layer number at which the Fermi level breaks on the left
% Fermi_layerbreak_R = 3;   %% it chooses the layer number at which the Fermi level breaks on the right
% 
% capa_layerbreak_L  = 3;   %% it chooses the layer number at which the charge starts to be counted (0 is OK)
% capa_layerbreak_R  = 5;   %% it chooses the layer number at which the charge stops to be counted

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.8V NMOS with 3.5nm oxide gate


% M = [
% Poly   3    50       0     5
% Poly   2    50       0     25
% Oxide    3   0       0      3
% Si    10    0      0.3     50
% Si    70    0      0.3     30
% ];
% 
% Fermi_layerbreak_L = 2;   %% it chooses the layer number at which the Fermi level breaks on the left
% Fermi_layerbreak_R = 3;   %% it chooses the layer number at which the Fermi level breaks on the right
% 
% capa_layerbreak_L  = 3;   %% it chooses the layer number at which the charge starts to be counted (0 is OK)
% capa_layerbreak_R  = 5;   %% it chooses the layer number at which the charge stops to be counted

