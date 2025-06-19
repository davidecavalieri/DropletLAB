%    ___                __    __  __   ___   ___
%   / _ \_______  ___  / /__ / /_/ /  / _ | / _ )
%  / // / __/ _ \/ _ \/ / -_) __/ /__/ __ |/ _  |
% /____/_/  \___/ .__/_/\__/\__/____/_/ |_/____/
%              /_/

% ------------------------------------------------------------------------%
% Contributors / Copyright Notice
% © 2025 Davide Cavalieri - davide.cavalieri@uniroma1.it
% Postdoctoral Researcher @ Sapienza University of Rome,
% Department of Mechanical and Aerospace Engineering (DIMA)
%
% © 2025 Jacopo Liberatori - jacopo.liberatori@centralesupelec.fr
% Postdoctoral Researcher @ CentraleSupélec, Laboratoire EM2C (CNRS)
%
% © 2025 Matteo Blandino - matteo.blandino@uniroma1.it
% Ph.D. Student @ Sapienza University of Rome,
% Department of Mechanical and Aerospace Engineering (DIMA)
%
% Reference:
% D. Cavalieri, J. Liberatori, M. Blandino, P.E. Lapenna, M. Valorani, and P.P. Ciottoli.
% Evaluation of non-ideal fluid modeling for droplet evaporation in jet-engine-like
% conditions. Flow, Turbulence and Combustion 114:3, pp. 857-885, 2024.
% DOI: 10.1007/s10494-024-00610-x.
% ------------------------------------------------------------------------%

function [Lvi] = YawsPoly_Hvap(Droplet,T)

%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routine returns the latent heat of evaporation of the droplet components
% using polynomial models based.

% Inputs:
% 1) Droplet    : Struct variable containing the droplet species names, critical
% 2) T          : Droplet Temperature   [K]

% Output:
% 1) Lvi        : Vaporization Enthalpy [J/kg]

% REFs:
% (i) The Yaws Handbook of Physical Properties for Hydrocarbons and Chemicals: 
% Physical Properties for More Than 54,000 Organic and Inorganic Chemical 
% Compounds,Coverage for C1 to C100 Organics and Ac to Zr Inorganics, 
% Gulf Professional Publishing (2015) 

    % NUMBER OF DROPLET SPECIES
    NsD = numel(cell2mat({Droplet.Tc}));

    % ALLOCATION
    Tr = zeros(1,NsD);
    Lvi = zeros(1,NsD);

for i = 1:NsD

    % REDUCED TEMPERATURE [-]
    Tr(i) = T/Droplet(i).Tc;

    Lvi(i) = 1e3*(Droplet(i).Ah*((1 - Tr(i))^Droplet(i).nh))/(Droplet(i).M);
    
end

end