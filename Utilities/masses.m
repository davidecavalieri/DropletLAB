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

function [mi] = masses(Droplet,d,Xli,T,liquidRho,Pamb, EoS)

%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routine returns the i-mass of the droplet components.

% Inputs:
% 1) Droplet    : Struct variable
% 2) d   : Droplet diameter           [m]
% 3) Xli : Liquid Mole Fractions      [-]
% 4) T   : Droplet Temperature        [K]

% Outputs:
% 1) mi : i-mass of the droplet components [kg]

% NUMBER OF DROPLET SPECIES
NsD = numel(Xli);

% DROPLET SPECIES MOLECULAR WEIGHT
MiD = cell2mat({Droplet.M});

% MOLE 2 MASS CONVERSION
[Yli] = Mole2Mass(MiD,Xli);

% RADIUS
r = d/2;

% DROPLET VOLUME
V = (4/3)*pi*(r^3);

% LIQUID DENSITIES (MASS)
if strcmp(liquidRho,'Coolprop')

    % RENAME THE DROPLET SPECIES
    Droplet_fluids = {Droplet(:).Species};

    % ALLOCATIONS
    rhoL_i = zeros(1,NsD);

    % ROUND OF TEMPERATURE TO MATCH COOLPROP'S TOLERANCE
    T = round(T,3);

    for i = 1:NsD
        rhoL_i(i) = py.CoolProp.CoolProp.PropsSI('DMASS','T',T,'P',Pamb,char(Droplet_fluids(i)));
    end

elseif strcmp(liquidRho,'Yaws')

    [rhoL_i,~,~,~] = YawsPoly_Liq(Droplet,T,Yli);

elseif strcmp(liquidRho,'EoS')

    % Note: Here we compute the single component properties (needed
    % for the initialization!)
    Composition = zeros(1, NsD);
    % ALLOCATIONS
    rhoL_i = zeros(1,NsD);
    for i = 1:NsD
        Composition(i) = 1;
        [rhoL_i(i),~,~,~] = RFproperties(Pamb,T,Composition,Droplet,EoS);
        Composition(i) = 0;
    end

else

    error('The selected Liquid properties model is not available')

end

% LIQUID DENSITIES (MOLAR)
rhomol_i = rhoL_i./(MiD);

% SPECIFIC VOLUME
vmi = 1./rhomol_i;

% MASSES OF EACH COMPONENT
mi = (V/dot(vmi,Xli))*Xli.*MiD;

end
