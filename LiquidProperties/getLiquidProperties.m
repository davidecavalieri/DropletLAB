%    ___                __    __  __   ___   ___
%   / _ \_______  ___  / /__ / /_/ /  / _ | / _ )
%  / // / __/ _ \/ _ \/ / -_) __/ /__/ __ |/ _  |
% /____/_/  \___/ .__/_/\__/\__/____/_/ |_/____/
%              /_/

% ------------------------------------------------------------------------%
% Contributors / Copyright Notice
% © 2025 Davide Cavalieri — davide.cavalieri@uniroma1.it
% Postdoctoral Researcher @ Sapienza University of Rome,
% Department of Mechanical and Aerospace Engineering (DIMA)
%
% © 2025 Jacopo Liberatori — jacopo.liberatori@centralesupelec.fr
% Postdoctoral Researcher @ CentraleSupélec, Laboratoire EM2C (CNRS)
%
% © 2025 Matteo Blandino — matteo.blandino@uniroma1.it
% Ph.D. Student @ Sapienza University of Rome,
% Department of Mechanical and Aerospace Engineering (DIMA)
%
% Reference:
% D. Cavalieri, J. Liberatori, M. Blandino, P.E. Lapenna, M. Valorani, and P.P. Ciottoli.
% Evaluation of non‑ideal fluid modeling for droplet evaporation in jet‑engine‑like
% conditions. Flow, Turbulence and Combustion 114:3, pp. 857-885, 2024.
% DOI: 10.1007/s10494-024-00610-x.
% ------------------------------------------------------------------------%

function [rhoL_i,rhoL,cpL_i,cpL] = ...
    getLiquidProperties(Droplet,liquidRho,liquidCp,Pamb,T,Yli,EoS)

%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routine returns the droplet's liquid properties.

% Inputs:
% 1) Droplet    : Struct variable containing the droplet inputs parameters
% 2) liquidProp : String containing the selected liquid properties model
% 3) Pamb       : Ambient Pressure      [Pa]
% 4) T          : Droplet Temperature   [K]
% 5) Yli        : Liquid Mass Fractions [-]

% Outputs:
% 1) rhoL_i : Liquid Density of each droplet's components                [kg/m^3]
% 2) rhoL   : Liquid Mixture density of the droplet                      [kg/m^3]
% 3) cpL_i  : Liquid Isobaric-Specific-Heat of each droplet's components [J/kgK]
% 4) cpL    : Liquid Mixture Isobaric-Specific-Heat of the droplet       [J/kgK]

% REF:
% (i) http://www.coolprop.org/index.html

% NUMBER OF DROPLET SPECIES
NsD = numel(Yli);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DENSITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COOLPROP LIQUID DENSITIES
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

    % LIQUID MIXTURE DENSITY
    rhoL = sum(Yli./rhoL_i)^-1;

elseif strcmp(liquidRho,'Yaws')

    [rhoL_i,rhoL,~,~] = YawsPoly_Liq(Droplet,T,Yli);

elseif strcmp(liquidRho,'EoS')

    % Note: Here we compute the liquid mixture density
    [rhoL,~,~,~] = RFproperties(Pamb,T,Yli,Droplet,EoS);

    % ALLOCATIONS
    rhoL_i = zeros(1,NsD);

    % Note: Here we compute the single component properties (needed
    % for the initialization!)
    Composition = zeros(1, NsD);
    for i = 1:NsD
        Composition(i) = 1;
        [rhoL_i(i),~,~,~] = RFproperties(Pamb,T,Composition,Droplet,EoS);
        Composition(i) = 0;
    end

else

    error('The selected Liquid density model is not available')

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COOLPROP LIQUID DENSITIES
if strcmp(liquidCp,'Coolprop')

    % RENAME THE DROPLET SPECIES
    Droplet_fluids = {Droplet(:).Species};

    % ALLOCATIONS
    cpL_i = zeros(1,NsD);

    % ROUND OF TEMPERATURE TO MATCH COOLPROP'S TOLERANCE
    T = round(T,3);

    for i = 1:NsD
        cpL_i(i)  = py.CoolProp.CoolProp.PropsSI('CPMASS','T',T,'P',Pamb,char(Droplet_fluids(i)));
    end

    % LIQUID MIXTURE ISOBARIC SPECIFIC HEAT
    cpL  = sum(Yli.*cpL_i);

elseif strcmp(liquidCp,'Yaws')

    [~,~,cpL_i,cpL] = YawsPoly_Liq(Droplet,T,Yli);

elseif strcmp(liquidCp,'EoS')

    % Note: Here we compute the liquid mixture cp
    [~,cpL,~,~] = RFproperties(Pamb,T,Yli,Droplet,EoS);

    % ALLOCATIONS
    cpL_i = zeros(1,NsD);

    % Note: Here we compute the single component properties (needed
    % for the initialization!)
    Composition = zeros(1, NsD);
    for i = 1:NsD
        Composition(i) = 1;
        [~,cpL_i(i),~,~] = RFproperties(Pamb,T,Composition,Droplet,EoS);
        Composition(i) = 0;
    end

else

    error('The selected Liquid cp model is not available')

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end