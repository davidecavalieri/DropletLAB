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

function [rhoL_i,rhoL,cpL_i,cpL] = YawsPoly_Liq(Droplet,T,Yli)

%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routine returns the droplet liquid properties using polynomial
% models

% INPUTS:
% 1) Droplet    : Struct variable containing the droplet inputs parameters
% 2) T          : Droplet Temperature   [K]
% 3) Yli        : Liquid Mass Fractions [-]

% OUTPUTS
% 1) rhoL_i : Liquid Density of each droplet's components                [kg/m^3]
% 2) rhoL   : Liquid Mixture density of the droplet                      [kg/m^3]
% 3) cpL_i  : Liquid Isobaric-Specific-Heat of each droplet's components [J/kgK]
% 4) cpL    : Liquid Mixture Isobaric-Specific-Heat of the droplet       [J/kgK]

% REFs:
% (i) The Yaws Handbook of Physical Properties for Hydrocarbons and Chemicals:
% Physical Properties for More Than 54,000 Organic and Inorganic Chemical
% Compounds,Coverage for C1 to C100 Organics and Ac to Zr Inorganics,
% Gulf Professional Publishing (2015)

% (ii) B. E. Poling, J. M. Prausnitz, J. P. O’Connell. Properties of Gases
% and Liquids, 5th Edition, 2001, McGraw-Hill Education

% NUMBER OF DROPLET SPECIES
NsD = numel(Yli);

% ALLOCATIONS
rhoL_i = zeros(1,NsD);
cpL_i  = zeros(1,NsD);

for i = 1:NsD

    rhoL_i(i) = 1000*Droplet(i).Arho*Droplet(i).Brho^(-(1-T/Droplet(i).Tc)^Droplet(i).n);

    cpL_i(i)  = (1/Droplet(i).M)*(Droplet(i).Acp + Droplet(i).Bcp*T ...
        + Droplet(i).Ccp*T*T + Droplet(i).Dcp*T*T*T);

end

% LIQUID MIXTURE DENSITY
rhoL = sum(Yli./rhoL_i).^-1;

% LIQUID MIXTURE ISOBARIC SPECIFIC HEAT
cpL  = sum(Yli.*cpL_i);

end