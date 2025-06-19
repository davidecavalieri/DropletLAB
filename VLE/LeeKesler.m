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

function [Psat_i] = LeeKesler(Droplet,T)

%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routine returns the saturated pressures of the droplet's components
% using the Lee-Kesler correlation

% Inputs:
% 1) Droplet  : Struct variable
% 2) T        : Droplet temperature [K]

% Output:
% 1) Psat_i   : Saturated Pressures [Pa]

%REF:
% (i) B. E. Poling, J. M. Prausnitz, J. P. O’Connell. Properties of Gases
% and Liquids, 5th Edition, 2001, McGraw-Hill Education

% CRITICAL PRESSURES [Pa]
Pci = cell2mat({Droplet.Pc});

% CRITICAL TEMPERAURES [K]
Tci = cell2mat({Droplet.Tc});

% ACENTRIC FACTORS [-]
omegai = cell2mat({Droplet.omega});

% REDUCED TEMPERATURES [-]
Tri = T./Tci;

% PITZER RELATIONS [-]
f0 = 5.92714 - 6.09648./Tri - 1.28862*log(Tri) + 0.169347*(Tri.^6);
f1 = 15.2518 - 15.6875./Tri - 13.4721*log(Tri) + 0.43577*(Tri.^6);

% SATURATED PRESSURES
Prsat = exp(f0 + omegai.*f1);
Psat_i = Prsat.*Pci;

end