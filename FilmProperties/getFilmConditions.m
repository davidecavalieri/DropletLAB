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

function [Tf,Yfi,Yfai] = getFilmConditions(Tamb,Ar,T,Ysi,Yvi_inf,YviA_inf)

%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routine returns the reference conditions of the film.

% Inputs:
% 1) Tamb   : Ambient temperature                        [K]
% 2) Ar     : Weight coefficient                         [-]
% 3) T      : Droplet Temperature                        [K]
% 4) Ysi    : Surface vapor mass fractions  [-]
% 5) Yv_inf : Ambient vapor fractions (Droplets species) [-]

% Outputs:
% 1) Tf   : Film Temperature                        [K]
% 2) Yfi  : Film vapor fractions (Droplets species) [-]
% 3) Yfai : Film vapor fractions (Ambient species)  [-]

NsD = numel(Yvi_inf);

% FILM TEMPERATURE
Tf   = T  + Ar*(Tamb-T);

% FILM MASS FRACTIONS (Droplet species)
Yfi  = Ysi(1:NsD) + Ar*(Yvi_inf -Ysi(1:NsD));

% FILM MASS FRACTIONS (Ambient species)
Yfai = Ysi(NsD+1:end) + Ar*(YviA_inf -Ysi(NsD+1:end));


end