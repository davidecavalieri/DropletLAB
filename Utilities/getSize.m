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

function [Vd,rd,d,Ad]= getSize(md,rhoL)

%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routine returns the droplet characteristic sizes.

% Inputs:
% 1) md   : Total mass of the droplet             [kg]
% 2) rhoL : Liquid Mixture density of the droplet [kg/m^3]

% Outputs:
% 1) Vd : Droplet Volume   [m^3]
% 2) rd : Droplet Radius   [m]
% 3) d  : Droplet Diameter [m]
% 4) Ad : Droplet Area     [m^2]

% DROPLET VOLUME   [m^3]
Vd = md/rhoL;

% DROPLET RADIUS   [m]
rd = ((3/4) * (Vd/pi))^(1/3);

% DROPLET DIAMETER [m]
d = 2 *rd;

% DROPLET AREA     [m^2]
Ad = pi * d^2;

end
