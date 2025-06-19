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

function [md,Yli,Xli,Ml] = getLiquidComposition(Droplet,mdi)

%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routine returns the liquid phase composition of the droplet.

% Inputs:
% 1) Droplet: Struct variable
% 2) mdi:     Masses of droplet components [kg]

% Outputs:
% 1) md  : Total mass of the droplet             [kg]
% 2) Yli : Liquid mass fractions                 [-]
% 3) Xli : Liquid mole fractions                 [-]
% 4) Ml  : Mean Molar weight of the liquid phase [Kg/mol]

% MOLAR WEIGHT OF THE DROPLET COMPONENTS [kg/mol]
MiD = cell2mat({Droplet.M});

% TOTAL MASS OF THE DROPLET   [kg]
md = sum(mdi);

% LIQUID PHASE MASS FRACTIONS [-]
Yli = mdi./md;

% LIQUID PHASE MOLE FRACTIONS [-]
[Xli,Ml] = Mass2Mol(MiD,Yli);

end