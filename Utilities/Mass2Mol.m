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

function [Xi,Mavg] = Mass2Mol(Mi,Yi)

%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routine performs the conversion of mass fractions to mole fractions.

% Inputs:
% 1) Wi: Molar weights  [kg/mol]
% 2) Yi: Mass fractions [-]

% Outputs:
% 1) Xi   : Mole fractions    [-]
% 2) Mavg : Mean Molar weight [Kg/mol]

% REVERSE MEAN MOLAR WEIGHT OF THE MIXTURE [mol/kg]
oneM  = sum(Yi./Mi);

% MEAN MOLAR WEIGHT OF THE MIXTURE [kg/mol]
Mavg  = 1/oneM;

% MOLE FRACTIONS OF THE MIXTURE [-]
Xi     = (Yi./Mi) * Mavg;

end