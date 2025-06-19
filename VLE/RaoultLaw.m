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

function [Xvsi,Xgsi,Yvsi,Ygsi] = RaoultLaw(Droplet,Ambient,Pamb,Xli,Psat_i,XviA_inf,YviA_inf) 

%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routine returns the droplet surface composition using Raoult's law

% Inputs:
% 1) Droplet : Struct variable containing the droplet inputs parameters
% 2) Ambient : Struct variable containing the droplet inputs parameters
% 3) Pamb    : Ambient pressure                             [Pa] 
% 3) Xli     : Liquid mole fractions of the droplet species [-]
% 4) Psati   : Saturated pressure of the droplet species    [Pa]

% Outputs:
% 1) Xvsi : Surface vapour mole fractions (droplet species)
% 2) Xgsi : Surface vapour mole fractions (Ambient species)
% 3) Yvsi : Surface vapour mass fractions (droplet species)
% 4) Ygsi : Surface vapour mass fractions (Ambient species)

% MOLAR WEIGHT (kg/mol)
MiD = cell2mat({Droplet.M});
MiA = cell2mat({Ambient.M});

% VAPOR MOLE FRACTIONS OF THE DROPLET COMPONENTS @SURFACE
Xvsi = (Xli.*Psat_i)./Pamb;

% TOTAL MOLE FRACTION OF THE DROPLET COMPONENTS @SURFACE
Xvs = sum(Xvsi);

% TOTAL MOLE FRACTION OF THE AMBIENT GAS  @SURFACE
Xgs  = 1 - Xvs; 

% MOLE FRACTIONS OF THE AMBIENT GAS MIXTURE @ SURFACE (IF MULTICOMPONENT e.g Air)
Xgsi = Xgs.*XviA_inf;

% VAPOUR MASS FRACTIONS OF THE DROPLET COMPONENTS @ SURFACE     
YvsiN = Xvsi.*MiD;
Yvsi  = YvsiN./(sum(YvsiN) + sum(Xgsi.*MiA));

% TOTAL AMBIENT VAPOUR MASS FRACTION @ SURFACE 
Ygs  = 1 - sum(Yvsi); 

% MOLE FRACTIONS OF THE AMBIENT GAS MIXTURE @ SURFACE (IF MULTICOMPONENT e.g Air)
Ygsi = Ygs.*YviA_inf;

end