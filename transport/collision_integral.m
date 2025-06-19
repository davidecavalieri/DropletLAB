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

function [SIGMA] = collision_integral(Tstar)

%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% collision_integral calculates the reduced collision integral with the
% Neufeld correlation (1972)

% Input:
% Tstar = Dimensionless temperature      [-]

% Output:
% SIGMA = reduced collision integral [-]

%% COEFFICIENTS.
%Chung et al. Generalized Multiparameter Correlation for Nonpolar and Polar Fluid
% Transport Properties. Ind. Eng. Chem. Res. 1988,27, 671-679

A = 1.16145;
B = 0.14874;
C = 0.52487;
D = 0.77320;
E = 2.16178;
F = 2.43787;
G = -6.435e-4;
H = 7.27371;
S = 18.0323;
W = -0.76830;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SIGMA =  (A/(Tstar.^B)) + C/(exp(D*Tstar)) + E/(exp(F*Tstar)) ...
    + (G*(Tstar.^B)).*sin(S*(Tstar^W) - H);

end
