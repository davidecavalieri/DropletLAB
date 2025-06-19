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

function [Shi,Sh_bulk,Nu] = getRanzMarshall(Ambient,NaturalConvection,d,Vrel,Tamb,Td,Sci,Sc_bulk,Pr,Pamb,YviA_inf)

%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routine returns the Sherwood and Nusselt numbers of a sphere in
% forced convection

% Inputs:
% 1) NaturalConvection  : Struct variable containing the droplet species names,
% 2) d          : Droplet diameter    [m]
% 3) Vrel       : Relative velocity   [m/s]
% 4) Tamb       : Ambient temperature [K]
% 5) Td         : Droplet Temperature [K]
% 6) Sci        : Schmidt numbers
% 7) Pr         : Prandtl number

% Outputs:
% 1) Shi : Sherwood numbers              [-]
% 2) Nu  : Nusselt number                [-]

R0 = 8.314472;

% GRAVITY ACCELERATION
g0 = 9.81;

% CONSTANTS
k1 = 1/2;
k2 = 1/3;


Mi = cell2mat({Ambient.M});

% FILM MOLAR FRACTIONS AND AVERAGE MIXTURE WEIGHT
[Xi,M] = Mass2Mol(Mi,YviA_inf);

% AMBIENT PROPERTIES
rho_g    = Pamb*M/(R0*Tamb);
[mu_g,~] = calcMuCHEMKIN(Xi,Tamb,Ambient);

% REYNOLDS NUMBER
ReD = (rho_g*d*Vrel)/mu_g;


if strcmp(NaturalConvection,'off')

    % SHERWOOD NUMBERS
    Shi = 2  + 0.6.*(sqrt(ReD)).*(Sci.^k2);

    % BULK SHERWOOD NUMBER
    Sh_bulk = 2  + 0.6.*(sqrt(ReD)).*(Sc_bulk.^k2);

    % NUSSELT NUMBER
    Nu  = 2  + 0.6*(sqrt(ReD))*(Pr.^k2);

elseif strcmp(NaturalConvection,'on')

    % GRASHOF NUMBER
    GrD = (g0*(rho_g^2)*(Tamb - Td)*(d^3))/(Tamb*mu_g*mu_g);

    Shi = 2.000 + 0.6.* ((max(ReD,max(GrD,0).^k1)).^k1).*(Sci.^k2);

    Sh_bulk = 2.000  + 0.6.* ((max(ReD,max(GrD,0).^k1)).^k1).*(Sc_bulk.^k2);

    Nu  = 2.000 + 0.6.* ((max(ReD,max(GrD,0).^k1)).^k1).*(Pr.^k2);


else

    error('The selected Natural convection model is not available')

end

end