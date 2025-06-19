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

function [lambdaMix] = lambdaMixWassiljewa(mu,lambda,W,Xi)

% DESCRIPTION:
% This routine calculates the gas phase mixture conductivity  employing the
% Wassiljewa method

% INPUTS:
% 1) mu     : Single component viscosities    [Pa s]
% 2) lambda : Single component conductivities [W/m/K]
% 2) W      : Molar Weights                   [kg/mol]
% 3) Xi     : Mole fractions                  [-]

% OUTPUTS
% 1) lambdaMix : Mixture conductivity [W/m/K]

eps = 1.065;

% NUMBER OF SPECIES
Ns = length(Xi);

% ALLOCATIONS
Ratio_ij = zeros(Ns,Ns);
Aij   = zeros(Ns,Ns);
AijN  = zeros(Ns,Ns);
AijD  = zeros(Ns,Ns);

for i = 1:Ns
    for j= 1:Ns

        if  i==j
            Ratio_ij(i,j) = 1;
            Aij(i,j)      = 1;
        else

            Ratio_ij(i,j) = (mu(i)/mu(j))*(W(j)/W(i));

            AijN(i,j) = eps *(1 + sqrt(Ratio_ij(i,j)) * ((W(i)/W(j))^0.25))^2;
            AijD(i,j) = sqrt(8 *(1 + W(i)/W(j)));
            Aij(i,j)  = AijN(i,j)/AijD(i,j);

        end

    end

end

lambdaMixN = Xi.*lambda;
lambdaMixD = Xi.*Aij;
lambdaMixD = (sum(lambdaMixD,2))';
lambdaMix = sum(lambdaMixN./lambdaMixD);

end