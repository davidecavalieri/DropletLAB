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

function [muMix] = muMixWilke(mu,W,Xi)

%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routine calculates the gas phase mixture viscosities  employing the
% Wilke method

% Inputs:
% 1) mu : Single component viscosities [-]
% 2) W  : Molar Weights                [kg/mol]
% 3) Xi : Mole fractions               [-]

% Outputs
% 1) muMix : Mixture viscosity [Pa*s]

% REFs.
% C. R. Wilke, Journal of Chemical Physics 18:517 (1950).
% R. B. Bird, W. E. Stewart, and E. N. Lightfoot, Transport Phenomena, John Wiley and Sons,
% New York, 2002.


% NUMBER OF SPECIES
Ns = length(Xi);

% ALLOCATIONS
Ratio_ij = zeros(Ns,Ns);
phi_ij   = zeros(Ns,Ns);
phi_ijN  = zeros(Ns,Ns);
phi_ijD  = zeros(Ns,Ns);

for i = 1:Ns
    for j= 1:Ns

        if  i==j
            Ratio_ij(i,j) = 1;
            phi_ij(i,j)   = 1;
        else

            Ratio_ij(i,j) = (mu(i)/mu(j));
            phi_ijN(i,j) = (1 + sqrt(Ratio_ij(i,j)) * ((W(i)/W(j))^0.25))^2;
            phi_ijD(i,j) = sqrt(8 *(1 + W(i)/W(j)));
            phi_ij(i,j)  = phi_ijN(i,j)/phi_ijD(i,j);

        end

    end

end

muMixN = Xi.*mu;
muMixD = Xi.*phi_ij;
muMixD = (sum(muMixD,2))';
muMix = sum(muMixN./muMixD);

end