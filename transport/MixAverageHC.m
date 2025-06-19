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

function [Di_mix] = MixAverageHC(Xi,Dij,Yi)

%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MixAverageHC calculates the mixture diffusion coeffiecient for species k
% into the mixture, using the mixture averaged formula.

% Inputs:
% 1) Xi              : Mole fractions [-]
% 2) Dij             : Matrix of binary mass diffusion coefficients [m^2/s]
% 3) Yi              : Mass fractions [-]

% Outputs:
% 1) Di_mix          : Mixture-averaged mass diffusion coefficients [m^2/s]

% Initialize output
Di_mix = zeros(1, numel(Yi));

for i = 1:numel(Yi)
    % SUM INITIALIZATION
    SUM = 0;
    for j = 1:numel(Xi)
        if j ~= i
            SUM = SUM + Xi(j) / Dij(j, i);
        end
    end
    Di_mix(i) = (1- Yi(i))/SUM;
end

end
