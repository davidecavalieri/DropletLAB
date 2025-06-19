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

function [Dij] = Fuller(Species,p, T)

%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fuller compute the binary mass diffusion coefficients using the Fuller
% method.

% Inputs:
% 1) Species             : Struct variable containing the Species inputs parameters
% 2) p                   : pressure [Pa]
% 3) T                   : Temperature [K]

% Output:
% Dij               : Matrix of the binary mass diffusion coefficients [m^2/s]

gram     = 1e3;        % Kg to grams
a        = 1/3;        % Atomic diffusion coefficient exponent
aT       = 1.75;       % Temperature exponent
PasTobar = 1e-05;      % Pascal to bar
cmTom    = 1e-04;      % cm2 to m2
Pbar     = p*PasTobar;
%--------------------------------------------------------------------------
Mi = cell2mat({Species.M});
viD = cell2mat({Species.vd}); % Già sommati
Migram   = Mi*gram;
viDexp   = viD.^a;

%% BINARY DIFFUSION COEFFICIENTS (Fuller method)
Mij      = zeros(numel(Mi),numel(Mi));
sum_vij  = zeros(numel(Mi),numel(Mi));
Dij_den  = zeros(numel(Mi),numel(Mi));
Dij      = zeros(numel(Mi),numel(Mi));

for i = 1:numel(Mi)
    for j= 1:numel(Mi)

        if  i==j
            Dij(i,j) = 1e20;
        else
            Mij(i,j)      = 2*((Migram(i) * Migram(j))/(Migram(i) + Migram(j)));
            sum_vij(i,j)  = (viDexp(i) + viDexp(j))^2;
            Dij_den(i,j)  = (Pbar) * (sqrt(Mij(i,j))) * (sum_vij(i,j));
            Dij(i,j)      = ((0.00143 * (T^aT))/(Dij_den(i,j)))*cmTom;
        end
    end
end

end