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

function [value, isterminal, direction] = combinedEvents(~, y, d0, Droplet, Pamb, liquidRho, liquidCp, EoS, min_fraction)
% Extract component masses (excluding droplet temperature)
mdi = y(2:end, 1);

% Condition 1: Check when any component's molar fraction is below the threshold
[~, Yli, Xli, ~] = getLiquidComposition(Droplet, mdi');

depletionCondition = any(Xli <= min_fraction);  % Any component's molar fraction < min_fraction

% Condition 2: Compute droplet diameter and check when it shrinks below 10% of initial diameter
[~, rhoL, ~, ~] = getLiquidProperties(Droplet, liquidRho, liquidCp, Pamb, y(1), Yli, EoS);
md = sum(mdi);  % Total droplet mass
[~, ~, d, ~] = getSize(md, rhoL);
shrinkageCondition = (d / d0) - 0.05;  % Becomes zero when d ≤ 5% d0

% The event is triggered when either condition is met
value = [depletionCondition; shrinkageCondition];

isterminal = [1; 1];  % Stop integration when either condition is met
direction  = [0; 0];  % Detect when mass is decreasing, diameter can cross either way
end