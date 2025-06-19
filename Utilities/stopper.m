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

function [value, isterminal, direction] = stopper(y, d0, Droplet, Pamb, liquidRho, liquidCp, EoS)
    % DESCRIPTION:
    % This routine stops the integration of the ODE system when the droplet diameter 
    % is 10% of the initial diameter.

    % STATE VECTOR VARIABLES
    Td  = y(1, 1);
    mdi = y(2:end, 1);

    % RETRIEVE THE LIQUID COMPOSITION 
    [~, Yli, ~, ~] = getLiquidComposition(Droplet, mdi');

    % COMPUTE THE LIQUID MIXTURE DENSITY                            
    [~, rhoL, ~, ~] = getLiquidProperties(Droplet, liquidRho, liquidCp, Pamb, Td, Yli, EoS);

    % COMPUTE THE TOTAL MASS OF THE DROPLET
    md = sum(y(2:end, 1));

    % COMPUTE THE DROPLET DIAMETER                                       
    [~, ~, d, ~] = getSize(md, rhoL);

    % THRESHOLD VALUE
    value = (d / d0 <= 0.1);

    isterminal = 1;   % Stop the integration
    direction  = 0;   % Detect zero-crossing in any direction
end