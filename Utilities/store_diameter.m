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

function status = store_diameter(t, w, flag, Droplet, liquidRho, liquidCp, Pamb, EoS)
    persistent diameters times  % Persistent storage across function calls

    if isempty(flag)  % Regular solver call
        % If multiple time steps are passed at once, loop through each step
        for i = 1:length(t)
            t_i = t(i);  % Single time point
            w_i = w(:, i);  % Corresponding state vector for that time step
            
            % Extract mass values from state vector (excluding temperature)
            mdi = w_i(2:end);  % Ensure correct shape

            % Compute liquid composition
            [~, Yli, ~, ~] = getLiquidComposition(Droplet, mdi');

            % Compute liquid density
            [~, rhoL, ~, ~] = getLiquidProperties(Droplet, liquidRho, liquidCp, Pamb, w_i(1), Yli, EoS);

            % Compute the droplet diameter
            [~, ~, d, ~] = getSize(sum(mdi), rhoL);

            % Store the results for each time step
            times = [times; t_i];  
            diameters = [diameters; d];
        end

    elseif strcmp(flag, 'init')  % Initialization
        diameters = [];
        times = [];
        
    elseif strcmp(flag, 'done')  % Finalize and store results
        assignin('base', 'd_times', times);
        assignin('base', 'd_values', diameters);
    end

    status = 0;  % Continue integration
end