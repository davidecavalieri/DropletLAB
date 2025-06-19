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

function [Lvi] = getVapEnthalpy(Droplet,Ambient,vapEnthalpy,Xli,Xvi,T,Pamb)

%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routine returns the vaporization enthalpies of the droplet components

% Inputs:
% 1) Droplet      : Struct variable containing the droplet inputs parameters
% 2) VapEnthalpy  : String containing the selected Vaporization-Enthalpy model
% 3) T            : Droplet Temperature [K]

% Output:
% 1) Lvi : Vaporization enthalpy of each droplet's components [J/kg]

if strcmp(vapEnthalpy,'Coolprop')

    % RENAME THE DROPLET SPECIES 
    Droplet_fluids = {Droplet(:).Species};

    % NUMBER OF DROPLET SPECIES
    NsD = numel(cell2mat({Droplet.Tc}));

    % ROUND OF TEMPERATURE TO MATCH COOLPROP'S TOLERANCE
    T = round(T,3);

    % ALLOCATION
    Lvi = zeros(1,NsD);
    
    for i = 1:NsD
        
        H_V  = py.CoolProp.CoolProp.PropsSI('HMASS','T',T,'Q',1,char(Droplet_fluids(i)));
        H_L  = py.CoolProp.CoolProp.PropsSI('HMASS','T',T,'Q',0,char(Droplet_fluids(i)));
        Lvi(i) = H_V-H_L;    
  
    end
    
elseif strcmp(vapEnthalpy,'Yaws')
    
        [Lvi] = YawsPoly_Hvap(Droplet,T);

else
    
    error('The selected Vaporization enthalpy model is not available')
    
end

end