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

function [Xvsi,Xgsi,Yvsi,Ygsi] = VLE_EoS(Droplet,Ambient,Pamb,Xli,T,Psati,XviA_inf,YviA_inf,EoS,activeComponents,t) 

%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routine returns the droplet surface composition using the VLE
% concept with a cubic EoS and neglecting the solubility effect of ambient
% inert gas at the droplet surface

% Inputs:
% 1) Droplet : Struct variable containing the droplet inputs parameters
% 2) Ambient : Struct variable containing the droplet inputs parameters
% 3) Pamb    : Ambient pressure                             [Pa] 
% 3) Xli     : Liquid mole fractions of the droplet species [-]
% 4) Psati   : Saturated pressure of the droplet species    [Pa]

% Outputs:
% 1) Xvsi : Surface vapour mole fractions (droplet species)
% 2) Xgsi : Surface vapour mole fractions (Ambient species)
% 3) Yvsi : Surface vapour mass fractions (droplet species)
% 4) Ygsi : Surface vapour mass fractions (Ambient species)

persistent Xvsi_prev

% Reset Xvsi_prev if necessary
if t == 0
    Xvsi_prev = [];
end

Species = [Droplet Ambient];

% INITIAL GUESS OF THE VAPOUR PHASE COMPOSITION
% First time step: estimate Xvsi_guess using Raoult’s law
if isempty(Xvsi_prev)
    [Xvsi_guess,~,~,~] = RaoultLaw(Droplet,Ambient,Pamb,Xli,Psati,XviA_inf,YviA_inf) ;
else
    if numel(Xvsi_prev) > numel(Xli)
        Xvsi_prev = Xvsi_prev(activeComponents);
    end
    Xvsi_guess = Xvsi_prev; % Use previous solution
end

% NUMBER OF DROPLET SPECIES = NUMBER OF EQUATION NEEDED FOR THE VLE, when the solubility effect is neglected
Ns = numel(Xli);

% STATE VECTOR FOR THE VLE PROBLEM:VAPOR PHASE COMPOSITION OF THE DROPLET
% SPECIES. 
W = Xvsi_guess;

options = optimset('TolFun',1e-12,'display','off','Algorithm','trust-region-reflective'); 

VLE     = @(W) objFunctionVLE(W,Xli,Species,T,Pamb,Ns,EoS,XviA_inf);
[W,~,~] = fsolve(VLE,W,options);

% VAPOR PHASE COMPOSITION (DROPLET SPECIES)
Xvsi = W;
% Store the solution for the next time step
Xvsi_prev = Xvsi;

% VAPOR PHASE COMPOSITION OF AMBIENT SPECIES
Xgsi = (1 - sum(Xvsi)).*XviA_inf;

% MOLAR WEIGHT (kg/mol)
MiD = cell2mat({Droplet.M});
MiA = cell2mat({Ambient.M});

% VAPOUR MASS FRACTIONS OF THE DROPLET COMPONENTS @ SURFACE     
YvsiN = Xvsi.*MiD;
Yvsi  = YvsiN./(sum(YvsiN) + sum(Xgsi.*MiA));

% TOTAL AMBIENT VAPOUR MASS FRACTION @ SURFACE 
Ygs  = 1 - sum(Yvsi); 

% MOLE FRACTIONS OF THE AMBIENT GAS MIXTURE @ SURFACE (IF MULTICOMPONENT e.g Air)
Ygsi = Ygs.*YviA_inf;

end