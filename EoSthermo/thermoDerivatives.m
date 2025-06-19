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

function [dpdT_v,dpdv_T] = thermoDerivatives(Species,am,bm,dadTm,v,T,X,EoS)

%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% thermoDerivatives calculates analytical thermodynamic derivatives for a
% given cubic equation of state
% Note: Most of them are commented (see below the code)

% Inputs:
% 1) Species     : Struct variable containing the Species inputs parameters
% 2) am     : Mixture attractive force term, am(X,T)
% 3) bm     : Mixture Covolume, bm(X,T)
% 4) dadTm  : Mixture first order partial derivative of am w.r.t T, dadTm(X,T)
% 5) v      : Molar volume
% 6) T      : Temperature [K]

% Outputs:
% 1) dpdT_v  : Partial derivatve of pressure w.r.t temperature @ fixed v
% 2) dpdv_T  : Partial derivatve of pressure w.r.t volume @ fixed T

R0 = 8.314472;

% NUMBER OF SPECIES (DROPLET+AMBIENT)
Ns = numel(X);

%%%%%%%%%%%%%%%%%%% GENERALIZED CUBIC EQUATION OF STATE %%%%%%%%%%%%%%%%%%%
if strcmp(EoS,'PR')
    delta1 = 1.0 + sqrt(2); delta2 = 1.0 - sqrt(2);
elseif strcmp(EoS,'SRK')
    delta1 = 1.0;           delta2 = 0;
elseif strcmp(EoS,'RKPR')

    d1 = 0.428363; d2 = 18.496215;  d3 = 0.338426;
    d4 = 0.660000; d5 = 789.723105; d6 = 2.512392;

    for k = 1:Ns
        if Species(k).Zc >= 0.289
            Species(k).Zc  = 0.289; % Only for the helium
        end
        delta1_i(k) = d1 + d2*((d3 - 1.168*Species(k).Zc)^d4) + d5*((d3 - 1.168*Species(k).Zc)^d6);
        delta2_i(k) = (1-delta1_i(k))/(1+delta1_i(k));
    end
    delta1 = sum(X.*delta1_i,2);
    delta2 = sum(X.*delta2_i,2);

else
    error('The selected EoS is not available')
end

D = (v + delta1*bm)*(v + delta2*bm);

dpdT_v = R0/(v-bm) - dadTm/D;

dpdv_T = - R0 * T / ( v - bm )^2  + (am*(2*v + (delta1 + delta2)*bm))/(D^2);

end