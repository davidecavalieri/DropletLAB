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

function [alpha, dalphadT, ddalphadT] = alphaFunction(T,Tc,cc,EoS)

%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alphaFunction calculates the non-ideal temperature dependend parameter
% and his derivatives of a cubic EoS

% Inputs:
% 1) T  : Temperature [K]
% 2) Tc : Critical temperature [K]
% 3) cc : noName parameter

% Outputs:
% 1) alpha     : Temperature dependent parameter
% 2) dalphadT  : First order partial derivative of alpha
% 3) ddalphadT : Second order partial derivative of alpha

if strcmp(EoS,'PR') || strcmp(EoS,'SRK')

    alpha     =   (1 + cc*(1 - sqrt(T/Tc))).^2;
    dalphadT  = - cc * sqrt(alpha)/ sqrt( Tc * T );
    ddalphadT =  (cc^2/(2*T*Tc)) + (cc/(2*sqrt(T^3*Tc))*(1 + cc*(1 - sqrt(T/Tc))));

elseif strcmp(EoS,'RKPR')

    alpha     =  (3/(2 + T/Tc))^cc;
    dalphadT  =  -(cc*(3^cc))/((Tc*(2 + T/Tc)^(cc+1)));
    ddalphadT =  ((3^cc)*cc*(cc+1))/((Tc*Tc)*((2+T/Tc)^(cc+2)));

else
    error('The selected EoS is not available')

end

end