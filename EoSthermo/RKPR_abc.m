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

function [ac,bc,cc] = RKPR_abc(Tc,Pc,omega,Zc)

%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RKPR_cbc calculates the non-ideal parameters of the Redlick-Kwong-Peng-Robinson EoS

% Inputs:
% 1) Tc    : Critical temperature [K]
% 2) Pc    : Critical pressure    [Pa]
% 3) omega : acentric factor      [-]
% 4) Zc    : Critical compressibility factor

% Outputs:
% ac : Attractive force term
% bc : Covolume
% cc : noName

R0 = 8.314472;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 1/3;

d1 = 0.428363;
d2 = 18.496215;
d3 = 0.338426;
d4 = 0.660000;
d5 = 789.723105;
d6 = 2.512392;

if Zc >= 0.289
    Zc = 0.289; % Only for the helium
end

delta1 = d1 + d2*((d3 - 1.168*Zc)^d4)  + d5*((d3 - 1.168*Zc)^d6);

d = (1 + delta1^2)/(1 + delta1);

y = 1 + (2*(1+delta1))^k  + (4/(1+delta1))^k;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A1 = -2.4407;
A0 = 0.0017;
B1 = 7.4513;
B0 = 1.9681;
C1 = 12.5040;
C0 = -2.7238;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc =   (1.168*Zc*A1 + A0)*(omega^2) ...
    + (1.168*Zc*B1 + B0)*omega ...
    + (1.168*Zc*C1 + C0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bc =  (1/(3*y + d -1))*((R0*Tc)/(Pc));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
acn = 3*y*y + 3*y*d + d*d + d -1;
acd = (3*y + d -1)^2;
ac =  (acn/acd)*((R0^2)*(Tc^2)/(Pc));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end