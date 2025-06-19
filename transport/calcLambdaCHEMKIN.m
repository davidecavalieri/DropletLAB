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

function [lambda] = calcLambdaCHEMKIN(P,Xi,T,cpi,mui,Species)

%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routine calculates the individual component conductivities using the
% standard kinetic theory formulation. The mixture conductivity is then calculated
% using the mixture averaged formula.

% Inputs:
% 1) Xi  : Mole Fractions
% 2) T   : Temperature
% 3) Species:  Struct variable containing the Species inputs parameters

% Output:
% 1) lambda  : Mixture thermal conductivity [w/mK]

R0 = 8.314472;
% LENNARD-JONES POTENTIAL WELL DEPTH [K]
eps_kb  = cell2mat({Species.eps_kb});
% REDUCED TEMPERATURE [-]
KT_eps  = T./eps_kb;
% MOLECULAR WEIGHT [Kg mol-1]
Mi      = cell2mat({Species.M});
% LENNARD-JONES COLLISION DIAMETER [m]
sigma   = (cell2mat({Species.sigma}))*1e-10;
% GEOMETRY TYPE (0,1 or 2)
GeoType = cell2mat({Species.geometry});
% GAS CONSTANT [J/kgK]
R       = R0./Mi;
% ROTATIONAL RELAXATION COLLISION NUMBER AT 298 K
Zrot0   = cell2mat({Species.Zrot});
% BOLTZAMANN CONSTANT [m2 kg s-2 K-1]
kb = 1.3806488e-23;
% AVOGADRO'S NUMBER [mol-1]
NA = 6.02214076e23;
% SPECIFIC HEAT AT CONSTANT VOLUME TRANSITIONAL [-]
cv_trans = 3/2;

% COLLISION INTEGRAL DATA, Bird et. al. (2002)
KT_eps_bird = [0.3:0.05:2 2.1:0.1:2.6 2.7:0.1:5 6:10 12:2:20 25:5:40 50 70 100];
OMEGA_Diff  = [2.649 2.468 2.314 2.182 2.066 1.965 1.877 1.799 1.729 1.667 ...
    1.612 1.562 1.517 1.477 1.44 1.406 1.375 1.347 1.32 1.296 ...
    1.274 1.253 1.234 1.216 1.199 1.183 1.168 1.154 1.141 1.128 ...
    1.117 1.105 1.095 1.085 1.075 1.058 1.042 1.027 1.013 1.0006 ...
    0.989 0.9782 0.9682 0.9588 0.9500 0.9418 0.9340 0.9267 0.9197 ...
    0.9131 0.9068 0.9008 0.8952 0.8897 0.8845 0.8796 0.8748 0.8703 ...
    0.8659 0.8617 0.8576 0.8537 0.8499 0.8463 0.8428 0.8129 0.7898 ...
    0.7711 0.7555 0.7422 0.7202 0.7025 0.6878 0.6751 0.6640 0.6414 ...
    0.6235 0.6088 0.5964 0.5763 0.5415 0.5180];

% ROTATIONAL RELAXATION COLLISION NUMBER [-]
eps_kb0 = eps_kb./298;
invTred = eps_kb./T;
% PARKER BRAU & JONKMAN
F_298 = 1+((pi^1.5)/2).*sqrt(eps_kb0)+((pi^2)/4 + 2).*eps_kb0+(pi^1.5).*(eps_kb0.^1.5);
F = 1+((pi^1.5)/2).*sqrt(invTred)+((pi^2)/4 + 2).*invTred+(pi^1.5).*(invTred.^1.5);
Zrot = Zrot0.*(F./F_298);

% SPECIFIC HEAT AT CONSTANT VOLUME ROTATIONAL [-]
for i = 1:numel(Mi)
    OMEGA_diff(i)= interp1(KT_eps_bird,OMEGA_Diff,max([real(KT_eps(i)) min(KT_eps_bird)]));
    if GeoType(i)     == 0
        cv_rot(i) = 0;
    elseif GeoType(i) == 1
        cv_rot(i) = 1;          % This is Cv/R
    elseif GeoType(i) == 2
        cv_rot(i) = 3/2;        % This is Cv/R
    else
        error('error in type')
    end
end

% SPECIFIC HEAT AT CONSTANT PRESSURE [-]
cp = cpi./R;
% SPECIFIC HEAT AT CONSTANT VOLUME []
cv = cp - 1;
% SPECIFIC HEAT AT CONSTANT VOLUME VIBRATIONAL [-]
cv_vib = cv - cv_trans - cv_rot;

% SELF DIFFUSION COEFFICIENTS [m2/s]
% Dii = 3/16*sqrt((2*pi*kb^3.*T.^3)./(Mi./NA))./(P*pi*(sigma.^2).*(OMEGA_diff));
Dii = 3/16*sqrt(2*pi*kb^3*(T^3)./(Mi./2./NA))./(P*pi*sigma.^2.*OMEGA_diff);

% VIBRATIONAL CONTRIBUTIONS [-]
fvib = ((P.*Mi./(R0.*T)).*Dii)./mui;

% ROTATIONAL CONTRIBUTIONS  [-]
A = 5/2-fvib;
B = Zrot + 2/pi*(5/3*cv_rot + fvib);
C = (2/pi)*(A./B);
frot = fvib.*(1+C);

% TRANSLATIONAL CONTRIBUTIONS [-]
ftrans = 5/2*(1-C.*cv_rot./cv_trans);

% SINGLE SPECIES CONDUCTIVITIES [W/mK]
lambdai = (mui./(Mi./R0)).*(ftrans.*cv_trans+frot.*cv_rot+fvib.*cv_vib);

% MIXTURE-AVERAGED THERMAL CONDUCTIVITY [W/mK]
term1 = sum(Xi.*lambdai);
term2 = 1./(sum(Xi./lambdai));
lambda = 0.5*(term1 + term2);
% MIXTURE THERMAL CONDUCTIVITY WITH WASSILJEWA MODEL
%[lambda] = lambdaMixWassiljewa(mui,lambdai,Mi*1e+03,Xi);

end

