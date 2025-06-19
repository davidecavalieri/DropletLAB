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

function [mu,lambda,mu0] = chungTransport(eps_m,Fcm,W_m,Vc_m,omega_m,Tc_m,cv0,v,T)

%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% chungTransport computes the viscosity and thermal conductivity using
% the Chung's model.

% Inputs:
% 1) eps_m    : Mixture Lennard-Jones potential
% 2) Fc_m     : Mixture
% 3) W_m      : Mixture molecular weight
% 4) Vc_m     : Mixture critical volume
% 5) omega_m  : Mixture acentric factor
% 6) Tc_m     : Mixture critical temperature
% 7) cv0      : Ideal gas isochoric specific heat
% 8) v        : Mixture molar volume
% 9) T        : Temperature

% Outputs:
% 1) mu      : Mixture viscosity [Pa*s]
% 2) lambda  : Mixture thermal conductivity [W/mK]
% 3) mu0     : Mixture low-pressure viscosity [Pa*s]

% Ref
% (i) Generalized Multiparameter Correlation for Nonpolar and Polar Fluid
% Transport Properties. Ind. Eng. Chem. Res. 1988,27, 671-679

%%%%%%%%%%%%%%%%%%%% DIMENSIONLESS TEMPERATURE [-] %%%%%%%%%%%%%%%%%%%%%%%%
Tstar_m = T/eps_m;
%%%%%%%%%%%%%%%%%%%%%%% COLLISION INTEGRAL [-] %%%%%%%%%%%%%%%%%%%%%%%%%OK!
[sigmaV] = collision_integral(Tstar_m);
%%%%%%%%%%%%%%%%%%%%%%% LOW-PRESSURE VISCOSITY [Pois] %%%%%%%%%%%%%%%%%%%%%
EXP=2/3;
mu0 = (4.0785*10^-5) * (sqrt(W_m*T) * Fcm ) /(sigmaV * Vc_m^EXP); % [P]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Y %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = (Vc_m)/(6*v*1e6);
G1 = (1 - 0.5*y)/((1-y)^3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHUNG'S COEFFICIENTS %%%%%%%%%%%%%%%%%%%%%%%%%%
a1 = [6.32402 0.0012102 5.28346 6.62263 19.74540 -1.89992 24.27450 0.79716 -0.23816 0.068629]; % i=1->10
a2 = [50.41190 -0.0011536 254.209 38.0957 7.63034 -12.53670 3.44945 1.11764 0.067695 0.34793];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%POLYNOMIAL CONSTRUCTION %%%%%%%%%%%%%%%%%%%%%%%%%%
E = zeros(1,10);
for k=1:10
    E(k)= a1(k) + a2(k)*omega_m;                                  %[-]
end
nG2 = (E(1)*(1 - exp(-E(4)*y)))/y  + E(2)*G1*exp(E(5)*y) + E(3)*G1;
dG2 = (E(1)*E(4) + E(2) + E(3));
G2  =  nG2/dG2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
muK = mu0 * ( 1/G2 + E(6)*y); % [Pois]
muP = ((36.344*10^-6)*((sqrt(Tc_m*W_m))/(Vc_m^EXP))) * E(7)*y*y*G2*(exp(E(8) + E(9)*(1/Tstar_m) + E(10)*(Tstar_m^-2)));
%%%%%%%%%%%%%%%%%%%%%% High-PRESSURE VISCOSITY [Pa*s] %%%%%%%%%%%%%%%%%%%%%
mu = 0.1*(muK + muP);

%%%%%%%%%%%%%%%%%%%%%%%%% THERMAL CONDUCTIVITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tr = T/Tc_m;
Rcal = 1.987; % [cal/mol*K]
cv0_mol = cv0*W_m*1e-03;
cv0_mol_cal = cv0_mol*0.239006;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = cv0_mol_cal/Rcal - 1.5;
beta = 0.7862  - 0.7109*omega_m   + 1.3168*omega_m*omega_m;
gamma = 2.0  + 10.5*Tr*Tr;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psiN = 0.215 + 0.28288*alpha - 1.061*beta + 0.26665*gamma;
psiD = 0.6366 + beta*gamma + 1.061*alpha*beta;
psi  = 1 + alpha*(psiN/psiD);
%%%%%%%%%%%%%%%%% LOW-PRESSURE CONDUCTIVITY [cal/cm*s*K] %%%%%%%%%%%%%%%%%%
lambda0 = 7.452 * (mu0/W_m) * psi; % cal/(cm*s*K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z1 = [2.41657 -0.50924 6.61069 14.54250 0.79274 -5.8634 81.171];
z2 = [0.74824 -1.50936 5.62073 -8.91387 0.82019 12.8005 114.1580];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B = zeros(1,7);
for i=1:7
    B(i) = z1(i) +z2(i)*omega_m;
end
nH2 = (B(1)*(1 - exp(-B(4)*y)))/y  + B(2)*G1*exp(B(5)*y)   +   B(3)*G1;
dH2 = B(1)*B(4) + B(2) + B(3);
H2  = nH2/dH2;
%%%%%%%%%%%%%%%%%%%%%% LOW-PRESSURE CONDUCTIVITY [W/m*K] %%%%%%%%%%%%%%%%%%
lambdak = lambda0*(1/H2  + B(6)*y);
lambdap = ((3.039*10^-4) * ((sqrt(Tc_m/W_m)) /(Vc_m^EXP))) * B(7)*y*y*H2*sqrt(Tr);
lambda = lambdak + lambdap;
lambda = (lambda/0.239)*100;

end
