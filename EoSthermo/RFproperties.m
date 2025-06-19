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

function [rho,cp,mu,k] = RFproperties(p,T,Y,Species,EoS)

%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RFproperties calculates thermodynamic and transport properties using
% cubic EoS based model.

% Inputs:
% 1) p                   : pressure [Pa]
% 2) T                   : Temperature [K]
% 3) Y                   : Mass fraction [-]
% 4) Species             : Struct variable containing the Species inputs parameters

% Outputs:
% 1) rho                 : Mixture density [kg/m^3]
% 2) cp                  : Mixture Isobaric specific heat [J/kgK]
% 3) mu                  : Mixture viscosity [Pa*s]
% 4) k                   : Mixture thermal conductivity [W/mK]

R0 = 8.314472;

%%%%%%%%%%%%%%%%%% Mass2Mol fraction conversion %%%%%%%%%%%%%%%%%%%%%%%%%%%
[X,M] = mass2mol(Species,Y);                                    %Check: ok!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% VdW Mixing-Rules %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[am, bm, dadTm,ddadTm] = VdWmixing(Species,X,T,EoS);                %Check: ok!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Z] = getZ(Species,am,bm,p,T,X,EoS); v=(Z * R0 * T)/p;              %Check: ok!
rho_mol = 1/v; rho = rho_mol*M;

%% CALORIC PROPERTIES -----------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%% Ideal-Thermo %%%%%%%%%%%%%%%%%%%%%%%%%Check: ok!
[~,cp0,~,~,~,~,~,~] = YawsPoly_Gas(T,X,Y,Species);
cv0 = cp0 - (R0/M);
%%%%%%%%%%%%%%%%%%%%%%%%%% Departure functions %%%%%%%%%%%%%%%%%%Check: ok!
[~,~,dcV] = departureFunctions(Species,am,bm,dadTm,ddadTm,M,v,T,p,X,EoS);
%%%%%%%%% Specific heat at constant volume %%%%%%%
cv = cv0 + dcV;
%%%%%%%%%%%%%%%%%%%%%%%% Thermo-Derivatives %%%%%%%%%%%%%%%%%%%%%%Check: ok!
[dpdT_v,dpdv_T] = thermoDerivatives(Species,am,bm,dadTm,v,T,X,EoS);
cp      = cv - (T * dpdT_v^2 / dpdv_T)/M;

%% TRANSPORT PROPERTIES ---------------------------------------------------
[Tc_m,Vc_m,omega_m,eps_m,M_m,Fcm]= chungMixing(Species,X);
%%%%%%% Chung's model problem for multicomponent mixture %%% (Hydrogen)
if omega_m < 0
    omega_m = 0;
end
if Fcm > 1
    Fcm = 1;
end

[mu,k,mu0]=chungTransport(eps_m,Fcm,M_m.*1e3,Vc_m*1e6,omega_m,Tc_m,cv0,v,T); %Check: ok!

end