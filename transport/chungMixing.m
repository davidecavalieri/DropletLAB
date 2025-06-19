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

function  [Tc_m,Vc_m,omega_m,eps_m,M_m,Fcm] = chungMixing(Species,X)

%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% chungMixing performs the non-ideal mixing rules of the Chung's Model

% Inputs:
% 1) Species     : Struct variable containing the Species inputs parameters
% 2) X  : Mole Fractions

% Outputs:
% 1) Tc_m     : Mixture critical temperature
% 2) Vc_m     : Mixture critical volume
% 3) omega_m  : Mixture acentric factor
% 4) eps_m    : Mixture Lennard-Jones potential
% 5) M_m      : Mixture molecular weight
% 6) Fc_m     : Mixture

% Ref
% (i) Generalized Multiparameter Correlation for Nonpolar and Polar Fluid
% Transport Properties. Ind. Eng. Chem. Res. 1988,27, 671-679

Ns   = numel(X);
EXP1 = 1/3;
EXP2 = 2/3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma3m = 0.;
omegaNm = 0.;
epsNm   = 0;
MNm    = 0.;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMBINING RULES %%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:Ns
    for j=1:Ns
        if i == j
            sigma_ij   = 0.809*((Species(i).Vc)^EXP1);
            eps_ij   = Species(i).Tc/1.2593;
            omega_ij = Species(i).omega;
            M_ij     = Species(i).M;
        else
            %%%%%%%%%%%%%%%%%%%%% ELEMENTI OFF-DIAGONAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            sigma_ij = sqrt((0.809*((Species(i).Vc)^EXP1))*(0.809*((Species(j).Vc)^EXP1)));
            eps_ij   = sqrt((Species(i).Tc/1.2593)*(Species(j).Tc/1.2593));
            omega_ij = (Species(i).omega + Species(j).omega)/2;
            M_ij     = 2*(Species(i).M*Species(j).M)/(Species(i).M + Species(j).M );
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end

        sigma3m   = sigma3m  + X(i) * X(j) * (sigma_ij^3);
        epsNm     = epsNm    + X(i) * X(j) * (eps_ij) * (sigma_ij^3);
        omegaNm   = omegaNm  + X(i) * X(j) * (omega_ij) * (sigma_ij^3);
        MNm       = MNm      + X(i) * X(j) * (eps_ij) * (sigma_ij^2) * (sqrt(M_ij));

    end

end

%% OUTPUTS-----------------------------------------------------------------
Vc_m    = sigma3m/(0.809^3);
eps_m   = epsNm/sigma3m;
omega_m = omegaNm/sigma3m;
M_m     = ((MNm)/((sigma3m^EXP2)*eps_m))^2;
Fcm     = 1 - 0.2756*omega_m;                                              %OK
Tc_m    = 1.2593*eps_m;

end