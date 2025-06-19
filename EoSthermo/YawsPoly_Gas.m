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

function [cpi,cp,mui,mu,lambdai,lambda,Di,D] = YawsPoly_Gas(T,Xi,Yi,Species)

%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YawsPoly_Gas calculates gas-phase properties using polynomial models

% Inputs:
% 1) T  : Temperature   [K]
% 2) Xi : Mole fractions 
% 3) Yi : Mass fractions
% 4) Species     : Struct variable containing the Species inputs parameters

% Outputs
% 1) cpi    : Gas Isobaric-Specific-Heat of each components [J/kgK]
% 2) cp     : Gas Mixture Isobaric-Specific-Heat            [J/kgK]
% 3) mui    : Gas viscosity of each components              [Pa*s]
% 4) mu     : Gas Mixture viscosity                         [Pa*s]
% 5) lambai : Gas conductivity of each components           [W/mK]
% 6) lambda : Gas Mixture conductivity                      [W/mK]
% 7) Di     : Gas diffusion coefficients of each components [m^2/s]
% 8) D      : Gas diffusion coefficient                     [m^2/s]

% REFs:
% (i) The Yaws Handbook of Physical Properties for Hydrocarbons and Chemicals: 
% Physical Properties for More Than 54,000 Organic and Inorganic Chemical 
% Compounds,Coverage for C1 to C100 Organics and Ac to Zr Inorganics, 
% Gulf Professional Publishing (2015) 

% (ii) B. E. Poling, J. M. Prausnitz, J. P. OConnell. Properties of Gases 
% and Liquids, 5th Edition, 2001, McGraw-Hill Education 

R0 = 8.314472;

% Allocations
cpi      = zeros(1,length(Yi));
mui      = zeros(1,length(Yi));
lambdai  = zeros(1,length(Yi));
Di       = zeros(1,length(Yi));

for i = 1:length(Yi)
    
    if T<=Species(i).Tmed
    cpi(i)     = R0*(Species(i).Acp_g + (Species(i).Bcp_g)*T + (Species(i).Ccp_g)*T*T ...
                 + (Species(i).Dcp_g)*T*T*T + (Species(i).Ecp_g)*T*T*T*T)/(Species(i).M);
    else
    cpi(i)     = R0*(Species(i).AcpH_g + (Species(i).BcpH_g)*T + (Species(i).CcpH_g)*T*T ...
                 + (Species(i).DcpH_g)*T*T*T + (Species(i).EcpH_g)*T*T*T*T)/(Species(i).M);
    end

    mui(i)     = (Species(i).Amu_g + Species(i).Bmu_g*T + Species(i).Cmu_g*T*T)*1e-7;

    lambdai(i) = Species(i).Ak_g + Species(i).Bk_g*T + Species(i).Ck_g*T*T;

    Di(i)      = (Species(i).Ad + Species(i).Bd*T + Species(i).Cd*T*T)*1e-4;
end

    % GAS MIXTURE ISOBARIC SPECIFIC HEAT [J/kgK]
    cp      = sum(Yi.*cpi);

    % GAS MIXTURE VISCOSITY [Pa*s]
    Mi = cell2mat({Species.M});
    [mu] = muMixWilke(mui,Mi*1e+03,Xi);

    % MIXTURE-AVERAGED THERMAL CONDUCTIVITY [W/mK]
    term1 = sum(Xi.*lambdai);
    term2 = 1./(sum(Xi./lambdai));
    lambda = 0.5*(term1 + term2);

    % GAS DIFFUSION COEFFICIENT IN AIR [m^2/s] (TODO: HC)
    D       = 1/sum(Xi./Di); % Blanc's Law

end