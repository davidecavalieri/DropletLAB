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

function [rho,cp,cpvap_f,mu,lambda,Di_f,Pr_f,Sci_f,Lei_f,Sc_f_bulk] = ...
    getFilmProperties(Droplet,Ambient,EoS,ViscosityModel,ConductivityModel, ...
    DiffusionModel,EffectiveDiffusionModel,Pamb,T,Yfi,Yfai)

%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inputs:
% 1) Droplet             : Struct variable containing the droplet inputs parameters
% 2) Ambient             : Struct variable containing the ambient inputs parameters
% 3) EoS                 : String containing the selected EoS model
% 4) ViscosityModel      : String containing the selected viscosity model
% 5) ConductivityModel   : String containing the selected conductivity model
% 6) DiffusionModel      : String containing the selected diffusion model
% 7) T                   : Film temperature [K]
% 8) Yfi                 : Film mass fraction (droplet species) [-]
% 9) Yfai                : Film mass fraction (ambient species) [-]

% Outputs:
% 1) rho                 : Mixture density [kg/m^3]
% 2) cp                  : Mixture Isobaric specific heat [J/kgK]
% 3) cpvap_f             : Isobaric specific heat (droplet species) [J/kgK]
% 4) mu                  : Mixture viscosity [Pa*s]
% 5) k                   : Mixture thermal conductivity [W/mK]
% 6) Di_f                : Mixture mass diffusion coefficient [m^2/s]
% 7) Pr_f                : Prandtl number
% 8) Sci_f               : Schmidt numbers
% 9) Lei_f               : Lewis numbers

R0 = 8.314472;

%% INPUTS REALLOCATIONS

% Molecular weigth array
MiD = cell2mat({Droplet.M});
MiA = cell2mat({Ambient.M});
Mi  = [MiD MiA];

% FILM MASS FRACTIONS (DROPLET+AMBIENT SPECIES)
Yi  = [Yfi Yfai];

% FILM MOLAR FRACTIONS AND AVERAGE MIXTURE WEIGHT
[Xi,Mf] = Mass2Mol(Mi,Yi);


%% THERMO-CALORIC PROPERTIES -------------------------------------------
if strcmp(EoS,'IdealGas')

    % IDEAL GAS MIXTURE DENSITY
    rho = Pamb*Mf/(R0*T);

    % IDEAL GAS MIXTURE ISOBARIC SPECIFIC HEAT (NASA7)
    Species = [Droplet Ambient];
    [cpi,cp,~,~,~,~,~,~] = YawsPoly_Gas(T,Xi,Yi,Species);

    % IDEAL GAS ISOBARIC SPECIFIC HEAT OF DROPLET SPECIES! (see Heat-correction)
    cpvap_f = cpi(1:numel(MiD));

elseif strcmp(EoS,'PR') || strcmp(EoS,'SRK') || strcmp(EoS,'RKPR')

    Species = [Droplet Ambient];
    % REAL-FLUID EOS MIXTURE DENSITY AND SPECIFIC HEAT
    [rho,cp,~,~] = RFproperties(Pamb,T,Yi,Species,EoS);

    % REAL-FLUID EOS DENSITY AND SPECIFIC HEAT OF DROPLET SPECIES
    Composition = zeros(1, numel(Species));
    cpi     = zeros(1, numel(Species));
    for i = 1:numel(Species)
        Composition(i) = 1;
        [~,cpi(i),~,~] = RFproperties(Pamb,T,Composition,Species,EoS);
        Composition(i) = 0;
    end
    cpvap_f = cpi(1:numel(MiD));

else
    error('The selected EoS is not available')
end

%% VISCOSITY [Pa*s]-----------------------------------------------------
% Note: The Yaws database does not contain the coefficients for computing
% transport properties for inorganic species.

if strcmp(ViscosityModel,'Chung') % Depends on the selected EoS

    [~,~,mu,~] = RFproperties(Pamb,T,Yi,Species,EoS);

elseif strcmp(ViscosityModel,'Sutherland')

    A1 = 1457e-09; A2 = 110; k = 3/2;

    mu = (A1*(T^k))/(T + A2);

elseif strcmp(ViscosityModel,'CHEMKIN')

    % [mu,~] = calcMuCHEMKIN(Xi,T,Species);
    [mu,~] = transport_mu_k_CHEMKIN(T,Xi,cpi,Species);

else

    error('The selected viscosity model is not available')

end

%% THERMAL CONDUCTIVITY [W/mK]------------------------------------------

if strcmp(ConductivityModel,'Chung') % Depends on the selected EoS

    [~,~,~,lambda] = RFproperties(Pamb,T,Yi,Species,EoS);

elseif strcmp(ConductivityModel,'Sutherland')

    A1 = 252e-05; A2 = 200; k = 3/2;
    lambda = (A1*(T^k))/(T + A2);

elseif strcmp(ConductivityModel,'CHEMKIN')

    [~,lambda] = transport_mu_k_CHEMKIN(T,Xi,cpi,Species);

else

    error('The selected Conductivity model is not available')

end

%% DIFFUSION COEFFICIENT [m^2/s]-------------------------------------------

if strcmp(DiffusionModel,'Fuller')

    [Dij] = Fuller(Species,Pamb, T);

    % BINARY SISTEM:
    if numel(Yi) == 2
        Di_f = Dij(1,2);
    else
        % MIXTURE AVERAGE
        if strcmp(EffectiveDiffusionModel,'HC')
            [Di_f] = MixAverageHC(Xi,Dij,Yfi);
        elseif strcmp(EffectiveDiffusionModel,'Wilke')
            [Di_f] = MixAverageWilke(Xi,Dij,Yfi);
        end

    end

elseif strcmp(DiffusionModel,'FullerRiazi')

    rho0 = Pamb*Mf/(R0*T);
    [rho,~,mu,~] = RFproperties(Pamb,T,Yi,Species,EoS);
    [mu0,~]  = calcMuCHEMKIN(Xi,T,Species);
    [Dij] = FullerRiazi(Species,Pamb, T, Xi,rho,rho0,mu,mu0);

    % BINARY SISTEM:
    if numel(Yi) == 2
        Di_f = Dij(1,2);
    else
        % MIXTURE AVERAGE
        if strcmp(EffectiveDiffusionModel,'HC')
            [Di_f] = MixAverageHC(Xi,Dij,Yfi);
        elseif strcmp(EffectiveDiffusionModel,'Wilke')
            [Di_f] = MixAverageWilke(Xi,Dij,Yfi);
        end

    end

else
    error('The selected diffusion model is not available')

end

%% ADIMENSIONAL NUMBERS ---------------------------------------------------

% PRANDTL NUMBER  (Film-Mixture)
Pr_f = mu*cp/lambda;

% SCHMIDT NUMBERS
Sci_f = mu./(rho.*Di_f);

% BULK SCHMIDT NUMBER
Sc_f_bulk = mu/(rho*sum(Di_f.*Yfi));

% LEWIS NUMBERS
Lei_f = lambda./(cp.*Di_f.*rho);

end
