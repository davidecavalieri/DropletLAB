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

function  [Droplet, Ambient] = databaseReader(Droplet,Ambient)

%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The databaseReader function performs a pre-process phase in which a struct
% variable is created that contains thermodynamic, critical and polynomial
% coefficient properties for both Droplet and Ambient species.

% Inputs:
% 1) Droplet: Struct. variable containing the droplet species names
% 2) Ambient: Struct. variable containing the ambient species names

% Outputs:
% -------Droplet and Ambient: Global struct. variable containing:-----------
% Tc    = Critical temperatures.
% Pc    = Critical pressures.
% omega = Acentric factors.
% M     = Molecular weights.
% Vc    = Critical volumes.
% Zc    = Critical compressibility factors.
% vd    = Sum of the atomic diffusion volumes.
% Arho  = A-coefficient for liquid density computation.
% Brho  = B-coefficient for liquid density computation.
% n     = n-coefficient for liquid density computation.
% Acp   = A-coefficient for liquid isobaric specific heat computation.
% Bcp   = B-coefficient for liquid isobaric specific heat computation
% Ccp   = C-coefficient for liquid isobaric specific heat computation.
% Dcp   = D-coefficient for liquid isobaric specific heat computation.
% Amu   = A-coefficient for liquid viscosity computation.
% Bmu   = B-coefficient for liquid viscosity computation.
% Cmu   = C-coefficient for liquid viscosity computation.
% Dmu   = D-coefficient for liquid viscosity computation.
% Ak    = A-coefficient for liquid thermal-conductivity computation.
% Bk    = B-coefficient for liquid thermal-conductivity computation.
% Ck    = C-coefficient for liquid thermal-conductivity computation.
% Ap    = A-coefficient for saturation pressure computation.
% Bp    = B-coefficient for saturation pressure computation.
% Cp    = C-coefficient for saturation pressure computation.
% Dp    = D-coefficient for saturation pressure computation.
% Ep    = E-coefficient for saturation pressure computation.
% Ah    = A-coefficient for vaporization enthalpy computation.
% nh    = B-coefficient for vaporization enthalpy computation.
% Acp_g = A-coefficient for gas isobaric specific heat computation.
% Bcp_g = B-coefficient for gas isobaric specific heat computation.
% Ccp_g = C-coefficient for gas isobaric specific heat computation.
% Dcp_g = D-coefficient for gas isobaric specific heat computation.
% Ecp_g = E-coefficient for gas isobaric specific heat computation.
% Amu_g = A-coefficient for gas viscosity computation.
% Bmu_g = B-coefficient for gas viscosity computation.
% Cmu_g = C-coefficient for gas viscosity computation.
% Ak_g  = A-coefficient for gas thermal-conductivity computation.
% Bk_g  = B-coefficient for gas thermal-conductivity computation.
% Ck_g  = C-coefficient for gas thermal-conductivity computation.
% Ad    = A-coefficient for diffusion coefficient computation.
% Bd    = B-coefficient for diffusion coefficient computation.
% Cd    = C-coefficient for diffusion coefficient computation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R0 = 8.314472;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Droplet_fluids   = {Droplet.Species};
Ambient_fluids   = {Ambient.Species};

%% LOAD INPUTS DATA (DROPLET AND AMBIENT): [Names Tc Pc omega M Vc vd]
data = readtable('database.txt', 'VariableNamingRule', 'preserve');
% Extract the Species names
species_names = data.Var1;
% Extract the critical temperature [K]
Tc            = data.Var2;
% Extract the critical pressure [Pa]
Pc            = (data.Var3).*1E6;
% Extract the acentric factor [-]
omega         = data.Var4;
% Extract the molecular weigth [Kg/mol]
M             = data.Var5.*1E-3;
% Extract the critical volume [m^3/mol]
Vc            = data.Var6;
% Extract the SUM of atomic diffusion volumes [-]
vd          =  data.Var7;

%%%%%%%%%%%%%%%%%%%%%%%%%% DROPLET SPECIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:numel(Droplet_fluids)
    speciesD_index = find(strcmp(species_names, Droplet_fluids(i)));

    % Droplet(i).Names = species_names(speciesD_index);
    Droplet(i).Tc    = Tc(speciesD_index);
    Droplet(i).Pc    = Pc(speciesD_index);
    Droplet(i).omega = omega(speciesD_index);
    Droplet(i).M     = M(speciesD_index);
    Droplet(i).Vc    = Vc(speciesD_index);
    Droplet(i).Zc    = Droplet(i).Pc*Droplet(i).Vc/(R0*Droplet(i).Tc);
    Droplet(i).vd    = vd(speciesD_index);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% AMBIENT SPECIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:numel(Ambient_fluids)
    speciesA_index = find(strcmp(species_names, Ambient_fluids(i)));

    Ambient(i).Tc    = Tc(speciesA_index);
    Ambient(i).Pc    = Pc(speciesA_index);
    Ambient(i).omega = omega(speciesA_index);
    Ambient(i).M     = (M(speciesA_index));
    Ambient(i).Vc    = Vc(speciesA_index);
    Ambient(i).Zc    = Ambient(i).Pc*Ambient(i).Vc/(R0*Ambient(i).Tc);
    Ambient(i).vd    = vd(speciesA_index);
end

%% LOAD LIQUID DATABASE (YAWS): [Arho Brho n Acp Bcp Ccp Dcp]
% Notes:
% (i) For lumped models, the liquid-phase properties needed are density
% and isoabaric specific heat!
% (ii) The only species in liquid phase are those in the droplet, for this
% reason the coefficients for the ambient species are sets to zero.

Yaws_l = readtable('Yaws_liquid.txt', 'VariableNamingRule', 'preserve');
% Extract the Species names
Names_Yaws = Yaws_l.ch_Species;
Arho       = Yaws_l.Arho;
Brho       = Yaws_l.Brho;
n          = Yaws_l.n;
Acp        = Yaws_l.Acp;
Bcp        = Yaws_l.Bcp;
Ccp        = Yaws_l.Ccp;
Dcp        = Yaws_l.Dcp;

%%%%%%%%%%%%%%%%%%%%%%%%%%% DROPLET SPECIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:numel(Droplet_fluids)
    Yaws_l_index = find(strcmp(Names_Yaws, Droplet_fluids(i)));
    Droplet(i).Arho  = Arho(Yaws_l_index);
    Droplet(i).Brho  = Brho(Yaws_l_index);
    Droplet(i).n     = n(Yaws_l_index);
    Droplet(i).Acp   = Acp(Yaws_l_index);
    Droplet(i).Bcp   = Bcp(Yaws_l_index);
    Droplet(i).Ccp   = Ccp(Yaws_l_index);
    Droplet(i).Dcp   = Dcp(Yaws_l_index);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% AMBIENT SPECIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:numel(Ambient_fluids)

    Ambient(i).Arho  = 0;
    Ambient(i).Brho  = 0;
    Ambient(i).n     = 0;
    Ambient(i).Acp   = 0;
    Ambient(i).Bcp   = 0;
    Ambient(i).Ccp   = 0;
    Ambient(i).Dcp   = 0;
end

%% LOAD SATURATION DATABASE (YAWS): [Ap Bp Cp Dp Ep Ah nh] (DROPLET)
% Note: See the previous notes!

Yaws_sat = readtable('Yaws_sat.txt', 'VariableNamingRule', 'preserve');
% Extract the Species names
Names_Yaws_sat = Yaws_sat.ch_Species;
Ap       = Yaws_sat.Ap;
Bp       = Yaws_sat.Bp;
Cp       = Yaws_sat.Cp;
Dp       = Yaws_sat.Dp;
Ep       = Yaws_sat.Ep;
Ah       = Yaws_sat.Ah;
nh       = Yaws_sat.nh;

%%%%%%%%%%%%%%%%%%%%%%%%%%% DROPLET SPECIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:numel(Droplet_fluids)
    Yaws_sat_index = find(strcmp(Names_Yaws_sat, Droplet_fluids(i)));
    Droplet(i).Ap  = Ap(Yaws_sat_index);
    Droplet(i).Bp  = Bp(Yaws_sat_index);
    Droplet(i).Cp  = Cp(Yaws_sat_index);
    Droplet(i).Dp  = Dp(Yaws_sat_index);
    Droplet(i).Ep  = Ep(Yaws_sat_index);
    Droplet(i).Ah  = Ah(Yaws_sat_index);
    Droplet(i).nh  = nh(Yaws_sat_index);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% AMBIENT SPECIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:numel(Ambient_fluids)

    Ambient(i).Ap  = 0;
    Ambient(i).Bp  = 0;
    Ambient(i).Cp  = 0;
    Ambient(i).Dp  = 0;
    Ambient(i).Ep  = 0;
    Ambient(i).Ah  = 0;
    Ambient(i).nh  = 0;
end

%% LOAD GAS-PHASE DATABASE (YAWS): [Acp_g Bcp_g Ccp_g Dcp_g Ecp_g Amu_g Bmu_g Cmu_g Ak_g Bk_g Ck_g]
% Note:
Yaws_gas = readtable('Yaws_gas.txt', 'VariableNamingRule', 'preserve');
% Extract the Species names
Names_Yaws_gas = Yaws_gas.ch_Species;
Acp_g  = Yaws_gas.Acp_g;
Bcp_g  = Yaws_gas.Bcp_g;
Ccp_g  = Yaws_gas.Ccp_g;
Dcp_g  = Yaws_gas.Dcp_g;
Ecp_g  = Yaws_gas.Ecp_g;
AcpH_g = Yaws_gas.AcpH_g;
BcpH_g = Yaws_gas.BcpH_g;
CcpH_g = Yaws_gas.CcpH_g;
DcpH_g = Yaws_gas.DcpH_g;
EcpH_g = Yaws_gas.EcpH_g;
Tmin   = Yaws_gas.Tmin;
Tmed   = Yaws_gas.Tmed;
Tmax   = Yaws_gas.Tmax;
Amu_g  = Yaws_gas.Amu_g;
Bmu_g  = Yaws_gas.Bmu_g;
Cmu_g  = Yaws_gas.Cmu_g;
Ak_g   = Yaws_gas.Ak_g;
Bk_g   = Yaws_gas.Bk_g;
Ck_g   = Yaws_gas.Ck_g;
Ad     = Yaws_gas.Ad;
Bd     = Yaws_gas.Bd;
Cd     = Yaws_gas.Cd;

%%%%%%%%%%%%%%%%%%%%%%%%%%% DROPLET SPECIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:numel(Droplet_fluids)
    Yaws_gasD_index = find(strcmp(Names_Yaws_gas, Droplet_fluids(i)));

    Droplet(i).Acp_g  = Acp_g(Yaws_gasD_index);
    Droplet(i).Bcp_g  = Bcp_g(Yaws_gasD_index);
    Droplet(i).Ccp_g  = Ccp_g(Yaws_gasD_index);
    Droplet(i).Dcp_g  = Dcp_g(Yaws_gasD_index);
    Droplet(i).Ecp_g  = Ecp_g(Yaws_gasD_index);
    Droplet(i).AcpH_g = AcpH_g(Yaws_gasD_index);
    Droplet(i).BcpH_g = BcpH_g(Yaws_gasD_index);
    Droplet(i).CcpH_g = CcpH_g(Yaws_gasD_index);
    Droplet(i).DcpH_g = DcpH_g(Yaws_gasD_index);
    Droplet(i).EcpH_g = EcpH_g(Yaws_gasD_index);

    Droplet(i).Tmin   = Tmin(Yaws_gasD_index);
    Droplet(i).Tmed   = Tmed(Yaws_gasD_index);
    Droplet(i).Tmax   = Tmax(Yaws_gasD_index);

    Droplet(i).Amu_g  = Amu_g(Yaws_gasD_index);
    Droplet(i).Bmu_g  = Bmu_g(Yaws_gasD_index);
    Droplet(i).Cmu_g  = Cmu_g(Yaws_gasD_index);

    Droplet(i).Ak_g   = Ak_g(Yaws_gasD_index);
    Droplet(i).Bk_g   = Bk_g(Yaws_gasD_index);
    Droplet(i).Ck_g   = Ck_g(Yaws_gasD_index);

    Droplet(i).Ad     = Ad(Yaws_gasD_index);
    Droplet(i).Bd     = Bd(Yaws_gasD_index);
    Droplet(i).Cd     = Cd(Yaws_gasD_index);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% AMBIENT SPECIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:numel(Ambient_fluids)
    Yaws_gasA_index = find(strcmp(Names_Yaws_gas, Ambient_fluids(i)));

    Ambient(i).Acp_g  = Acp_g(Yaws_gasA_index);
    Ambient(i).Bcp_g  = Bcp_g(Yaws_gasA_index);
    Ambient(i).Ccp_g  = Ccp_g(Yaws_gasA_index);
    Ambient(i).Dcp_g  = Dcp_g(Yaws_gasA_index);
    Ambient(i).Ecp_g  = Ecp_g(Yaws_gasA_index);
    Ambient(i).AcpH_g = AcpH_g(Yaws_gasA_index);
    Ambient(i).BcpH_g = BcpH_g(Yaws_gasA_index);
    Ambient(i).CcpH_g = CcpH_g(Yaws_gasA_index);
    Ambient(i).DcpH_g = DcpH_g(Yaws_gasA_index);
    Ambient(i).EcpH_g = EcpH_g(Yaws_gasA_index);

    Ambient(i).Tmin   = Tmin(Yaws_gasA_index);
    Ambient(i).Tmed   = Tmed(Yaws_gasA_index);
    Ambient(i).Tmax   = Tmax(Yaws_gasA_index);

    Ambient(i).Amu_g  = Amu_g(Yaws_gasA_index);
    Ambient(i).Bmu_g  = Bmu_g(Yaws_gasA_index);
    Ambient(i).Cmu_g  = Cmu_g(Yaws_gasA_index);

    Ambient(i).Ak_g   = Ak_g(Yaws_gasA_index);
    Ambient(i).Bk_g   = Bk_g(Yaws_gasA_index);
    Ambient(i).Ck_g   = Ck_g(Yaws_gasA_index);

    Ambient(i).Ad     = Ad(Yaws_gasA_index);
    Ambient(i).Bd     = Bd(Yaws_gasA_index);
    Ambient(i).Cd     = Cd(Yaws_gasA_index);

end

%% CHEMKIN TRANSPORT DATABASE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD CHEMKIN TRANSPORT DATABASE
trans = readtable('Chemkintrans.txt', 'VariableNamingRule', 'preserve');
% Extract the species names
species_CHEMKIN = trans.Var1;
% Extract the geometry
geometry      = trans.Var2;
% Extract the Lennard-Jones potential [K]
eps_kb        = trans.Var3;
% Extract the Lennard-Jones collision diameter [Angstroms] 1e-10m
sigma         = trans.Var4;
% Extract the dipole moment [Debye]
mu            = trans.Var5;
% extract the polarizability [Angstrom^3]
alpha         = trans.Var6;
% Extract the rotational relaxation collision number at T = 298 K
Zrot          = trans.Var7;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% DROPLET SPECIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:numel(Droplet_fluids)

    speciesCHEMKIN_index = find(strcmp(species_CHEMKIN, Droplet_fluids(i)));

    Droplet(i).geometry = geometry(speciesCHEMKIN_index);
    Droplet(i).eps_kb   = eps_kb(speciesCHEMKIN_index);
    Droplet(i).sigma    = sigma(speciesCHEMKIN_index);
    Droplet(i).mu       = mu(speciesCHEMKIN_index);
    Droplet(i).alpha    = alpha(speciesCHEMKIN_index);
    Droplet(i).Zrot     = Zrot(speciesCHEMKIN_index);

end

%%%%%%%%%%%%%%%%%%%%%%%%%% AMBIENT SPECIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:numel(Ambient_fluids)

    speciesCHEMKIN_index = find(strcmp(species_CHEMKIN, Ambient_fluids(i)));

    Ambient(i).geometry = geometry(speciesCHEMKIN_index);
    Ambient(i).eps_kb   = eps_kb(speciesCHEMKIN_index);
    Ambient(i).sigma    = sigma(speciesCHEMKIN_index);
    Ambient(i).mu       = mu(speciesCHEMKIN_index);
    Ambient(i).alpha    = alpha(speciesCHEMKIN_index);
    Ambient(i).Zrot     = Zrot(speciesCHEMKIN_index);

end

end


