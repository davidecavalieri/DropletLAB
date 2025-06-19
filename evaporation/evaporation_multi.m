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

function [dw, d, mdot_i, Qg, Ql, Bm, Xvsi, Yvs, Tf, Yf, rho_f, cpvap_f, cp_f, Le_f, ...
    f, D_f, mu_f, lambda_f, Lv, Sh_fAS, Nu_fAS] = evaporation_multi(t, w, activeComponents, ...
    Droplet, Ambient, liquidRho, liquidCp, SatPressure, vapEnthalpy, Pamb, ...
    Tamb, VLE, Ar, EoS, DiffusionModel, EffectiveDiffusionModel, StefanFlow, AbramzonSirignano, NaturalConvection, ...
    Yvi_inf, XviA_inf, YviA_inf, Model, ViscosityModel, ConductivityModel, Vrel)

persistent x_prev

%% STATE VECTOR VARIABLES--------------------------------------------------
% DROPLET TEMPERATURE [K]
Td = w(1);
Td = real(Td);

% DROPLET MASS OF EACH COMPONENTS [kg]
w(2:end) = max(real(w(2:end)), 0);
mdi = w(2:end)';

%% LIQUID PHASE -----------------------------------------------------------
% COMPUTE THE LIQUID PHASE COMPOSITION
[md, Yli, Xli, ~] = getLiquidComposition(Droplet, mdi);
Xli = max(Xli, 0);
Yli = max(Yli, 0);

% COMPUTE DROPLET LIQUID PROPERTIES
[~, rhoL, ~, cpL] = getLiquidProperties(Droplet, liquidRho, liquidCp, Pamb, Td, Yli, EoS);

% COMPUTE DROPLET SIZE
[~, r, d, ~] = getSize(md, rhoL);

%% VAPOUR-LIQUID-EQUILIBRIUM (Td = Tsurface)
% COMPUTE SATURATION PRESSURES
[Psati] = getSatPressure(Droplet, SatPressure, Td);

if strcmp(VLE,'Raoult')
    % DROPLET SURFACE MOLAR AND MASS FRACTIONS (Xvsi, Yvsi) and AMBIENT SURFACE MOLAR AND MASS FRACTIONS (Xgsi, Ygsi)
    [Xvsi, Xgsi, Yvsi, Ygsi] = RaoultLaw(Droplet,Ambient,Pamb,Xli,Psati,XviA_inf,YviA_inf);
    % VAPOR PHASE MOLE FRACTIONS @ DROPLET SURFACE
    Xsi = [Xvsi Xgsi];
    % VAPOR PHASE MASS FRACTIONS @ DROPLET SURFACE
    Ysi = [Yvsi Ygsi];
    % COMPUTE VAPORIZATION ENTHALPIES
    [Lvi] = getVapEnthalpy(Droplet, Ambient, vapEnthalpy, Xli, Xsi, Td, Pamb);
elseif strcmp(VLE,'VLE-EoS')
    % DROPLET SURFACE MOLAR AND MASS FRACTIONS (Xvsi, Yvsi) and AMBIENT SURFACE MOLAR AND MASS FRACTIONS (Xgsi, Ygsi)
    [Xvsi, Xgsi, Yvsi, Ygsi] = VLE_EoS(Droplet,Ambient,Pamb,Xli,Td,Psati,XviA_inf,YviA_inf,EoS,activeComponents,t);
    % VAPOR PHASE MOLE FRACTIONS @ DROPLET SURFACE
    Xsi = [Xvsi Xgsi];
    Xsi = max(Xsi,0);
    % VAPOR PHASE MASS FRACTIONS @ DROPLET SURFACE
    Ysi = [Yvsi Ygsi];
    Ysi = max(Ysi,0);
    % COMPUTE VAPORIZATION ENTHALPIES
    [Lvi] = getVapEnthalpy(Droplet, Ambient, vapEnthalpy, Xli, Xsi, Td, Pamb);
elseif strcmp(VLE,'VLE-EoS-solub')
    % DROPLET SURFACE MOLAR AND MASS FRACTIONS (Xvsi, Yvsi) and AMBIENT SURFACE MOLAR AND MASS FRACTIONS (Xgsi, Ygsi)
    [Xlsi,Xsi,Xvsi,~,Yvsi,Ygsi] = VLE_EoS_solubility(Droplet,Ambient,Pamb,Xli,Td,Psati,XviA_inf,YviA_inf,EoS);
    % VAPOR PHASE MASS FRACTIONS @ DROPLET SURFACE
    Ysi = [Yvsi Ygsi];
    % COMPUTE VAPORIZATION ENTHALPIES
    [Lvi] = getVapEnthalpy(Droplet, Ambient, vapEnthalpy, Xlsi, Xsi, Td, Pamb);
else
    error('The selected VLE model is not available');
end

%% GAS PHASE --------------------------------------------------------------
% COMPUTE THE FILM COMPOSITION
[Tf, Yf, Yfa] = getFilmConditions(Tamb, Ar, Td, Ysi, Yvi_inf, YviA_inf);

% COMPUTE THE FILM PROPERTIES
[rho_f, cp_f, cpvap_f, mu_f, lambda_f, D_f, Pr_f, Sc_f, Le_f, Sc_f_bulk] = ...
    getFilmProperties(Droplet, Ambient, EoS, ViscosityModel, ConductivityModel, ...
    DiffusionModel, EffectiveDiffusionModel, Pamb, Tf, Yf, Yfa);

% RANZ-MARSHALL
[Sh_f0, Sh_f0_bulk, Nu_f0] = getRanzMarshall(Ambient, NaturalConvection, d, Vrel, Tamb, Td, Sc_f, Sc_f_bulk, Pr_f,Pamb,YviA_inf);

%% MODELS FOR SINGLE-COMPONENT DROPLET-------------------------------------
if strcmp(Model, 'CEM-B')

    tau_d = 4 * rhoL * (r^2) / (18 * mu_f);

    % CUMULATIVE SPALDING MASS TRANSFER NUMBER [-]
    Bm = (sum(Yvsi) - sum(Yvi_inf)) / (1 - sum(Yvsi));

    % SPECIFIC HEAT OF A PSEUDO-FUEL SPECIES AT FILM CONDITIONS
    cpvap_f_mix = sum(Yf./sum(Yf).*cpvap_f);

    % ABRAMZON-SIRIGNANO HEAT AND MASS TRANSFER CORRECTIONS-------------------------------------
    if strcmp(AbramzonSirignano, 'True')

        % MODIFIED NUSSELT NUMBER AND AVERAGED SHERWOOD NUMBER
        [Sh_fAS_i,Sh_fAS,Nu_fAS] = AScorrection_multi(Bm,cpvap_f_mix,cp_f,Pr_f,Sh_f0_bulk,Nu_f0,Sh_f0,Sc_f_bulk,Sc_f);

        Sh_f0_bulk = Sh_fAS;
        Nu_f0 = Nu_fAS;
        Sh_f0 = Sh_fAS_i;

    end

    if numel(Yvsi) > 1
        % INDIVIDUAL SPALDING MASS TRANSFER NUMBERS
        eta_i = (Sh_f0_bulk.*Sc_f)./(Sh_f0.*Sc_f_bulk);
        Bm_i = (1 + Bm).^eta_i - 1;
        % NON-DIMENSIONAL PARTIAL EVAPORATION RATES
        eps_i = ((1 + Bm_i).*Yvsi - Yvi_inf) ./ Bm_i;
    else
        eps_i = 1;
    end

    % TOTAL MASS VAPORIZATION RATE
    dmdt = -(Sh_f0_bulk / (3 * Sc_f_bulk)) * (md / tau_d) * log(1 + Bm) .* eps_i;

    if strcmp(StefanFlow, 'True')

        % COMPUTE STEFAN-FLOW CORRECTION FACTOR FOR HEAT TRANSFER
        beta = 0;
        if (sum(mdi) > 0 && r > 0)
            beta = -3 * cpvap_f_mix/cp_f * Pr_f/Nu_f0 * sum(dmdt) / (md / tau_d);
        end

        f2 = 0;
        if beta ~= 0
            f2 = beta / (exp(beta) - 1);
        end

    else

        f2 = 1;

    end

    % COMPUTE TEMPERATURE DERIVATIVE
    dTdt = (Nu_f0 / (3 * Pr_f)) * (cp_f / cpL) * (f2 / tau_d) * (Tamb - Td) + ...
        dot(dmdt, Lvi) / (md * cpL);

else

    error('The selected evaporation model is not available');

end

%% ODE---------------------------------------------------------------------
dw = [dTdt; dmdt'];
dw = real(dw);

end