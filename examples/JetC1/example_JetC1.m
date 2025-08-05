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

% DESCRIPTION:
% The present example calculates the isolated droplet vaporization process
% under forced convection conditions for a two-component Jet C-1 surrogate mixture
% within air ambient at 10 bar.

%% ----- UTILITIES ----- %%

clear all;
addpath('../../Utilities/');
addpath('../../LiquidProperties');
addpath('../../EoSthermo');
addpath('../../transport');
addpath('../../VLE');
addpath('../../FilmProperties');
addpath('../../HeatCorrection');
addpath('../../evaporation');

%% ----- LIBRARIES ----- %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Model             = 'CEM-B';
%%%%%%%%%%%%%%%%%%%%%%%%% LIQUID PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
liquidRho        = 'Yaws'; % options: Yaws, Coolprop (check species availability), EoS
liquidCp        = 'Yaws'; % options: Yaws, Coolprop (check species availability), EoS
%%%%%%%%%%%%%%%%%%%%%%%%% SATURATED PRESSURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SatPressure       = 'Yaws'; % options: Yaws, Coolprop (check species availability), LeeKesler
%%%%%%%%%%%%%%%%%%%%%%%%% VAPORIZATION ENTHALPY %%%%%%%%%%%%%%%%%%%%%%%%%%%
vapEnthalpy       = 'Yaws'; % options: Yaws, Coolprop (check species availability)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VLE LAWS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VLE               = 'VLE-EoS'; % options: Raoult, VLE-EoS, VLE-EoS-solub
%%%%%%%%%%%%%%%% WEIGHT COEFFICIENT FOR THE FILM LAW %%%%%%%%%%%%%%%%%%%%%%
Ar                = 1/3;
%%%%%%%%%%%%%%%%%%%%%%%%%%% EQUATION OF STATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EoS               = 'PR'; % options: IdealGas, PR, SRK, RKPR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VISCOSITY MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ViscosityModel    = 'Chung'; % options: Chung, Sutherland, CHEMKIN
%%%%%%%%%%%%%%%%%%%%%%%%% CONDUCTIVITY MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ConductivityModel = 'Chung'; % options: Chung, Sutherland, CHEMKIN
%%%%%%%%%%%%%%%%%%%% COEFFICIENTS DIFFUSION MODEL %%%%%%%%%%%%%%%%%%%%%%%%%
DiffusionModel    = 'FullerRiazi'; % options: Fuller, FullerRiazi
%%%%%%%%%%%%%%%%%%%% EFFECTIVE DIFFUSION MODEL %%%%%%%%%%%%%%%%%%%%%%%%%
EffectiveDiffusionModel    = 'HC'; % options: HC, Wilke
%%%%%%%%%%%%%%%%%%%% STEFAN-FLOW CORRECTION %%%%%%%%%%%%%%%%%%%%%%%%%
StefanFlow    = 'True';
%%%%%%%%%%%%%%%%%%%% ABRAMZON-SIRIGNANO HEAT AND MASS TRANSFER CORRECTIONS %%%%%%%%%%%%%%%%%%%%%%%%%
AbramzonSirignano    = 'False';
%%%%%%%%%%%%%%%%%%%% NATURAL CONVECTION CORRECTION %%%%%%%%%%%%%%%%%%%%%%%%
NaturalConvection = 'on'; % 1) On; 2) Off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ----- DROPLET VAPORIZATION ----- %%

% +++ Initial conditions (ICs) +++ %
% Initial droplet diameter [m]
d0 = 30e-06; %
% Initial droplet temperature [K]
T0 = 300;

% +++ Boundary conditions (BCs) +++ %
% Ambient pressure [Pa]
Pamb = 10*1e+05;
% Ambient temperature [K]
Tamb = 800;
% Relative velocity between ambient and droplet [m/s]
Vrel = 20;
% +++ Gaseous ambient composition (far away from the droplet) +++ %
Ambient(1).Species = 'nitrogen';
Ambient(2).Species = 'oxygen';
% Molar fractions of gaseous ambient species (nitrogen, oxygen, water, carbon-monoxide, carbon-dioxide)
XviA_inf = [0.79 0.21];
% Mass fractions of gaseous ambient species (nitrogen, oxygen, water, carbon-monoxide, carbon-dioxide)
YviA_inf = [0.76 0.24];

% +++ Droplet composition +++ %
Droplet(1).Species = '2,2,4,6,6-pentamethylheptane';
Droplet(2).Species = '2,2,4,4,6,8,8-heptamethylnonane';
Xli0(1) = 0.81658;
Xli0(2) = 0.18342;
% +++ Vapor mass fractions of liquid-phase species far away from the droplet +++ %
Yvi_inf = zeros(1, numel(Droplet));

% Create a containers.Map from liquid species names to global indices
species_l_map = containers.Map({Droplet.Species}, 1:numel(Droplet));

% +++ CREATE STRUCTURE VARIABLE FOR DROPLET AND AMBIENT +++ %
[Droplet, Ambient] = databaseReader(Droplet, Ambient);

% +++ Inputs for the ODE system +++ %
[mdi0] = masses(Droplet, d0, Xli0, T0, liquidRho, Pamb, EoS);
initialConditions = [T0 mdi0];
tend = 100;

min_fraction = 1e-03; % threshold molar fraction below which one component is depleted

% +++ Solve the equations using ODE45 +++ %

% Number of components in the initial state
Ns = length(initialConditions) - 1; % Excluding droplet temperature
% Start with all components active
activeComponents = 1:Ns;
active_species_names = {Droplet.Species};

% Set event options for both component depletion and droplet shrinkage along with output function
options = odeset('Events', @(t, w) combinedEvents(t, w, d0, Droplet, Pamb, liquidRho, liquidCp, EoS, min_fraction), ...
    'OutputFcn', @(t, w, flag) store_diameter(t, w, flag, Droplet, liquidRho, liquidCp, Pamb, EoS), ...
    'RelTol', 1e-12, 'AbsTol', ones(numel(Yvi_inf) + 1, 1) * 1e-12);

% Initialize time, state, and diameter (and diameter times) vectors
t = [];
w = [];
d = [];
t_d = [];

% Initial run
tspan = [0, tend];
currentConditions = initialConditions;

while true

    % Solve ODE system
    if isempty(w)
        [t1, w1, te, we, ie] = ode45(@(t, w) evaporation_multi(t, w, activeComponents, ...
            Droplet, Ambient, liquidRho, liquidCp, SatPressure, vapEnthalpy, ...
            Pamb, Tamb, VLE, Ar, EoS, DiffusionModel, EffectiveDiffusionModel, StefanFlow, AbramzonSirignano, NaturalConvection, ...
            Yvi_inf, XviA_inf, YviA_inf, Model, ...
            ViscosityModel, ConductivityModel, Vrel), tspan, currentConditions, options);
    else
        w_red = expandStateVector(w(end,:), activeComponents);
        Droplet = Droplet(activeComponents);
        Yvi_inf = Yvi_inf(activeComponents);
        % Update active components
        activeComponents = 1:numel(Droplet);
        active_species_names = {Droplet.Species};
        % Set event options for both component depletion and droplet shrinkage
        options = odeset('Events', @(t, w_red) combinedEvents(t, w_red, d0, Droplet, Pamb, liquidRho, liquidCp, EoS, min_fraction), ...
            'OutputFcn', @(t, w_red, flag) store_diameter(t, w_red, flag, Droplet, liquidRho, liquidCp, Pamb, EoS), ...
            'RelTol', 1e-12, 'AbsTol', ones(numel(Yvi_inf) + 1, 1) * 1e-12);
        [t1, w1, te, we, ie] = ode15s(@(t, w_red) evaporation_multi(t, w_red, activeComponents, ...
            Droplet, Ambient, liquidRho, liquidCp, SatPressure, vapEnthalpy, ...
            Pamb, Tamb, VLE, Ar, EoS, DiffusionModel, EffectiveDiffusionModel, StefanFlow, AbramzonSirignano, NaturalConvection, ...
            Yvi_inf, XviA_inf, YviA_inf, Model, ...
            ViscosityModel, ConductivityModel, Vrel), tspan, currentConditions, options);
    end

    % +++ Start appending results +++ %
    t = [t; t1];

    w_frac = zeros(size(w1,1), Ns+1);
    if size(w1,2) == Ns+1
        w_frac = w1;
    else
        w_frac(:,1) = w1(:,1);
        for j = 1:length(active_species_names)
            global_idx = species_l_map(active_species_names{j});
            w_frac(:, global_idx + 1) = w1(:, j + 1);  % +1 for temperature offset
        end
    end

    w = [w; w_frac];

    % +++ End +++ %

    % Retrieve computed diameter values from base workspace
    d1 = evalin('base', 'd_values');
    t_d1 = evalin('base', 'd_times');

    d = [d; d1];
    t_d = [t_d; t_d1];

    % Check which event stopped the integration
    if isempty(ie)  % No event occurred
        break;
    elseif ie == 2  % Droplet shrank below 5%
        break;
    elseif ie == 1  % A component reached the threshold molar fraction
        mdi = we(2:end);  % Extract component masses
        % Extract droplet composition
        [~, ~, Xli, ~] = getLiquidComposition(Droplet, we(2:end));
        depleted_idx = find(Xli <= min_fraction, 1);  % Find first depleted component

        if isscalar(mdi)  % Only one component left
            break;
        end

        % Remove the depleted component
        activeComponents(depleted_idx) = [];  % Remove from active set
        currentConditions = we;
        currentConditions(depleted_idx + 1) = [];  % Remove corresponding mdi

        % Restart integration from the last event time
        tspan = [te, tend];
    end
end

%% --- POST-PROCESSING --- %%

% Diameter evolution
figure;
til = tiledlayout(1,1);
til.Padding = 'none';
til.TileSpacing = 'none';
nexttile
plot(t_d.*1e+03, (d./d0).^2, 'Color', 'red', 'LineWidth', 2);
box on;
xlabelName = '$t \ [ms]$';
xlabelHandle = xlabel(xlabelName);
set(xlabelHandle, 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex');
ylabelName = '$(d/d_0)^2 \ [-]$';
ylabelHandle = ylabel(ylabelName);
set(ylabelHandle, 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex');
set(0,'defaultfigurecolor',[1 1 1]);
axesHandle = gca;
axesHandle.Box = 'on';          % draws the box
axesHandle.Layer = 'top';       % ensures axes lines are on top
set(axesHandle, 'FontSize', 20, 'FontName', 'Times', 'Linewidth', 1.2);
set(gcf, 'Position', [100, 100, 550, 550/1.3991]);

% Temperature evolution
figure;
til = tiledlayout(1,1);
til.Padding = 'none';
til.TileSpacing = 'none';
nexttile
plot(t.*1e+03, w(:,1), 'Color', 'red', 'LineWidth', 2);
box on;
xlabelName = '$t \ [ms]$';
xlabelHandle = xlabel(xlabelName);
set(xlabelHandle, 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex');
ylabelName = '$T \ [K]$';
ylabelHandle = ylabel(ylabelName);
set(ylabelHandle, 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex');
set(0,'defaultfigurecolor',[1 1 1]);
axesHandle = gca;
axesHandle.Box = 'on';          % draws the box
axesHandle.Layer = 'top';       % ensures axes lines are on top
set(axesHandle, 'FontSize', 20, 'FontName', 'Times', 'Linewidth', 1.2);
set(gcf, 'Position', [100, 100, 550, 550/1.3991]);

