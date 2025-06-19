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

function [xi,yi,Xvsi,Xgsi,Yvsi,Ygsi] = VLE_EoS_solubility(Droplet,Ambient,Pamb,Xli,T,Psati,XviA_inf,YviA_inf,EoS) 

%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routine returns the droplet surface composition using the VLE
% concept with a cubic EoS

% Inputs:
% 1) Droplet : Struct variable containing the droplet inputs parameters
% 2) Ambient : Struct variable containing the droplet inputs parameters
% 3) Pamb    : Ambient pressure                             [Pa] 
% 3) Xli     : Liquid mole fractions of the droplet species [-]
% 4) Psati   : Saturated pressure of the droplet species    [Pa]

% Outputs:
% 1) Xvsi : Surface vapour mole fractions (droplet species)
% 2) Xgsi : Surface vapour mole fractions (Ambient species)
% 3) Yvsi : Surface vapour mass fractions (droplet species)
% 4) Ygsi : Surface vapour mass fractions (Ambient species)

Species = [Droplet Ambient];

% NUMBER OF DROPLET SPECIES
NsD = numel(Xli);

% NUMBER OF AMBIENT SPECIES
NsA = numel(XviA_inf);

% PERTURB THE INITIAL COMPOSITION OF THE DROPLET COMPOSITION
perturb = 1e-5;
XliD = Xli - perturb;

% INITAL LIQUID PHASE COMPOSITION OF AMBIENT SPECIES
XliA = (1 - sum(XliD)).*XviA_inf + perturb;

% INITAL LIQUID PHASE COMPOSITION OF DROPLET+AMBIENT SPECIES
Xli = [XliD XliA];

% INITIAL GUESS OF THE VAPOUR PHASE COMPOSION USING RAOULT'S LAW
[XvsiD_guess,~,~,~] = RaoultLaw(Droplet,Ambient,Pamb,XliD,Psati,XviA_inf,YviA_inf) ;
XvsiA = (1 - XvsiD_guess(1)).*XviA_inf;
Xvsi_guess = [XvsiD_guess XvsiA];

% STATE VECTOR FOR THE VLE PROBLEM: LIQUID VAPOR PHASE COMPOSITION.
W = [Xli Xvsi_guess];

options = optimset('TolFun',1e-06,'display','on','Algorithm','Levenberg-Marquardt'); 

VLE     = @(W) objFunctionVLE_solubility(W,Species,T,Pamb,NsD+NsA,EoS);
[W,~,~] = fsolve(VLE,W,options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LIQUID PHASE COMPOSITION (TOTAL) (MOLE FRACTION)
xi  = W(1:NsD+NsA);

% LIQUID PHASE COMPOSITION OF DROPLET SPECIES (MOLE FRACTION)
xiD = xi(1:NsD);

% LIQUID PHASE COMPOSITION OF AMBIENT SPECIES (MOLE FRACTION)
xiA = xi(NsD+1:end);

% VAPOR PHASE COMPOSITION (TOTAL) (MOLE FRACTION)
yi = W(NsD+NsA+1:end);

% VAPOR PHASE COMPOSITION OF DROPLET SPECIES (MOLE FRACTION)
Xvsi = yi(1:NsD);

% VAPOR PHASE COMPOSITION OF AMBIENT SPECIES (MOLE FRACTION)
Xgsi = yi(NsD+1:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MOLAR WEIGHT (kg/mol)
MiD = cell2mat({Droplet.M});
MiA = cell2mat({Ambient.M});

Mi = [MiD MiA];

% AFTER THE THERMODYNAMIC EQUILIBRIUM COMPUTATION THE LIQUID AND VAPOR PHASE
% COMPOSITION AT DROPLET INTERFACE IS KNOWN, WE CAN'T CALCULATE THE MASS FRACTION IN
% THE SAME WAY AS DONE WITH THE RAOULT LAW, IN THIS WAY WE ARE IMPOSING A
% WRONG RE-DISTRIBUTION OF SPECIES (E.G AIR)

[Yi,~] = Mole2Mass(Mi,yi);

% VAPOR PHASE COMPOSITION OF DROPLET SPECIES (MASS FRACTION)
Yvsi = Yi(1:NsD);

% VAPOR PHASE COMPOSITION OF AMBIENT SPECIES (MASS FRACTION)
Ygsi = Yi(NsD+1:end);


end