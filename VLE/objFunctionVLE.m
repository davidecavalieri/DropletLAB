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

function [f] = objFunctionVLE(W,XliD,Species,T,Pamb,Ns,EoS,XviA_inf)

%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routine construct the objective function

% HERE WE IMPOSE THAT THE LIQUID COMPOSITION OF THE AMBIENT AT THE DROPLET
% SURFACE IS NULL, I.E. NO SOLUBILITY EFFECT IS ACCOUNTED!
XliA = zeros(1,numel(XviA_inf));

% LIQUID PHASE COMPOSITION AT THE DROPLET SURFACE
Xli = [XliD XliA];

% VAPOR PHASE COMPOSITION OF THE DROPLET SPECIES.
Xvsi = W;

% VAPOR PHASE COMPOSITION OF THE AMBIENT 
XvA = 1 - sum(Xvsi);

% VAPOR PHASE COMPOSITION OF THE AMBIENT SPECIES (e.g AIR)
XviA = XvA.*XviA_inf;

% VAPOR PHASE COMPOSITION OF DROPLET+AMBIENT SPECIES
Xvi = [Xvsi XviA];

% LIQUID PHASE NON IDEAL PARAMETERS and COMPRESSIBILITY FACTOR
[am_L, bm_L, ~,~,sum_ref_i_L,bc] = VdWmixing(Species,Xli,T,EoS);
[Z_L] = getZ(Species,am_L,bm_L,Pamb,T,Xli,EoS);

% VAPOUR PHASE NON IDEAL PARAMETERS and COMPRESSIBILITY FACTOR
[am_V, bm_V, ~,~,sum_ref_i_V,~] = VdWmixing(Species,Xvi,T,EoS);
[Z_V] = getZ(Species,am_V,bm_V,Pamb,T,Xvi,EoS);

lnphi_L = zeros(1, Ns);
lnphi_V = zeros(1, Ns);
K = zeros(1, Ns);
f = zeros(1, Ns);

for i = 1:Ns

   lnphi_L(i) = fug_coef(Species,Xli,bc(i), bm_L, sum_ref_i_L(i), am_L, Z_L, T, Pamb,EoS );
   lnphi_V(i) = fug_coef(Species,Xvi,bc(i), bm_V, sum_ref_i_V(i), am_V, Z_V, T, Pamb,EoS );
   K(i) = Xvi(i) / Xli(i);

end

f(1:Ns) = lnphi_L - lnphi_V - log(K);

end