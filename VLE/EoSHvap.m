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

function [Lvi] = EoSHvap(Droplet,Ambient,Xli,Xvi,T,Pamb)

%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R0 = 8.314;

Species = [Droplet Ambient];

MiD = cell2mat({Droplet.M});

% LIQUID PHASE NON IDEAL PARAMETERS and COMPRESSIBILITY FACTOR
[am_L, bm_L, ~,~,sum_ref_i_L,bc] = VdWmixing(Species,Xli,T);
[Z_L] = getZ(Species,am_L,bm_L,Pamb,T,Xli);

% VAPOUR PHASE NON IDEAL PARAMETERS and COMPRESSIBILITY FACTOR
[am_V, bm_V, ~,~,sum_ref_i_V,~] = VdWmixing(Species,Xvi,T);
[Z_V] = getZ(Species,am_V,bm_V,Pamb,T,Xvi);

deltaT = 1;

for i = 1:1

   lnphi_L_DeltaT(i) = fug_coef(Species,Xli,bc(i), bm_L, sum_ref_i_L(i), am_L, Z_L, T+deltaT, Pamb );
   lnphi_V_DeltaT(i) = fug_coef(Species,Xvi,bc(i), bm_V, sum_ref_i_V(i), am_V, Z_V, T+deltaT, Pamb );
 
   lnphi_L_mDeltaT(i) = fug_coef(Species,Xli,bc(i), bm_L, sum_ref_i_L(i), am_L, Z_L, T-deltaT, Pamb );
   lnphi_V_mDeltaT(i) = fug_coef(Species,Xvi,bc(i), bm_V, sum_ref_i_V(i), am_V, Z_V, T-deltaT, Pamb );
end

dlnphidT_L = (lnphi_L_DeltaT - lnphi_L_mDeltaT)./(2*deltaT);

dlnphidT_V = (lnphi_V_DeltaT - lnphi_V_mDeltaT)./(2*deltaT);

Lvi = -(R0*T*T.*(dlnphidT_V - dlnphidT_L))./MiD;

end
