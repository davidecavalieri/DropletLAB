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

function [am, bm, dadTm,ddadTm,sum_ref_i,bc] = VdWmixing(Species,X,T,EoS)

%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VdWmixing calculates the non-ideal coefficients and their derivatives of a cubic
% equation of state, using the single-fluid-mixing model and pseudo-critical
% combination rules

% Inputs:
% 1) Species     : Struct variable containing the Species inputs parameters
% 2) X  : Mole Fractions
% 3) T  : Temperature [K]

% Outputs:
% 1) am     : Mixture attractive force term, am(X,T)
% 2) bm     : Mixture Covolume, bm(X,T)
% 3) dadTm  : Mixture first order partial derivatives of am w.r.t T, dadTm(X,T)
% 4) ddadTm : Mixture second order partial derivatives of am w.r.t T, ddadTm(X,T)

R0 = 8.314472;
%%%%%%%%%%%%%%%% TODO Matrix of binary interaction parameters  %%%%%%%%%%%
kij = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% A,B,C EoS Parameters %%%%%%%%%%%%%%%%%%%%%%%%%
Ns = numel(X);
ac = zeros(1,Ns); bc = zeros(1,Ns); cc =  zeros(1,Ns);

for k = 1:Ns
    if strcmp(EoS,'PR')
        [ac(k),bc(k),cc(k)] = PR_abc(Species(k).Tc,Species(k).Pc,Species(k).omega);
    elseif strcmp(EoS,'SRK')
        [ac(k),bc(k),cc(k)] = SRK_abc(Species(k).Tc,Species(k).Pc,Species(k).omega);
    elseif strcmp(EoS,'RKPR')
        [ac(k),bc(k),cc(k)] = RKPR_abc(Species(k).Tc,Species(k).Pc,Species(k).omega,Species(k).Zc);
    else
        error('The selected EoS is not available')
    end
end

%%  van der Waals one-fluid mixing rules %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
am          = 0.;
bm          = 0.;
dadTm       = 0.;
ddadTm      = 0.;
sum_ref_i   = zeros(1,Ns);
aij   = zeros(Ns,Ns);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:Ns
    ref = 0.;
    for j=1:Ns
        if i == j
            a_ij = ac(i);
            [alpha_ij, dalphadT_ij, ddalphadT_ij] = alphaFunction(T, Species(i).Tc, cc(i), EoS );
        else
            %%%%%%%%%%%%%%%%%%%%% ELEMENTI OFF-DIAGONAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            vc_ij    = 0.125 * ( Species(i).Vc^(1/3) + Species(j).Vc^(1/3))^3;
            Zc_ij    = 0.5   * ( Species(i).Zc + Species(j).Zc);
            Tc_ij    = sqrt( Species(i).Tc * Species(j).Tc)*(1 - kij);
            pc_ij    = Zc_ij * R0 * Tc_ij / vc_ij;
            omega_ij = 0.5 * ( Species(i).omega + Species(j).omega);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EoS choice %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if strcmp(EoS,'PR')
                [a_ij,~,c_ij] = PR_abc(Tc_ij,pc_ij,omega_ij);
            elseif strcmp(EoS,'SRK')
                [a_ij,~,c_ij] = SRK_abc(Tc_ij,pc_ij,omega_ij);
            elseif strcmp(EoS,'RKPR')
                [a_ij,~,c_ij] = RKPR_abc(Tc_ij,pc_ij,omega_ij,Zc_ij);
            else
                error('The selected EoS is not available')
            end
            [alpha_ij, dalphadT_ij, ddalphadT_ij] = alphaFunction( T, Tc_ij, c_ij, EoS );
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        am     = am     + X(i) * X(j) * a_ij * alpha_ij;
        dadTm  = dadTm  + X(i) * X(j) * a_ij * dalphadT_ij;
        ddadTm = ddadTm + X(i) * X(j) * a_ij * ddalphadT_ij;
        ref    = ref    +          X(j) * a_ij * alpha_ij;
        aij(i,j) = a_ij * alpha_ij;
    end
    bm             = bm + X(i) * bc(i);
    sum_ref_i(i)  = ref;

end

end