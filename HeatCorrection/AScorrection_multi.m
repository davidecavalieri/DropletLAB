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

function [Sh_AS_i,Sh_AS,Nu_AS] = AScorrection_multi(Bm,cpvap_f_mix,cp_f,Pr_f,Sh_f0_bulk,Nu_f0,Sh_f0,Sc_f_bulk,Sc_f)

%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routine returns the modified Sherwood and Nusselt numbers using the
% Abramzon & Sirignano Model

% Inputs:
% 1) Bm   : Spalding mass number
% 2) cpvap_f_mix  : Specific heat of droplet species
% 3) cp_f   : Specific heat of droplet + ambient species
% 4) Pr_f   : Prandtl number
% 5) Sh_f0_bulk : Uncorrected averaged Sherwood number
% 6) Nu_f0 : Uncorrected Nusselt number
% 7) Sh_f0 : Uncorrected individual Sherwood numbers
% 8) Sc_f_bulk : Averaged Schmidt number
% 9) Sc_f : Individual Schmidt numbers

% Outputs:
% 1) Sh_AS : Corrected averaged Sherwood number
% 2) Nu_AS : Corrected Nusselt number

% CLIPPING OF THE SPALDING MASS TRANSFER NUMBER
if Bm >= 20
    Bm = 20;
end

% FILM THICKNESS CORRECTION FACTOR (MASS)
Fm = ((1 + Bm).^0.7).*(log(1 + Bm)./Bm);

% EFFECTIVE SHERWOOD NUMBER
Sh_AS = 2 + (Sh_f0_bulk - 2)./Fm;

% UNCORRECTED PHI
phi = (cpvap_f_mix./cp_f).*(Pr_f/Sc_f_bulk).*(Sh_f0_bulk/Nu_f0);

% OBJECTIVE FUNCTION FOR THE MODIFIED NUSSELT NUMBER
NUFUN = @(Nu_AS) objNu_star(Nu_AS,Sh_AS,Sh_f0_bulk,Nu_f0,phi,Bm);

options=optimset('display','off');
Nu_AS_guess = Nu_f0;
[Nu_AS,~,~] = fzero(NUFUN,Nu_AS_guess,options);

Sh_AS_i = zeros(1,numel(Sh_f0));

for i = 1:numel(Sh_f0)

    % GUESS SHERWOOD NUMBER FOR THE i-th DROPLET SPECIES
    Sh_AS_i_0 = 2 + (Sh_f0(i) - 2)./Fm;

    % OBJECTIVE FUNCTION FOR THE MODIFIED i-th SHERWOOD NUMBER
    SHFUN = @(Sh_AS_i_f) objSh_star(Sh_AS_i_f,Sh_AS,Sc_f(i),Sc_f_bulk,Sh_AS_i_0,Bm);

    options=optimset('display','off');
    Sh_AS_i_guess = Sh_AS_i_0;
    [Sh_AS_i_f,~,~] = fzero(SHFUN,Sh_AS_i_guess,options);

    Sh_AS_i(i) = Sh_AS_i_f;

end



    function Fun = objNu_star(Nuf_fAS,Sh_AS,Sh_f0_bulk,Nu_f0,phi,Spalding)

        %%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % This routine constructs the objective function to minimize for the
        % computation of the modified Nusselt number.

        % CORRECTED PHI
        phi_star = phi.*(Sh_AS./Sh_f0_bulk).*(Nu_f0./Nuf_fAS);

        % CORRECTED SPALDING THERMAL NUMBER
        Bt_star = ((1 + Spalding).^phi_star) - 1;

        % FILM THICKNESS CORRECTION FACTOR (THERMAL)
        Ft_star = ((1 + Bt_star).^0.7).*(log(1 + Bt_star)./Bt_star);

        % OBJECTIVE FUNCTION
        Fun = 2 + (Nu_f0 - 2)./(Ft_star) - Nuf_fAS;

    end

    function Fun = objSh_star(Sh_AS_i,Sh_AS,Sc_f_i,Sc_f_bulk,Sh_AS_i_0,Spalding)

        %%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % This routine constructs the objective function to minimize for the
        % computation of the modified Sherwood number for the i-th droplet species.

        % NON-DIMENSIONAL PARTIAL EVAPORATION RATE FOR THE i-th DROPLET SPECIES
        eta_i = (Sh_AS*Sc_f_i)/(Sh_AS_i*Sc_f_bulk);

        % SPALDING MASS TRANSFER NUMBER FOR THE i-th DROPLET SPECIES
        Bm_i = (1 + Spalding).^eta_i - 1;

        % FILM THICKNESS CORRECTION FACTOR (MASS) FOR THE i-th DROPLET SPECIES
        Fm_i = ((1 + Bm_i).^0.7).*(log(1 + Bm_i)./Bm_i);

        % OBJECTIVE FUNCTION
        Fun = 2 + (Sh_AS_i_0 - 2)./(Fm_i) - Sh_AS_i;

    end

end