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

function lnphi_i = fug_coef(Species,X,bi, bm, a_jk, am, Z, T, P, EoS )

R0 = 8.314472;

%% NORMALIZED PARAMETER
A = ( am * P ) / ( R0 * T )^2;
B = ( bm * P ) / ( R0 * T );

if strcmp(EoS,'PR')

    delta1 = 1 + sqrt(2); delta2 = 1 - sqrt(2);

    lnphi_i =  (bi / bm) * (Z - 1) - log(Z - B) ...
        - A / ( B * (delta1 - delta2) ) ...
        * ( 2 * a_jk / am - bi / bm  ) ...
        * log((Z + B * delta1)/(Z + B * delta2));

elseif strcmp(EoS,'SRK')

    delta1 = 1.0;         delta2 = 0;

    lnphi_i =  (bi / bm) * (Z - 1) - log(Z - B) ...
        - A / ( B * (delta1 - delta2) ) ...
        * ( 2 * a_jk / am - bi / bm  ) ...
        * log((Z + B * delta1)/(Z + B * delta2));

elseif strcmp(EoS,'RKPR')

    d1 = 0.428363; d2 = 18.496215;  d3 = 0.338426;
    d4 = 0.660000; d5 = 789.723105; d6 = 2.512392;

    for k = 1:numel(Species)
        delta1_i(k) = d1 + d2*((d3 - 1.168*Species(k).Zc)^d4) + d5*((d3 - 1.168*Species(k).Zc)^d6);
        delta2_i(k) = (1-delta1_i(k))/(1+delta1_i(k));
    end
    delta1 = sum(X.*delta1_i,2);
    delta2 = sum(X.*delta2_i,2);

    u  = delta1 + delta2;
    uw = sqrt( (delta1 + delta2)^2 - 4*delta1*delta2 );

    lnphi_i =  (bi / bm) * (Z - 1) - log(Z - B) ...
        - A / ( B * uw) * ( 2 * a_jk / am - bi / bm  ) ...
        *log((2*Z + B*(u + uw))/(2*Z + B*(u - uw)));
else
    error('The selected EoS is not available')
end

end
