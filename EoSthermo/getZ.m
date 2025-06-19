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

function [Z] = getZ(Species,am,bm,P,T,X,EoS)

%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getZ calculate the mixture compressibility factor for a cubic EoS
% using Cardano's Methos and the Mixing Gibbs Criterion.

% Inputs:
% 1) Species     : Struct variable containing the Species inputs parameters
% 2) am     : Mixture attractive force term, am(X,T)
% 3) bm     : Mixture Covolume, bm(X,T)
% 4) T      : Temperature [K]
% 5) P      : Pressure   [Pa]

% Outputs:
% 1) Z      : Mixture compressibility factor

R0 = 8.314472;

% NUMBER OF SPECIES (DROPLET+AMBIENT)
Ns = numel(X);

%%%%%%%%%%%%%%%%%%% GENERALIZED CUBIC EQUATION OF STATE %%%%%%%%%%%%%%%%%%%
if strcmp(EoS,'PR')
    delta1 = 1 + sqrt(2); delta2 = 1 - sqrt(2);
elseif strcmp(EoS,'SRK')
    delta1 = 1.0;         delta2 = 0;
elseif strcmp(EoS,'RKPR')

    d1 = 0.428363; d2 = 18.496215;  d3 = 0.338426;
    d4 = 0.660000; d5 = 789.723105; d6 = 2.512392;

    for k = 1:Ns
        if Species(k).Zc >= 0.289
            Species(k).Zc  = 0.289; % Only for the helium
        end
        delta1_i(k) = d1 + d2*((d3 - 1.168*Species(k).Zc)^d4) + d5*((d3 - 1.168*Species(k).Zc)^d6);
        delta2_i(k) = (1-delta1_i(k))/(1+delta1_i(k));
    end
    delta1 = sum(X.*delta1_i,2);
    delta2 = sum(X.*delta2_i,2);

else
    error('The selected EoS is not available')
end

%%%%%%%%%%%%%%%% NON IDEAL PARAMETERS NORMALIZATION %%%%%%%%%%%%%%%%%%%%%%%
ad = 1./((R0^2)*(T.^2)); bd = 1./(R0*T);
A=(am.*P)*(ad); B=(bm*P)*(bd);
%%%%%%%%%%%%%%%%%%%%% COEFFICIENTS Z^3 +C2*Z^2 + C1*Z +C0 = 0 %%%%%%%%%%%%%
C2 = B*( delta1 + delta2 - 1) - 1;
C1 = A  + (delta1*delta2*B*B) - (delta1 + delta2)*B*(B+1);
C0 = -B*(delta1*delta2*B*B  + delta1*delta2*B  + A);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q = ( C2^2 - 3*C1 ) / 9.;
R = ( 2*C2^3 - 9*C2*C1 +27*C0 ) / 54.;

%%  THREE SOLUTION ZONE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if R^2 < Q^3
    sqrtQ = sqrt(Q);
    phi = acos( R / ( sqrtQ * Q ) );
    z1 = -2. * sqrtQ * cos(phi/3)        - C2/3;
    z2 = -2. * sqrtQ * cos((phi+2*pi)/3) - C2/3;
    z3 = -2. * sqrtQ * cos((phi-2*pi)/3) - C2/3;
    Z3roots = [z1, z2, z3];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if min(Z3roots) < B
        Z = max(Z3roots);
    else
        %%%%%%%%%%%%% Liquid to Vapor---> small to large Root %%%%%%%%%%%%%%%%%%%%%
        Zsort   = sort(Z3roots);
        Z2roots = [Zsort(1), Zsort(3)];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GIBBS CRITERION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ZL = Z2roots(1); % Liquid root
        ZV = Z2roots(2); % Vapor root
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GIBBS ENERGY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        delta_G = (A/(B*(delta1 -delta2))) ...
            * log(((ZL + delta1*B)*(ZV + delta2*B))/((ZL + delta2*B)*(ZV + delta1*B))) ...
            - ( ZL - ZV ) ...
            + log( ( ZL - B ) / ( ZV - B ) );
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if delta_G < 0
            Z = Z2roots(2); % Vapor root is the stable solution
        else
            Z = Z2roots(1); % Liquid root is the stable solution
        end
    end
    %% ONE NON IMMAGINARY SOLUTION ZONE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
else
    E = - R / abs(R) * ( abs(R) + sqrt( R^2 -Q^3 ) )^(1./3.);

    if E == 0
        F = 0.;
    else
        F = Q/E;
    end

    Z = ( E + F - C2/3 );
end

end




