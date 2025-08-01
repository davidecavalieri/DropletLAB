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

function [mu,mui] = calcMuCHEMKIN(Xi,T,Species)

%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routine calculates the individual component viscosities using the
% standard kinetic theory formulation. The mixture viscosity is then calculated
% using the formula of Wilke et al.

% Inputs:
% 1) Xi  : Mole Fractions
% 2) T   : Temperature
% 3) Species:  Struct variable containing the Species inputs parameters

% Outputs:
% 1) mu      : Mixture viscosity [Pa*s]
% 2) mui     : Single component viscosity [Pa*s]

% REFs:
% (i) Kee, R J, Dixon-Lewis, G, Warnatz, J, Coltrin, M E, and Miller, J A.
% A Fortran computer code package for the evaluation of gas-phase, multicomponent
% transport properties. United States: N. p., 1986.

% (ii) Bird, R., Stewart, W. and Lightfoot, E. (2002) Transport Phenomena.
% 2nd Edition, John Wiley and Sons, New York.

% LENNARD-JONES POTENTIAL WELL DEPTH [K]
eps_kb = cell2mat({Species.eps_kb});
% REDUCED TEMPERATURE [-]
KT_eps = T./eps_kb;
% LENNARD-JONES COLLISION DIAMETER [m]
sigma  = (cell2mat({Species.sigma}))*1e-10;
% MOLECULAR WEIGHT [Kg mol-1]
Mi     = cell2mat({Species.M});
% BOLTZAMANN CONSTANT [m2 kg s-2 K-1]
kb = 1.3806488e-23;
% AVOGADRO'S NUMBER [mol-1]
NA = 6.02214076e23;
% DIPOLE MOMENT
dipM  = (cell2mat({Species.mu}))*1e-21/299792458;

delta = 1/2*(dipM./sqrt(eps_kb.*kb.*sigma.^3)).^2;
Omega = collInt(KT_eps,delta);

mui = 5/16*sqrt(kb*pi*Mi.*T/NA)./(pi*(sigma.^2).*(Omega));

%%%%%%%%%%%%%%%%%%%%%%%%%%% MIXTURE VISCOSITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mu] = muMixWilke(mui,Mi*1e+03,Xi);

%% COLLISION INTEGRAL DATA, Bird et. al. (2002)
    function out = collInt(Tred,delta)
        out = zeros(1,length(delta));
        for k = 1:length(delta)

            if delta(k) == 0
                exponents = [1 2 3 4 5 6];
                cc2A = -0.92032979;
                cc2B = [2.3508044 0.50110649 -4.7193769/10 1.5806367/10 -2.6367184/100 1.8120118/1e3];
                cc2C = [1.6330213 -6.9795156/10 1.6096572/10 -2.2109440/100 1.7031434/1e3 -0.56699986/1e4];
                out(k) = cc2A + sum(cc2B./Tred(k).^exponents + cc2C.*log(Tred(k)).^exponents);
            else
                deltaV = [0.0, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 2.5];

                TredV = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0...
                    1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0...
                    5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0...
                    18.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 75.0, 100.0];

                coll22 = [4.1005, 4.266,  4.833,  5.742,  6.729,  8.624,  10.34,  11.89;
                    3.2626, 3.305,  3.516,  3.914,  4.433,  5.57,   6.637,  7.618;
                    2.8399, 2.836,  2.936,  3.168,  3.511,  4.329,  5.126,  5.874;
                    2.531,  2.522,  2.586,  2.749,  3.004,  3.64,   4.282,  4.895;
                    2.2837, 2.277,  2.329,  2.46,   2.665,  3.187,  3.727,  4.249;
                    2.0838, 2.081,  2.13,   2.243,  2.417,  2.862,  3.329,  3.786;
                    1.922,  1.924,  1.97,   2.072,  2.225,  2.614,  3.028,  3.435;
                    1.7902, 1.795,  1.84,   1.934,  2.07,   2.417,  2.788,  3.156;
                    1.6823, 1.689,  1.733,  1.82,   1.944,  2.258,  2.596,  2.933;
                    1.5929, 1.601,  1.644,  1.725,  1.838,  2.124,  2.435,  2.746;
                    1.4551, 1.465,  1.504,  1.574,  1.67,   1.913,  2.181,  2.451;
                    1.3551, 1.365,  1.4,    1.461,  1.544,  1.754,  1.989,  2.228;
                    1.28,   1.289,  1.321,  1.374,  1.447,  1.63,   1.838,  2.053;
                    1.2219, 1.231,  1.259,  1.306,  1.37,   1.532,  1.718,  1.912;
                    1.1757, 1.184,  1.209,  1.251,  1.307,  1.451,  1.618,  1.795;
                    1.0933, 1.1,    1.119,  1.15,   1.193,  1.304,  1.435,  1.578;
                    1.0388, 1.044,  1.059,  1.083,  1.117,  1.204,  1.31,   1.428;
                    0.99963, 1.004, 1.016,  1.035,  1.062,  1.133,  1.22,   1.319;
                    0.96988, 0.9732, 0.983,  0.9991, 1.021,  1.079,  1.153,  1.236;
                    0.92676, 0.9291, 0.936,  0.9473, 0.9628, 1.005,  1.058,  1.121;
                    0.89616, 0.8979, 0.903,  0.9114, 0.923,  0.9545, 0.9955, 1.044;
                    0.87272, 0.8741, 0.878,  0.8845, 0.8935, 0.9181, 0.9505, 0.9893;
                    0.85379, 0.8549, 0.858,  0.8632, 0.8703, 0.8901, 0.9164, 0.9482;
                    0.83795, 0.8388, 0.8414, 0.8456, 0.8515, 0.8678, 0.8895, 0.916;
                    0.82435, 0.8251, 0.8273, 0.8308, 0.8356, 0.8493, 0.8676, 0.8901;
                    0.80184, 0.8024, 0.8039, 0.8065, 0.8101, 0.8201, 0.8337, 0.8504;
                    0.78363, 0.784,  0.7852, 0.7872, 0.7899, 0.7976, 0.8081, 0.8212;
                    0.76834, 0.7687, 0.7696, 0.7712, 0.7733, 0.7794, 0.7878, 0.7983;
                    0.75518, 0.7554, 0.7562, 0.7575, 0.7592, 0.7642, 0.7711, 0.7797;
                    0.74364, 0.7438, 0.7445, 0.7455, 0.747,  0.7512, 0.7569, 0.7642;
                    0.71982, 0.72,   0.7204, 0.7211, 0.7221, 0.725,  0.7289, 0.7339;
                    0.70097, 0.7011, 0.7014, 0.7019, 0.7026, 0.7047, 0.7076, 0.7112;
                    0.68545, 0.6855, 0.6858, 0.6861, 0.6867, 0.6883, 0.6905, 0.6932;
                    0.67232, 0.6724, 0.6726, 0.6728, 0.6733, 0.6743, 0.6762, 0.6784;
                    0.65099, 0.651,  0.6512, 0.6513, 0.6516, 0.6524, 0.6534, 0.6546;
                    0.61397, 0.6141, 0.6143, 0.6145, 0.6147, 0.6148, 0.6148, 0.6147;
                    0.5887, 0.5889, 0.5894, 0.59,   0.5903, 0.5901, 0.5895, 0.5885];

                [X,Y] = meshgrid(TredV,deltaV);
                out(k) = interp2(X,Y,coll22',Tred(k),delta(k),"spline");
            end
        end
    end
end