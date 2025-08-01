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

function [mu,lambda] = transport_mu_k_CHEMKIN(T,Xi,cpi,Species)

%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routine calculates the individual component viscosities and
% conductivities using the standard kinetic theory formulation.

% Inputs:
% 1) T     : Temperature
% 2) cpi   : Species isobaric specific heat  [J/kgK]
% 3) Species:  Struct variable containing the Species inputs parameters

% Output:
% 1) lambdai  : Single species thermal conductivity [w/mK]

% UNIVERSAL GAS CONSTANT
R0       = 8.314472;
% BOLTZAMANN CONSTANT [m2 kg s-2 K-1]
kb       = 1.3806488e-23;
% AVOGADRO'S NUMBER [mol-1]
NA       = 6.02214076e23;
% SPECIFIC HEAT AT CONSTANT VOLUME TRANSITIONAL [-]
cv_trans = 3/2*R0;


% MOLECULAR WEIGHT [Kg mol-1]
Mi      = (cell2mat({Species.M}));
% GEOMETRY TYPE (0,1 or 2)
GeoType = cell2mat({Species.geometry});
% LENNARD-JONES POTENTIAL WELL DEPTH [K]
eps_kb  = cell2mat({Species.eps_kb});
% LENNARD-JONES COLLISION DIAMETER [m]
sigma   = (cell2mat({Species.sigma}))*1e-10;
% DIPOLE MOMENT
dipM    = (cell2mat({Species.mu}))*1e-21/299792458;
% ROTATIONAL RELAXATION COLLISION NUMBER AT 298 K
Zrot0   = cell2mat({Species.Zrot});
% REDUCED TEMPERATURE [-]
Tred    = T./eps_kb;
% DELTA
delta = 1/2*(dipM./sqrt(eps_kb.*kb.*sigma.^3)).^2;

cv_rot = zeros(1,length(Mi));
% SPECIFIC HEAT AT CONSTANT VOLUME ROTATIONAL [-]
for i = 1:numel(Mi)
    if GeoType(i)     == 0
        cv_rot(i) = 0;
    elseif GeoType(i) == 1
        cv_rot(i) = 1;
    elseif GeoType(i) == 2
        cv_rot(i) = 3/2*R0;
    else
        error('error in type')
    end
end

% SPECIFIC HEAT AT CONSTANT VOLUME [J/molK]
cvi     = cpi.*Mi - R0;
% VIBRATIONAL SPECIFIC HEAT AT CONSTANT VOLUME [J/molK]
cv_vib  = cvi - cv_trans - cv_rot;
% COLLISION INTEGRAL 1
OMEGA_1 = collInt(1,real(Tred),delta);
% COLLISION INTEGRAL 2
OMEGA_2 = collInt(2,real(Tred),delta);
% SELF DIFFUSION COEFFICIENTS [m2/s]
pDii    = 3/16*sqrt(2*pi*kb^3*(T.^3)./(Mi./2./NA))./(pi.*sigma.^2.*OMEGA_1);
% SPECIES VISCOSITIES [Pa s]
mui     = 5/16*sqrt(pi.*Mi.*kb.*T./NA)./(pi*sigma.^2.*OMEGA_2);
% VIBRATIONAL CONTRIBUTIONS [-]
f_vib   = Mi./R0./T.*pDii./mui;
% PARKER BRAU & JONKMAN
Tred298 = 298./eps_kb;
F298    = 1+(pi^1.5)./sqrt(Tred298).*(0.5 + 1./Tred298) + (0.25*pi*pi+2)./Tred298;
F       = 1+(pi^1.5)./sqrt(Tred).*(0.5 + 1./Tred) + (0.25*pi*pi+2)./Tred;
Zrot    = Zrot0.*F298./F;
% ROTATIONAL CONTRIBUTIONS  [-]
A       = 5/2 - f_vib;
B       = Zrot + 2/pi*(5/3*cv_rot/R0 + f_vib);
C       = (2/pi)*(A./B);
f_rot   = f_vib.*(1+C);
% TRANSLATIONAL CONTRIBUTIONS [-]
f_trans = 5/2*(1 - 2/pi*cv_rot./cv_trans.*A./B);
% SINGLE SPECIES CONDUCTIVITIES [W/mK]
lambdai = mui./Mi.*(f_trans.*cv_trans + f_rot.*cv_rot + f_vib.*cv_vib);

%%%%%%%%%%%%%%%%%%%%%%%%%%% MIXTURE VISCOSITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu = muMixWilke(mui,Mi*1e+03,Xi);

% MIXTURE-AVERAGED THERMAL CONDUCTIVITY [W/mK]
term1 = sum(Xi.*lambdai);
term2 = 1./(sum(Xi./lambdai));
lambda = 0.5*(term1 + term2);

end


function OMEGA = collInt(n,Tred,delta)
OMEGA = zeros(1,length(delta));
for k = 1:length(delta)
    if n==1
        if delta(k) == 0
            cc1A = -1.1036729;
            cc1B = [2.6431984 0.0060432255 -1.5158773/10 0.54237938/10 -0.90468682/100 0.61742007/1e3];
            cc1C = [1.6690746 -6.9145890/10 1.5502132/10 -2.0642189/100 1.5402077/1e3 -0.49729535/1e4];
            exponents = [1 2 3 4 5 6];
            OMEGA(k) = cc1A + sum(cc1B./Tred(k).^exponents + cc1C.*log(Tred(k)).^exponents);
        else
            deltaV = [0.0, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 2.5];

            TredV = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0...
                1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0...
                5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0...
                18.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 75.0, 100.0];

            coll11 = [4.00791711465155	4.00187617260788	4.65606936416185	5.52115384615385	6.45158197507191	8.21333333333333	9.82889733840304	11.3130352045671
                3.12989255564083	3.16267942583732	3.35496183206107	3.72053231939164	4.19791666666667	5.23004694835681	6.22607879924953	7.15977443609023
                2.64940759399198	2.65791940018744	2.76981132075472	3.00284360189574	3.31852551984877	4.05337078651685	4.78618113912232	5.48459383753501
                2.31437454279444	2.32014719411224	2.40111420612813	2.57156220767072	2.81273408239700	3.38604651162791	3.97217068645640	4.54081632653061
                2.06613589070841	2.07377049180328	2.14062500000000	2.27777777777778	2.47217068645640	2.94547134935305	3.43819188191882	3.91974169741697
                1.87662103746398	1.88496376811594	1.94343065693431	2.05968778696051	2.22559852670350	2.62809917355372	3.05412844036697	3.47339449541284
                1.72935036890409	1.73803071364047	1.79090909090909	1.89223744292237	2.03568161024703	2.38721461187215	2.76277372262774	3.13698630136986
                1.61221181556196	1.62149954832882	1.66969147005445	1.75978161965423	1.88524590163934	2.19727272727273	2.53454545454545	2.87170154686078
                1.51749954898070	1.52712477396022	1.57259528130672	1.65304268846503	1.76566757493188	2.04343891402715	2.34932126696833	2.65670289855072
                1.43984452680105	1.45018115942029	1.49048050770626	1.56391659111514	1.66485507246377	1.91696750902527	2.19567177637511	2.47833935018051
                1.32041742286751	1.32940108892922	1.36355394378966	1.42443438914027	1.50858175248419	1.72032374100719	1.95605381165919	2.19820627802691
                1.23359126081020	1.24203821656051	1.27157129881926	1.32336956521739	1.39350180505415	1.57309417040359	1.77747989276139	1.98928571428571
                1.16788321167883	1.17609489051095	1.20200181983621	1.24569356300997	1.30595667870036	1.46057347670251	1.63960749330955	1.82651245551601
                1.11660422187700	1.12420091324201	1.14558689717925	1.18511796733212	1.23646209386282	1.37153088630260	1.52983081032947	1.69804618117229
                1.07526980062191	1.08226691042048	1.10209662716500	1.13520871143376	1.18066847335140	1.30017921146953	1.44078361531612	1.59131205673759
                1.00064067362255	1.00548446069470	1.02005469462170	1.04640582347589	1.07963800904977	1.16950672645740	1.27782724844167	1.39646017699115
                0.950064020486556	0.953424657534247	0.965360072926162	0.985441310282075	1.01177536231884	1.08176100628931	1.16755793226381	1.26483613817538
                0.913070880526124	0.916058394160584	0.925318761384335	0.940909090909091	0.962828649138713	1.01888489208633	1.09025915996425	1.17036379769299
                0.884523483812129	0.887146763901550	0.894449499545041	0.907447774750227	0.924818840579710	0.972072072072072	1.03130590339893	1.09769094138544
                0.842738928798763	0.844636363636364	0.850136239782017	0.859618874773140	0.871312217194570	0.905405405405405	0.948028673835126	0.998219056099733
                0.812843537414966	0.814052583862194	0.817934782608696	0.824796380090498	0.834538878842676	0.859909909909910	0.892825112107623	0.931311329170384
                0.789791855203620	0.791040723981901	0.793851717902351	0.799006323396567	0.806407942238267	0.826372637263727	0.852466367713005	0.883303571428571
                0.771125361271676	0.772267389340560	0.774368231046931	0.779061371841155	0.784761045987376	0.800449640287770	0.821883408071749	0.847363717605005
                0.755522495717248	0.756357078449053	0.758701532912534	0.761801801801802	0.766426642664266	0.779694519317161	0.797757847533632	0.818588025022341
                0.742189610155758	0.742664266426643	0.744644464446445	0.747794779477948	0.751438848920863	0.762387791741472	0.777419354838710	0.795442359249330
                0.720237132848289	0.720287253141831	0.722282120395328	0.723967684021544	0.727199281867145	0.735515695067265	0.746374216651746	0.759964253798034
                0.702555137170522	0.703139013452915	0.703584229390681	0.705376344086022	0.707795698924731	0.714055505819158	0.722808586762075	0.733214285714286
                0.687737200143215	0.688182632050134	0.688988361683080	0.689803220035778	0.691681574239714	0.697137745974955	0.704021447721180	0.712767857142857
                0.675111746826390	0.675067024128686	0.675781948168007	0.676943699731904	0.678462913315460	0.682931188561215	0.688482142857143	0.695539696699376
                0.664023573533351	0.664107142857143	0.664732142857143	0.665625000000000	0.666964285714286	0.670115967885816	0.675200713648528	0.681105169340463
                0.641379310344828	0.641711229946524	0.642067736185383	0.642691622103387	0.643582887700535	0.645592163846839	0.649065004452360	0.652935943060498
                0.623472382815974	0.623754448398576	0.624021352313167	0.624466192170818	0.625088967971530	0.626957295373665	0.628977777777778	0.632177777777778
                0.608801847410960	0.608792184724689	0.609058614564831	0.609325044404974	0.609857904085258	0.611278863232682	0.613232682060391	0.615630550621670
                0.596398474230462	0.596628216503993	0.596805678793256	0.596983141082520	0.597426796805679	0.598314108251997	0.600000000000000	0.601418439716312
                0.576250331946535	0.576106194690266	0.576283185840708	0.576371681415929	0.576637168141593	0.577345132743363	0.578230088495575	0.579805137289637
                0.541467501543346	0.541534391534392	0.541710758377425	0.541409691629956	0.541585903083701	0.542151675485009	0.542151675485009	0.543021201413428
                0.518039422738472	0.517941952506596	0.518381706244503	0.518453427065026	0.518261633011414	0.518541300527241	0.518469656992085	0.518502202643172];

            [X,Y] = meshgrid(TredV,deltaV);
            OMEGA(k) = interp2(X,Y,coll11',Tred(k),delta(k),"spline");
        end
    elseif n==2
        if delta(k) == 0
            exponents = [1 2 3 4 5 6];
            cc2A = -0.92032979;
            cc2B = [2.3508044 0.50110649 -4.7193769/10 1.5806367/10 -2.6367184/100 1.8120118/1e3];
            cc2C = [1.6330213 -6.9795156/10 1.6096572/10 -2.2109440/100 1.7031434/1e3 -0.56699986/1e4];
            OMEGA(k) = cc2A + sum(cc2B./Tred(k).^exponents + cc2C.*log(Tred(k)).^exponents);
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
            OMEGA(k) = interp2(X,Y,coll22',Tred(k),delta(k),"spline");
        end
    end
end
end
