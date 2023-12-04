%% ===== SET UP ===========================================================
clear('all'); close('all'); addpath('lib');

% --- Back code up --------------------------------------------------------
fname = 'EFAST_ode2MEv1strain';
mkdir(fname);
backupcode('RUN_efast_ode2MEv1strain', fname);

% --- Turn off warnings ---------------------------------------------------
warning('off','all');

% --- Open parallel pool --------------------------------------------------
if isunix
    parpool('local', str2num(getenv('SLURM_CPUS_PER_TASK')))
end

% --- Initiate random number seed -----------------------------------------
rngsettings = rng('shuffle');

%% ===== NOMINAL PARAMETER SET ============================================
% --- Time span -----------------------------------------------------------
tmax = 7*24;

% --- Initial conditions --------------------------------------------------
G0 = 20; % g per L
aA0 = 0; % mg per L
aB0 = 0; % mg per L
r0A = 0.5;
r0B = 0.5;

% --- Biomass parameters --------------------------------------------------
gamma_gluc   = 0.05; % Biomass yield of S. cerevisiae grown on glucose
gamma_A      = 0.5; % Biomass yield of S. cerevisiae grown on A
gamma_B      = 0.5; % Biomass yield of S. cerevisiae grown on B
% --- Strain A parameters (Genotype A+/B-) --------------------------------
Vmax_gluc_yA = 7.2;  % Maximum glucose uptake rate
Km_gluc_yA   = 5;    % Michaelis constant for glucose uptake
eta_yA       = 1e-3; % Mortality rate constant of strain A
Vmax_yA_B    = 30;   % Maximum B uptake rate by strain A
Km_yA_B      = 50;   % Michaelis constant for B uptake
delta_yA_A   = 1;    % Number of B molecules produced per molecule of glucose is consumed
% --- Strain B parameters (Genotype A-/B+) --------------------------------
Vmax_gluc_yB = 7.2;  % Maximum glucose uptake rate
Km_gluc_yB   = 5;    % Michaelis constant for glucose uptake
eta_yB       = 1e-3; % Mortality rate constant of strain A
Vmax_yB_A    = 30;   % Maximum B uptake rate by strain A
Km_yB_A      = 50;   % Michaelis constant for B uptake
delta_yB_B   = 1;    % Number of B molecules produced per molecule of glucose is consumed

% --- Metabolic engineering parameters ------------------------------------
delta_yA_X   = 1;    % Number of X molecules produced per molecule of glucose is consumed
phi_yA_X     = 0.5;  % Train to new product
Vmax_yB_XY   = 30;   % Maximum X uptake rate by strain B
Km_yB_XY     = 50;   % Michaelis constant for X uptake

% ----- Leaky vectors -----------------------------------------------------
omega(1) = 0.1;      % phi_yA_A
omega(2) = 0.1;      % phi_yB_B

% --- Build vectors -------------------------------------------------------
theta_biomass = [gamma_gluc gamma_A gamma_B];
theta_yA      = [Vmax_gluc_yA Km_gluc_yA eta_yA Vmax_yA_B Km_yA_B delta_yA_A];
theta_yB      = [Vmax_gluc_yB Km_gluc_yB eta_yB Vmax_yB_A Km_yB_A delta_yB_B];
theta_meta    = [delta_yA_X phi_yA_X Vmax_yB_XY Km_yB_XY];

% --- Set up options ------------------------------------------------------
odeoptions = odeset('NonNegative', 1:9);

%% ===== E-FAST SETTINGS ==================================================
% --- e-FAST settings -----------------------------------------------------
initalpha = 0.01; % p-value threshold
NR = 100;         % number of resampling
NS = 5*257;       % samples per search
MI = 4;           % maximum number of fourier coeffs retained
delta = 0.5;      % lb/ub bounds

% --- Initial conditons ---------------------------------------------------
dN0 = 0.03;

% --- Build k vector ------------------------------------------------------
kmidpoint     = zeros(1,23);
kmidpoint( 1) = dN0; % total initial population
kmidpoint( 2) = r0A;  % A ratio
kmidpoint( 3) = r0B;  % B ratio
kmidpoint( 4) = G0;  % glucose conc.
kmidpoint( 5) = aA0; % amino acid A conc
kmidpoint( 6) = aB0; % amino acid B conc
kmidpoint( 7: 9) = theta_biomass;
kmidpoint(10:15) = theta_yA;
kmidpoint(16:21) = theta_yB;
kmidpoint(22:23) = omega;
kmidpoint(24:27) = theta_meta;

% --- Vector names --------------------------------------------------------
klistnames{ 1} = 'dN0';          kc0( 1) = 0; klb( 1) = 0.01; kub( 1) = 1;     klist( 1,:) = [0];
klistnames{ 2} = 'r0A';          kc0( 2) = 1; klb( 2) = 1e-2; kub( 2) = 1;     klist( 2,:) = [1];
klistnames{ 3} = 'r0B';          kc0( 3) = 1; klb( 3) = 1e-2; kub( 3) = 1;     klist( 3,:) = [1];
klistnames{ 4} = 'G0';           kc0( 4) = 0; klb( 4) = 0;    kub( 4) = 1;     klist( 4,:) = [0];
klistnames{ 5} = 'aA0';          kc0( 5) = 0; klb( 5) = 0;    kub( 5) = 75;    klist( 5,:) = [1];
klistnames{ 6} = 'aB0';          kc0( 6) = 0; klb( 6) = 0;    kub( 6) = 75;    klist( 6,:) = [1];
klistnames{ 7} = 'gamma_gluc';   kc0( 7) = 0; klb( 7) = 1e-2; kub( 7) = 1;     klist( 7,:) = [0];
klistnames{ 8} = 'gamma_A';      kc0( 8) = 0; klb( 8) = 1e-2; kub( 8) = 1;     klist( 8,:) = [0];
klistnames{ 9} = 'gamma_B';      kc0( 9) = 0; klb( 9) = 1e-2; kub( 9) = 1;     klist( 9,:) = [0];
% --- Strain A ------------------------------------------------------------
klistnames{10} = 'Vmax_gluc_yA'; kc0(10) = 0; klb(10) = 1;    kub(10) = 30;    klist(10,:) = [1];
klistnames{11} = 'Km_gluc_yA';   kc0(11) = 0; klb(11) = 1;    kub(11) = 100;   klist(11,:) = [0];
klistnames{12} = 'eta_yA';       kc0(12) = 0; klb(12) = 1e-4; kub(12) = 1e-2;  klist(12,:) = [0];
klistnames{13} = 'Vmax_yA_B';    kc0(13) = 0; klb(13) = 1;    kub(13) = 120;   klist(13,:) = [1];
klistnames{14} = 'Km_yA_B';      kc0(14) = 0; klb(14) = 1;    kub(14) = 1000;  klist(14,:) = [0];
klistnames{15} = 'delta_yA_A';   kc0(15) = 0; klb(15) = 1e-2; kub(15) = 1;     klist(15,:) = [0];
% --- Strain B ------------------------------------------------------------
klistnames{16} = 'Vmax_gluc_yB'; kc0(16) = 0; klb(16) = 1;    kub(16) = 30;    klist(16,:) = [1];
klistnames{17} = 'Km_gluc_yB';   kc0(17) = 0; klb(17) = 1;    kub(17) = 100;   klist(17,:) = [0];
klistnames{18} = 'eta_yB';       kc0(18) = 0; klb(18) = 1e-4; kub(18) = 1e-2;  klist(18,:) = [0];
klistnames{19} = 'Vmax_yB_A';    kc0(19) = 0; klb(19) = 1;    kub(19) = 120;   klist(19,:) = [1];
klistnames{20} = 'Km_yB_A';      kc0(20) = 0; klb(20) = 1;    kub(20) = 1000;  klist(20,:) = [0];
klistnames{21} = 'delta_yB_B';   kc0(21) = 0; klb(21) = 1e-2; kub(21) = 1;     klist(21,:) = [0];
% --- Leaks ---------------------------------------------------------------
klistnames{22} = 'phi_yA_A';     kc0(22) = 0; klb(22) = 1e-2; kub(22) = 0.5;   klist(22,:) = [1];
klistnames{23} = 'phi_yB_B';     kc0(23) = 0; klb(23) = 1e-2; kub(23) = 0.5;   klist(23,:) = [1];
% --- Metabolic engineering parameters ------------------------------------
klistnames{24} = 'delta_yA_X';   kc0(24) = 0; klb(24) = 1e-2; kub(24) = 1;     klist(24,:) = [0];
klistnames{25} = 'phi_yA_X';     kc0(25) = 0; klb(25) = 1e-2; kub(25) = 0.5;   klist(25,:) = [0];
klistnames{26} = 'Vmax_yB_XY';   kc0(26) = 0; klb(26) = 1;    kub(26) = 120;   klist(26,:) = [0];
klistnames{27} = 'Km_yB_XY';     kc0(27) = 0; klb(27) = 1;    kub(27) = 1000;  klist(27,:) = [0];

% --- Add the dummy parameter ---------------------------------------------
kmidpoint(end + 1) = 1;
klist(end + 1,:) = 1;
kc0(end + 1) = 0;
klb(end + 1) = 0;
kub(end + 1) = 1;
klistnames{end + 1} = 'dmy';

% --- Choose distribution for parameters ----------------------------------
choosedist = 'unif';

% --- Save set up ---------------------------------------------------------
save([fname,'/efast_setup.mat']);

%% ===== SET UP eFAST ANALYSIS ============================================

% --- list of parameters indecies -----------------------------------------
kidx = find(klist == 1);
lb = klb(kidx);
ub = kub(kidx);
knames = klistnames(kidx);
constrainedkidx = find(kc0(kidx));

% --- Calculate Bonerroni corrected alpha ---------------------------------
alpha = initalpha./length(knames);
    
% --- Build function ------------------------------------------------------
odefun = @(X) efast2MEv1strains(X, kmidpoint, kidx, tmax, odeoptions);

% ---- Simulate model -----------------------------------------------------
nominalX = kmidpoint(kidx);
[performancevector, performancenames, T, yNa, Z] = odefun(nominalX);

% --- Number of outputs from odefun ---------------------------------------
noutputs = length(performancevector);

%% ===== CARRYOUT ANALYSIS ================================================

% --- Calculate the number of parameters ----------------------------------
K = length(kidx);

% --- Carryout GSA using eFAST approach -----------------------------------
[youtput, xparameters, N, OMi, NS] = efastRunAnalysisWithConstraintsV2b(odefun, lb, ub, [], [], choosedist, noutputs, K, NR, NS, MI, constrainedkidx);
    
%% ===== SAVE RESULTS =====================================================
save([fname,'/efast_results.mat'], ...
    'alpha', 'choosedist', 'K', 'kidx', 'kmidpoint', 'knames', 'lb', ...
    'MI', 'N', 'noutputs', 'NR', 'NS', 'odefun', 'odeoptions', 'OMi', ...
    'performancenames', 'performancevector', ...
    'T', 'nominalX', 'tmax', 'ub', ...
    'xparameters', 'Z', 'yNa', 'youtput');
   