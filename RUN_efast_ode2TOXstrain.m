%% ===== SET UP ===========================================================
clear('all'); close('all');

% --- Back code up --------------------------------------------------------
fname = 'EFAST_ode2TOXstrain';
mkdir(fname);
backupcode('RUN_efast_ode2TOXstrain', fname);

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
I_yA_A       = 1e3;
I_yA_B       = 1e3;
% --- Strain B parameters (Genotype A-/B+) --------------------------------
Vmax_gluc_yB = 7.2;  % Maximum glucose uptake rate
Km_gluc_yB   = 5;    % Michaelis constant for glucose uptake
eta_yB       = 1e-3; % Mortality rate constant of strain A
Vmax_yB_A    = 30;   % Maximum B uptake rate by strain A
Km_yB_A      = 50;   % Michaelis constant for B uptake
delta_yB_B   = 1;    % Number of B molecules produced per molecule of glucose is consumed
I_yB_A       = 1e3;
I_yB_B       = 1e3;

% ----- Leaky vectors -----------------------------------------------------
omega(1) = 0.1;      % phi_yA_A
omega(2) = 0.1;      % phi_yB_B

% --- Build vectors -------------------------------------------------------
theta_biomass = [gamma_gluc gamma_A gamma_B];
theta_yA      = [Vmax_gluc_yA Km_gluc_yA eta_yA Vmax_yA_B Km_yA_B delta_yA_A I_yA_A I_yA_B];
theta_yB      = [Vmax_gluc_yB Km_gluc_yB eta_yB Vmax_yB_A Km_yB_A delta_yB_B I_yB_A I_yB_B];

% --- Set up options ------------------------------------------------------
odeoptions = odeset('NonNegative', 1:7);

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
kmidpoint(10:17) = theta_yA;
kmidpoint(18:25) = theta_yB;
kmidpoint(26:27) = omega;

% --- Vector names --------------------------------------------------------
klistnames{ 1} = 'dN0';          kc0( 1) = 0; klb( 1) = 0.01; kub( 1) = 1;     klist( 1,:) = [0];
klistnames{ 2} = 'r0A';          kc0( 2) = 1; klb( 2) = 1e-2; kub( 2) = 1;     klist( 2,:) = [1];
klistnames{ 3} = 'r0B';          kc0( 3) = 1; klb( 3) = 1e-2; kub( 3) = 1;     klist( 3,:) = [1];
klistnames{ 4} = 'G0';           kc0( 4) = 0; klb( 4) = 0;    kub( 4) = 20;    klist( 4,:) = [0];
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
klistnames{16} = 'I_yA_A';       kc0(16) = 0; klb(16) = 1;    kub(16) = 1000;  klist(16,:) = [1];
klistnames{17} = 'I_yA_B';       kc0(17) = 0; klb(17) = 1;    kub(17) = 1000;  klist(17,:) = [1];
% --- Strain B ------------------------------------------------------------
klistnames{18} = 'Vmax_gluc_yB'; kc0(18) = 0; klb(18) = 1;    kub(18) = 30;    klist(18,:) = [1];
klistnames{19} = 'Km_gluc_yB';   kc0(19) = 0; klb(19) = 1;    kub(19) = 100;   klist(19,:) = [0];
klistnames{20} = 'eta_yB';       kc0(20) = 0; klb(20) = 1e-4; kub(20) = 1e-2;  klist(20,:) = [0];
klistnames{21} = 'Vmax_yB_A';    kc0(21) = 0; klb(21) = 1;    kub(21) = 120;   klist(21,:) = [1];
klistnames{22} = 'Km_yB_A';      kc0(22) = 0; klb(22) = 1;    kub(22) = 1000;  klist(22,:) = [0];
klistnames{23} = 'delta_yB_B';   kc0(23) = 0; klb(23) = 1e-2; kub(23) = 1;     klist(23,:) = [0];
klistnames{24} = 'I_yB_A';       kc0(24) = 0; klb(24) = 1;    kub(24) = 1000;  klist(24,:) = [1];
klistnames{25} = 'I_yB_B';       kc0(25) = 0; klb(25) = 1;    kub(25) = 1000;  klist(25,:) = [1];
% --- Leaks ---------------------------------------------------------------
klistnames{26} = 'phi_yA_A';     kc0(26) = 0; klb(26) = 1e-2; kub(26) = 0.5;   klist(26,:) = [1];
klistnames{27} = 'phi_yB_B';     kc0(27) = 0; klb(27) = 1e-2; kub(27) = 0.5;   klist(27,:) = [1];

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
odefun = @(X) efast2TOXstrainsV2b(X, kmidpoint, kidx, tmax, odeoptions);

% ---- Simulate model -----------------------------------------------------
nominalX = kmidpoint(kidx);
[performancevector, performancenames, T, yNa, Y] = odefun(nominalX);

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
    'xparameters', 'Y', 'yNa', 'youtput');