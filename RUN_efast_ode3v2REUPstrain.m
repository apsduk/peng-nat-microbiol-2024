%% ===== SET UP ===========================================================
clear('all'); close('all')

% --- Back code up --------------------------------------------------------
fname = 'EFAST_ode3v2REUPstrain';
mkdir(fname);
backupcode('RUN_efast_ode3v2REUPstrain', fname);

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
aC0 = 0; % mg per L

% --- Biomass parameters --------------------------------------------------
gamma_gluc   = 0.05; % Biomass yield of S. cerevisiae grown on glucose
gamma_A      = 0.5;  % Biomass yield of S. cerevisiae grown on A
gamma_B      = 0.5;  % Biomass yield of S. cerevisiae grown on B
gamma_C      = 0.5;  % Biomass yield of S. cerevisiae grown on C
% --- Strain A parameters ----- (A+/B-/C-) --------------------------------
Vmax_gluc_yA = 7.2;  % Maximum glucose uptake rate
Km_gluc_yA   = 5;    % Michaelis constant for glucose uptake
eta_yA       = 1e-3; % Mortality rate constant of strain A
Vmax_yA_A    = 30;   % Maximum A uptake rate by the strain A
Km_yA_A      = 50;   % Michaelis constant for A uptake
Vmax_yA_B    = 30;   % Maximum B uptake rate by the strain A
Km_yA_B      = 50;   % Michaelis constant for B uptake
Vmax_yA_C    = 30;   % Maximum C uptake rate by the strain A
Km_yA_C      = 50;   % Michaelis constant for C uptake
delta_yA_A   = 1;    % Number of A molecules produced per molecule of glucose is consumed
% --- Strain B parameters ----- (A-/B+/C-) --------------------------------
Vmax_gluc_yB = 7.2;  % Maximum glucose uptake rate
Km_gluc_yB   = 5;    % Michaelis constant for glucose uptake
eta_yB       = 1e-3; % Mortality rate constant of strain B
Vmax_yB_A    = 30;   % Maximum A uptake rate by the strain B
Km_yB_A      = 50;   % Michaelis constant for A uptake
Vmax_yB_B    = 30;   % Maximum B uptake rate by the strain B
Km_yB_B      = 50;   % Michaelis constant for B uptake
Vmax_yB_C    = 30;   % Maximum C uptake rate by the strain B
Km_yB_C      = 50;   % Michaelis constant for C uptake
delta_yB_B   = 1;    % Number of B molecules produced per molecule of glucose is consumed
% --- Strain C parameters ----- (A-/B-/C+) --------------------------------
Vmax_gluc_yC = 7.2;  % Maximum glucose uptake rate
Km_gluc_yC   = 5;    % Michaelis constant for glucose uptake
eta_yC       = 1e-3; % Mortality rate constant of strain C
Vmax_yC_A    = 30;   % Maximum A uptake rate by the strain C
Km_yC_A      = 50;   % Michaelis constant for A uptake
Vmax_yC_B    = 30;   % Maximum A uptake rate by the strain C
Km_yC_B      = 50;   % Michaelis constant for A uptake
Vmax_yC_C    = 30;   % Maximum C uptake rate by the strain C
Km_yC_C      = 50;   % Michaelis constant for C uptake
delta_yC_C   = 1;    % Number of C molecules produced per molecule of glucose is consumed

% ----- Leaky vectors -----------------------------------------------------
phi_yA_A = 0.1;
phi_yB_B = 0.1;
phi_yC_C = 0.1;

% --- Build vectors -------------------------------------------------------
theta_biomass = [gamma_gluc gamma_A gamma_B gamma_C];
theta_yA      = [Vmax_gluc_yA Km_gluc_yA eta_yA Vmax_yA_A Km_yA_A Vmax_yA_B Km_yA_B Vmax_yA_C Km_yA_C delta_yA_A];
theta_yB      = [Vmax_gluc_yB Km_gluc_yB eta_yB Vmax_yB_A Km_yB_A Vmax_yB_B Km_yB_B Vmax_yB_C Km_yB_C delta_yB_B];
theta_yC      = [Vmax_gluc_yC Km_gluc_yC eta_yC Vmax_yC_A Km_yC_A Vmax_yC_B Km_yC_B Vmax_yC_C Km_yC_C delta_yC_C];
omega         = [phi_yA_A phi_yB_B phi_yC_C];

% --- Set up options ------------------------------------------------------
odeoptions = odeset('NonNegative', 1:10);

%% ===== E-FAST SETTINGS ==================================================
% --- e-FAST settings -----------------------------------------------------
initalpha = 0.01;  % p-value threshold
NR = 100;      % number of resampling
NS = 5*257;    % samples per search
MI = 4;        % maximum number of fourier coeffs retained
delta = 0.5;   % lb/ub bounds

% --- Initial conditons ---------------------------------------------------
dN0 = 0.03;
r0A = 1/3;
r0B = 1/3;
r0C = 1/3;

% --- Build k vector ------------------------------------------------------
kmidpoint     = zeros(1,45);
kmidpoint( 1) = dN0; % total initial population
kmidpoint( 2) = r0A;  % A:N ratio
kmidpoint( 3) = r0B;  % B:N ratio
kmidpoint( 4) = r0C;  % C:N ratio
kmidpoint( 5) = G0;  % glucose conc.
kmidpoint( 6) = aA0; % amino acid A conc
kmidpoint( 7) = aB0; % amino acid B conc
kmidpoint( 8) = aC0; % amino acid C conc
kmidpoint( 9:12) = theta_biomass;
kmidpoint(13:22) = theta_yA;
kmidpoint(23:32) = theta_yB;
kmidpoint(33:42) = theta_yC;
kmidpoint(43:45) = omega;

% --- Vector names --------------------------------------------------------
klistnames{ 1} = 'dN0';          kc0( 1) = 0; klb( 1) = 0.01; kub( 1) = 1;     klist( 1,:) = [0];
klistnames{ 2} = 'r0A';          kc0( 2) = 1; klb( 2) = 1e-2; kub( 2) = 1;     klist( 2,:) = [1];
klistnames{ 3} = 'r0B';          kc0( 3) = 1; klb( 3) = 1e-2; kub( 3) = 1;     klist( 3,:) = [1];
klistnames{ 4} = 'r0C';          kc0( 4) = 1; klb( 4) = 1e-2; kub( 4) = 1;     klist( 4,:) = [1];
klistnames{ 5} = 'G0';           kc0( 5) = 0; klb( 5) = 0;    kub( 5) = 20;    klist( 5,:) = [0];
klistnames{ 6} = 'aA0';          kc0( 6) = 0; klb( 6) = 0;    kub( 6) = 75;    klist( 6,:) = [1];
klistnames{ 7} = 'aB0';          kc0( 7) = 0; klb( 7) = 0;    kub( 7) = 75;    klist( 7,:) = [1];
klistnames{ 8} = 'aC0';          kc0( 8) = 0; klb( 8) = 0;    kub( 8) = 75;    klist( 8,:) = [1];
klistnames{ 9} = 'gamma_gluc';   kc0( 9) = 0; klb( 9) = 1e-2; kub( 9) = 1;     klist( 9,:) = [0];
klistnames{10} = 'gamma_A';      kc0(10) = 0; klb(10) = 1e-2; kub(10) = 1;     klist(10,:) = [0];
klistnames{11} = 'gamma_B';      kc0(11) = 0; klb(11) = 1e-2; kub(11) = 1;     klist(11,:) = [0];
klistnames{12} = 'gamma_C';      kc0(12) = 0; klb(12) = 1e-2; kub(12) = 1;     klist(12,:) = [0];
% --- Strain A ------------------------------------------------------------
klistnames{13} = 'Vmax_gluc_yA'; kc0(13) = 0; klb(13) = 1;    kub(13) = 30;    klist(13,:) = [1];
klistnames{14} = 'Km_gluc_yA';   kc0(14) = 0; klb(14) = 1;    kub(14) = 100;   klist(14,:) = [0];
klistnames{15} = 'eta_yA';       kc0(15) = 0; klb(15) = 1e-4; kub(15) = 1e-2;  klist(15,:) = [0];
klistnames{16} = 'Vmax_yA_A';    kc0(16) = 0; klb(16) = 1;    kub(16) = 120;   klist(16,:) = [1];
klistnames{17} = 'Km_yA_A';      kc0(17) = 0; klb(17) = 1;    kub(17) = 1000;  klist(17,:) = [0];
klistnames{18} = 'Vmax_yA_B';    kc0(18) = 0; klb(18) = 1;    kub(18) = 120;   klist(18,:) = [1];
klistnames{19} = 'Km_yA_B';      kc0(19) = 0; klb(19) = 1;    kub(19) = 1000;  klist(19,:) = [0];
klistnames{20} = 'Vmax_yA_C';    kc0(20) = 0; klb(20) = 1;    kub(20) = 120;   klist(20,:) = [1];
klistnames{21} = 'Km_yA_C';      kc0(21) = 0; klb(21) = 1;    kub(21) = 1000;  klist(21,:) = [0];
klistnames{22} = 'delta_yA_A';   kc0(22) = 0; klb(22) = 1e-2; kub(22) = 1;     klist(22,:) = [0];
% --- Strain B ------------------------------------------------------------
klistnames{23} = 'Vmax_gluc_yB'; kc0(23) = 0; klb(23) = 1;    kub(23) = 30;    klist(23,:) = [1];
klistnames{24} = 'Km_gluc_yB';   kc0(24) = 0; klb(24) = 1;    kub(24) = 100;   klist(24,:) = [0];
klistnames{25} = 'eta_yB';       kc0(25) = 0; klb(25) = 1e-4; kub(25) = 1e-2;  klist(25,:) = [0];
klistnames{26} = 'Vmax_yB_A';    kc0(26) = 0; klb(26) = 1;    kub(26) = 120;   klist(26,:) = [1];
klistnames{27} = 'Km_yB_A';      kc0(27) = 0; klb(27) = 1;    kub(27) = 1000;  klist(27,:) = [0];
klistnames{28} = 'Vmax_yB_B';    kc0(28) = 0; klb(28) = 1;    kub(28) = 120;   klist(28,:) = [1];
klistnames{29} = 'Km_yB_B';      kc0(29) = 0; klb(29) = 1;    kub(29) = 1000;  klist(29,:) = [0];
klistnames{30} = 'Vmax_yB_C';    kc0(30) = 0; klb(30) = 1;    kub(30) = 120;   klist(30,:) = [1];
klistnames{31} = 'Km_yB_C';      kc0(31) = 0; klb(31) = 1;    kub(31) = 1000;  klist(31,:) = [0];
klistnames{32} = 'delta_yB_B';   kc0(32) = 0; klb(32) = 1e-2; kub(32) = 1;     klist(32,:) = [0];
% --- Strain B ------------------------------------------------------------
klistnames{33} = 'Vmax_gluc_yC'; kc0(33) = 0; klb(33) = 1;    kub(33) = 30;    klist(33,:) = [1];
klistnames{34} = 'Km_gluc_yC';   kc0(34) = 0; klb(34) = 1;    kub(34) = 100;   klist(34,:) = [0];
klistnames{35} = 'eta_yC';       kc0(35) = 0; klb(35) = 1e-4; kub(35) = 1e-2;  klist(35,:) = [0];
klistnames{36} = 'Vmax_yC_A';    kc0(36) = 0; klb(36) = 1;    kub(36) = 120;   klist(36,:) = [1];
klistnames{37} = 'Km_yC_A';      kc0(37) = 0; klb(37) = 1;    kub(37) = 1000;  klist(37,:) = [0];
klistnames{38} = 'Vmax_yC_B';    kc0(38) = 0; klb(38) = 1;    kub(38) = 120;   klist(38,:) = [1];
klistnames{39} = 'Km_yC_B';      kc0(39) = 0; klb(39) = 1;    kub(39) = 1000;  klist(39,:) = [0];
klistnames{40} = 'Vmax_yC_C';    kc0(40) = 0; klb(40) = 1;    kub(40) = 120;   klist(40,:) = [1];
klistnames{41} = 'Km_yC_C';      kc0(41) = 0; klb(41) = 1;    kub(41) = 1000;  klist(41,:) = [0];
klistnames{42} = 'delta_yC_C';   kc0(42) = 0; klb(42) = 1e-2; kub(42) = 1;     klist(42,:) = [0];
% --- Leaks ---------------------------------------------------------------
klistnames{43} = 'phi_yA_A';     kc0(43) = 0; klb(43) = 1e-2; kub(43) = 0.5;   klist(43,:) = [1];
klistnames{44} = 'phi_yB_B';     kc0(44) = 0; klb(44) = 1e-2; kub(44) = 0.5;   klist(44,:) = [1];
klistnames{45} = 'phi_yC_C';     kc0(45) = 0; klb(45) = 1e-2; kub(45) = 0.5;   klist(45,:) = [1];

% --- Add the dummy parameter ---------------------------------------------
kmidpoint(end + 1) = 1;
klist(end + 1,:) = 1;
klb(end + 1) = 0;
kub(end + 1) = 1;
kc0(end + 1) = 0;
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

% --- Build function --------------------------------------------------
odefun = @(X) efast3v2REUPstrainsV2b(X, kmidpoint, kidx, tmax, odeoptions);

% ---- Simulate model -------------------------------------------------
nominalX = kmidpoint(kidx);
[performancevector, performancenames, T, yNa, Y] = odefun(nominalX);

% --- Number of outputs from odefun -----------------------------------
noutputs = length(performancevector);

% ===== CARRYOUT ANALYSIS =============================================

% --- Calculate the number of parameters ------------------------------
K = length(kidx);

% --- Carryout GSA using eFAST approach -------------------------------
[youtput, xparameters, N, OMi, NS] = efastRunAnalysisWithConstraintsV2b(odefun, lb, ub, [], [], choosedist, noutputs, K, NR, NS, MI, constrainedkidx);

%% ===== SAVE RESULTS =================================================
save([fname,'/efast_results.mat'], ...
    'alpha', 'choosedist', 'K', 'kidx', 'kmidpoint', 'knames', 'lb', ...
    'MI', 'N', 'noutputs', 'NR', 'NS', 'odefun', 'odeoptions', 'OMi', ...
    'performancenames', 'performancevector', ...
    'T', 'nominalX', 'tmax', 'ub', ...
    'xparameters', 'Y', 'yNa', 'youtput');