%% ===== SET UP ===========================================================
clear('all'); close('all');

% --- Choose model --------------------------------------------------------
fname = 'EFAST_ode3v2strain';

%% ===== LOAD RESULTS =====================================================
disp(['Loading results'])
load([fname,'/efast_results.mat'], ...
    'alpha', 'knames', 'MI', 'noutputs', 'OMi', 'youtput', 'performancenames', 'odefun');

%% ===== CALCULATE SENSITIVITY INDICES ====================================
% Si :: first order sensitivity indices
% Sti :: total effect sensitivity indices
disp('Calculating sensitivity indices');
[Si, Sti, rangeSi, rangeSti] = efastSD(youtput, OMi, MI, 1:noutputs);

%% ===== CALCULATE STATISTICS USING SIMPLE T-TEST =========================
disp('Calculating statistics');
stats = efastTTest(Si, rangeSi, Sti, rangeSti, 1:noutputs, alpha);

%% ===== SAVE RESULTS =====================================================
save([fname,'/efast_stats.mat'], ...
    'odefun', 'knames', 'performancenames', 'noutputs', 'rangeSi', 'rangeSti', 'Si', 'Sti', 'stats', 'alpha');
