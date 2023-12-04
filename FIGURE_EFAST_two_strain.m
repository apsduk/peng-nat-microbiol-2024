%% ===== SET UP ===========================================================
clear('all'); close('all'); addpath('lib');

% --- Choose model --------------------------------------------------------
fname = 'RUN_EFAST_ALL_ode2strain';

% --- Back code up --------------------------------------------------------
% backupcode(mfilename, fname);

% --- Significance start offset  ------------------------------------------
staroffset = 0.05;

% --- Save set up ---------------------------------------------------------
load([fname,'/efast_setup.mat'], 'klist');

%% ===== LOAD RESULTS =====================================================
disp(['Loading results ...'])
load([fname,'/efast_stats.mat'], ...
    'knames', 'MI', 'noutputs', 'OMi', 'youtput', 'performancenames', 'odefun', ...
    'rangeSi', 'rangeSti', 'Si', 'Sti', 'originalstats', 'stats', 'alpha');

% --- Extract title ------------------------------------------------------
str = func2str(odefun);
str = split(str,')');
str = split(str{end},'(');
str = str{1};

%% ===== UPDATE KNAMES ====================================================
knames = strrep(knames,'A','1');
knames = strrep(knames,'B','2');
knames = strrep(knames,'phi','\phi');
knames = strrep(knames,'\phi_y1_1','\phi_1');
knames = strrep(knames,'\phi_y2_2','\phi_2');
knames = strrep(knames,'Vmax_gluc_y1','V_{max,G}^{y_1}');
knames = strrep(knames,'Vmax_gluc_y2','V_{max,G}^{y_2}');
knames = strrep(knames,'Vmax_y1_2','V_{max,2}^{y_1}');
knames = strrep(knames,'Vmax_y2_1','V_{max,1}^{y_2}');
knames = strrep(knames,'a10','x_{0,1}');
knames = strrep(knames,'a20','x_{0,2}');
knames = strrep(knames,'r01','r_{0,1}');
knames = strrep(knames,'r02','r_{0,2}');
knames = strrep(knames,'dmy','\delta');
knames

%% ===== PLOT OUTPUT ======================================================
fsort = figure('Units', 'centimeters', 'Position', [5 5 20 12], 'Visible', 'off');

n_k = length(knames);

outputname = {'Population', 'Batch time', 'y_1/(y_1 + y_2)', 'J_{grow}^{y_1}'};
v_y = [1 2 3 5];

ymax = [0.6 0.85 0.6 0.8]; % length(v_y);

% --- Plot N --------------------------------------------------------------
for yid = 1:length(v_y)
    
    y = v_y(yid);
    
    performancenames{y}
    
    % --- Test stats ------------------------------------------------------
    sig_Si  = find(stats.p_Si(:,y) < alpha);
    sig_Sti = find(stats.p_Sti(:,y) < alpha);

    % --- Order points by Si ----------------------------------------------
    [~, idx] = sort(stats.avg_Si(:,y), 'descend');
    for i = 1:length(idx); sortlist{i} = knames{idx(i)}; end

    % --- Find significant points -----------------------------------------
    [~, sig_Si_loc] = ismember(sig_Si, idx);
    [~, sig_Sti_loc] = ismember(sig_Sti, idx);

    % --- Get list of all points ------------------------------------------
    sig_loc = unique([sig_Si_loc; sig_Sti_loc]);

    % --- Plot outputs ----------------------------------------------------
    figure(fsort.Number);
    subplot(2, 2, yid);
    hold('on');
    bar(1:n_k, stats.avg_Sti(idx, y), 'EdgeColor', 'none');
    bar(1:n_k, stats.avg_Si(idx, y), 'EdgeColor', 'none');
    errorbar(1:n_k, stats.avg_Sti(idx, y), stats.std_Sti(idx, y), 'k.', 'LineWidth', 1);
    errorbar(1:n_k, stats.avg_Si(idx, y),  stats.std_Si(idx, y),  'k.', 'LineWidth', 1);
    plot(sig_loc, stats.avg_Sti(idx(sig_loc), y) + stats.std_Sti(idx(sig_loc), y) + staroffset, '*k', 'MarkerSize', 6);
    legend('total-order','first-order')
    set(gca, 'XTick', 1:n_k, 'XTickLabel', sortlist, 'Box', 'on', 'FontSize', 10, 'FontWeight', 'bold', 'LineWidth', 1);
    ylabel([outputname{yid},char(10),'eFAST sensitivity'], 'FontSize', 10);
    xtickangle(-45);
    ylim([0 ymax(yid)]);

end
saveas(fsort.Number, 'FIG_1_E_H', 'svg');

