%% ===== SET UP ===========================================================
clear('all'); close('all');

% --- Choose model --------------------------------------------------------
fname = 'EFAST_ode3v2REUPstrain';

% --- Save set up ---------------------------------------------------------
load([fname,'/efast_setup.mat'], 'klist');

%% ===== LOAD RESULTS =====================================================
disp(['Loading results'])
load([fname,'/efast_stats.mat'], ...
    'alpha', 'odefun', 'knames', 'MI', 'noutputs', 'OMi', 'youtput', 'performancenames', 'odefun', ...
    'rangeSi', 'rangeSti', 'stats', 'Si', 'Sti');

% --- Extract title ------------------------------------------------------
str = func2str(odefun);
str = split(str,')');
str = split(str{end},'(');
str = str{1};

%% ===== UPDATE KNAMES ====================================================
knames = strrep(knames,'aA0', 'x_{0,1}');
knames = strrep(knames,'aB0', 'x_{0,2}');
knames = strrep(knames,'aC0', 'x_{0,3}');
knames = strrep(knames,'delta_yA_A', '\delta_1');
knames = strrep(knames,'delta_yA_C', '\delta^{y_1}_3');
knames = strrep(knames,'delta_yB_B', '\delta_2');
knames = strrep(knames,'delta_yC_C', '\delta_3');
knames = strrep(knames,'dN0', 'N_0');
knames = strrep(knames,'eta_yA', '\eta_{y_1}');
knames = strrep(knames,'eta_yB', '\eta_{y_2}');
knames = strrep(knames,'eta_yC', '\eta_{y_3}');
knames = strrep(knames,'G', 'G_0');
knames = strrep(knames,'gamma_A', '\gamma_1');
knames = strrep(knames,'gamma_B', '\gamma_2');
knames = strrep(knames,'gamma_C', '\gamma_3');
knames = strrep(knames,'gamma_gluc', '\gamma_{gluc}');
knames = strrep(knames,'Km_gluc_yA', 'K_{M,G}^{y_1}');
knames = strrep(knames,'Km_gluc_yB', 'K_{M,G}^{y_2}');
knames = strrep(knames,'Km_gluc_yC', 'K_{M,G}^{y_3}');
knames = strrep(knames,'Km_yA_B', 'K_{M,2}^{y_1}');
knames = strrep(knames,'Km_yA_C', 'K_{M,3}^{y_1}');
knames = strrep(knames,'Km_yB_A', 'K_{M,1}^{y_2}');
knames = strrep(knames,'Km_yB_C', 'K_{M,3}^{y_2}');
knames = strrep(knames,'Km_yC_A', 'K_{M,1}^{y_3}');
knames = strrep(knames,'Km_yC_B', 'K_{M,2}^{y_3}');
knames = strrep(knames,'phi_yA_A', '\phi_1');
knames = strrep(knames,'phi_yA_C', '\phi_3^{y_1}');
knames = strrep(knames,'phi_yB_B', '\phi_2');
knames = strrep(knames,'phi_yC_C', '\phi_3');
knames = strrep(knames,'r0A', 'r_{0,1}');
knames = strrep(knames,'r0B', 'r_{0,2}');
knames = strrep(knames,'r0C', 'r_{0,3}');
knames = strrep(knames,'Vmax_gluc_yA', 'V_{max,G}^{y_1}');
knames = strrep(knames,'Vmax_gluc_yB', 'V_{max,G}^{y_2}');
knames = strrep(knames,'Vmax_gluc_yC', 'V_{max,G}^{y_3}');
knames = strrep(knames,'Vmax_yA_B', 'V_{max,2}^{y_1}');
knames = strrep(knames,'Vmax_yA_C', 'V_{max,3}^{y_1}');
knames = strrep(knames,'Vmax_yB_A', 'V_{max,1}^{y_2}');
knames = strrep(knames,'Vmax_yB_C', 'V_{max,3}^{y_2}');
knames = strrep(knames,'Vmax_yC_A', 'V_{max,1}^{y_3}');
knames = strrep(knames,'Vmax_yC_B', 'V_{max,2}^{y_3}');
knames = strrep(knames,'dmy','\delta_{efast}');

performancenames = strrep(performancenames,'max(J_grow_yA)', 'y_1 growth');
performancenames = strrep(performancenames,'max(J_grow_yB)', 'y_2 growth');
performancenames = strrep(performancenames,'max(J_grow_yC)', 'y_3 growth');
performancenames = strrep(performancenames,'max(J_leak_yA_A)', 'x_1 production');
performancenames = strrep(performancenames,'max(J_leak_yB_B)', 'x_2 production');
performancenames = strrep(performancenames,'max(J_leak_yC_C)', 'x_3 production');
performancenames = strrep(performancenames,'max(J_upt_yX_A)', 'x_1 uptake');
performancenames = strrep(performancenames,'max(J_upt_yX_B)', 'x_2 uptake');
performancenames = strrep(performancenames,'max(J_upt_yX_C)', 'x_3 uptake');
performancenames = strrep(performancenames,'N(end)', 'Population size');
performancenames = strrep(performancenames,'T(tdx)', 'Batch time');
performancenames = strrep(performancenames,'yAa(end)_by_yNa(end)', 'y_1 ratio');
performancenames = strrep(performancenames,'yBa(end)_by_yNa(end)', 'y_2 ratio');
performancenames = strrep(performancenames,'yCa(end)_by_yNa(end)', 'y_3 ratio');


%% ===== PLOT OUTPUT ======================================================
fsort = figure('Units', 'centimeters', 'Position', [5 5 30 30], 'Visible', 'off');

n_k = length(knames);

v_fid = [1 2, 4 5 6, 7 8 9, 10 11 12, 13 14 15];
v_y   = [1 2, 3 4 5, 6 7 8,  9 10 11, 12 13 14];

yoffset = 0.05;
cmap = lines(8);

% --- Plot N --------------------------------------------------------------
for yid = 1:length(v_y)
    
    y = v_y(yid);
    
    performancenames{y}
    
    [~, idx] = sort(stats.avg_Si(:,y), 'descend');
    for i = 1:length(idx); sortlist{i} = knames{idx(i)}; end
    
    % --- Test stats -------------------------------------------------
    sig_Si  = find(stats.p_Si(:,y) < alpha);
    sig_Sti = find(stats.p_Sti(:,y) < alpha);
    
    % --- Find significant points ------------------------------------
    [~, sig_Si_loc] = ismember(sig_Si, idx);
    [~, sig_Sti_loc] = ismember(sig_Sti, idx);
    
    figure(fsort.Number);
    subplot(5, 3, v_fid(yid));
    hold('on');
    bar(1:n_k, stats.avg_Sti(idx, y), 'EdgeColor', 'none');
    bar(1:n_k, stats.avg_Si(idx, y), 'EdgeColor', 'none');
    errorbar(1:n_k, stats.avg_Sti(idx, y), stats.std_Sti(idx, y), 'k.', 'LineWidth', 1);
    errorbar(1:n_k, stats.avg_Si(idx, y),  stats.std_Si(idx, y),  'k.', 'LineWidth', 1);
    set(gca, 'XTick', 1:n_k, 'XTickLabel', sortlist, 'Box', 'on', 'FontSize', 8, 'LineWidth', 1);
    xtickangle(-90);
    ylabel([performancenames{yid},char(10),' eFAST sensitivity'], 'FontSize', 10);
    plot( sig_Si_loc, stats.avg_Sti(idx( sig_Si_loc), y) + stats.std_Sti( idx(sig_Si_loc), y) + yoffset, '+', 'Color', cmap(1,:),  'MarkerSize', 6);
    plot(sig_Sti_loc, stats.avg_Sti(idx(sig_Sti_loc), y) + stats.std_Sti(idx(sig_Sti_loc), y) + yoffset, 'x', 'Color', cmap(2,:), 'MarkerSize', 6);

    [~,ymaxid] = max(stats.avg_Sti(:, y));
    ymax = stats.avg_Sti(ymaxid, y) + stats.std_Sti(ymaxid, y) + yoffset;
    ylim([0, 1.1*round(ymax,2)]);

    legend('total-order','first-order')

end
saveas(fsort.Number, fname, 'svg');


