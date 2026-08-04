x_len_all       = [];
y_vol_len_all   = [];
y_flat_len_all  = [];

for li = 1:n_L
    for di = 1:n_D
        for ti = 1:n_trials
            val_v = vol_far(li, di, ti);
            val_f = flat_far(li, di, ti);
            if ~isnan(val_v)
                x_len_all      = [x_len_all;      lengths_mm(li)];
                y_vol_len_all  = [y_vol_len_all;   val_v];
                y_flat_len_all = [y_flat_len_all;  val_f];
            end
        end
    end
end

% For vs-DIAMETER plots: pool all lengths and all trials
x_dia_all       = [];
y_vol_dia_all   = [];
y_flat_dia_all  = [];

for li = 1:n_L
    for di = 1:n_D
        for ti = 1:n_trials
            val_v = vol_far(li, di, ti);
            val_f = flat_far(li, di, ti);
            if ~isnan(val_v)
                x_dia_all      = [x_dia_all;      diameters_mm(di)];
                y_vol_dia_all  = [y_vol_dia_all;   val_v];
                y_flat_dia_all = [y_flat_dia_all;  val_f];
            end
        end
    end
end

%% ======================== HELPER: FIT AND PLOT ===========================
function plot_fit_with_CI(x, y, cap_val, cap_label, x_label, y_label, fig_title, fit_type_str, start_point)
    % Fits data, shades 95% prediction interval, plots all data points.

    figure('Position', [100 100 700 500]);
    hold on; grid on;

    % --- Curve fit ---
    ft = fittype(fit_type_str);
    linear_fits = {'poly1', 'poly2', 'poly3', 'poly4'};
    if any(strcmp(fit_type_str, linear_fits))
        [fit_obj, gof, ~] = fit(x, y, ft);
    else
        [fit_obj, gof, ~] = fit(x, y, ft, 'StartPoint', start_point);
    end

    % Significance check on parameters
    param_CI = confint(fit_obj, 0.95);
    significance_check = (param_CI(1,:) < 0 & param_CI(2,:) < 0) | ...
                         (param_CI(1,:) > 0 & param_CI(2,:) > 0);

    % Fine x vector for smooth curves
    x_ = linspace(min(x), max(x), 200)';

    % 95% functional prediction interval
    y_predint = predint(fit_obj, x_, 0.95, 'functional', 'off');
    y_upper   = y_predint(:, 2);
    y_lower   = y_predint(:, 1);

    % --- Shade confidence interval (plot first so it sits behind everything) ---
    x2         = [x_', fliplr(x_')];
    in_between = [y_upper; flipud(y_lower)]';
    hf = fill(x2, in_between, 'b', 'DisplayName', '95% prediction interval');
    hf.FaceAlpha = 0.2;
    hf.FaceColor = [0.2 0.4 0.8];
    hf.EdgeColor = 'none';

    % --- Plot fit line ---
    if min(significance_check > 0)
        plot(x_, feval(fit_obj, x_), '-', 'Color', [0.2 0.4 0.8], ...
            'LineWidth', 2, 'DisplayName', 'LSQ fit');
    else
        plot(x_, feval(fit_obj, x_), 'r--', ...
            'LineWidth', 2, 'DisplayName', 'Fit not statistically significant');
    end

    % --- Raw data points on top ---
    plot(x, y, 'o', 'MarkerFaceColor', [0.1 0.1 0.35], ...
        'MarkerEdgeColor', [0.1 0.1 0.35], 'MarkerSize', 5, ...
        'DisplayName', 'Data (all trials & configurations)');

    % --- Capped reference line ---
    yline(cap_val, '--k', 'LineWidth', 1.5, ...
        'Label', cap_label, 'LabelHorizontalAlignment', 'left');

    xlabel(x_label);
    ylabel(y_label);
    title(fig_title);

    % Get fit formula and coefficients as a string
    fit_str = formula(fit_obj);
    param_names = coeffnames(fit_obj);
    param_vals = coeffvalues(fit_obj);
    for k = 1:length(param_names)
        fit_str = strrep(fit_str, param_names{k}, sprintf('%.4g', param_vals(k)));
    end
    
    % Add R-squared
    gof_text = sprintf('R² = %.4f', gof.rsquare);
    
    % Display on figure
    annotation_str = sprintf('Fit: y = %s\n%s', fit_str, gof_text);
    text(0.05, 0.95, annotation_str, ...
        'Units', 'normalized', ...
        'VerticalAlignment', 'top', ...
        'FontSize', 9, ...
        'BackgroundColor', 'white', ...
        'EdgeColor', 'black');
   
    legend('Location', 'best');

    % Print fit summary
    fprintf('\n=== %s ===\n', fig_title);
    disp(fit_obj);
    fprintf('Parameter 95%% CI:\n'); disp(param_CI);
    fprintf('Statistically significant parameters: '); disp(significance_check);
end

%% ======================== FIT TYPE SETTINGS ==============================
% Change these strings to try different models.
% Common options:
%   'poly1'   — linear:     a*x + b
%   'poly2'   — quadratic:  a*x^2 + b*x + c
%   'power1'  — power:      a*x^b
%   'exp1'    — exponential: a*exp(b*x)
%
% Matching start points:
%   poly1  -> [1, 1]
%   poly2  -> [1, 1, 1]
%   power1 -> [1, 1]
%   exp1   -> [1, 0.01]
%
% To explore interactively first, uncomment the cftool lines at the bottom.

fit_vol_len  = 'poly1';  sp_vol_len  = [1, 1];
fit_flat_len = 'poly1';  sp_flat_len = [1, 1];
fit_vol_dia  = 'poly1';  sp_vol_dia  = [1, 1];
fit_flat_dia = 'poly1';  sp_flat_dia = [1, 1];

%% ======================== GENERATE PLOTS =================================

plot_fit_with_CI(...
    x_len_all, y_vol_len_all, ...
    cap_vol_far, sprintf('Capped (%.1f dB)', cap_vol_far), ...
    'Port Length (mm)', 'Average Volume (dB)', ...
    'Volume vs. Port Length — Far Mic (1 m)', ...
    fit_vol_len, sp_vol_len);

plot_fit_with_CI(...
    x_len_all, y_flat_len_all, ...
    cap_flat_far, sprintf('Capped (%.2f dB)', cap_flat_far), ...
    'Port Length (mm)', 'Flatness — RMS Deviation (dB)', ...
    'Flatness vs. Port Length — Far Mic (1 m)', ...
    fit_flat_len, sp_flat_len);

plot_fit_with_CI(...
    x_dia_all, y_vol_dia_all, ...
    cap_vol_far, sprintf('Capped (%.1f dB)', cap_vol_far), ...
    'Port Diameter (mm)', 'Average Volume (dB)', ...
    'Volume vs. Port Diameter — Far Mic (1 m)', ...
    fit_vol_dia, sp_vol_dia);

plot_fit_with_CI(...
    x_dia_all, y_flat_dia_all, ...
    cap_flat_far, sprintf('Capped (%.2f dB)', cap_flat_far), ...
    'Port Diameter (mm)', 'Flatness — RMS Deviation (dB)', ...
    'Flatness vs. Port Diameter — Far Mic (1 m)', ...
    fit_flat_dia, sp_flat_dia);

%% ======================== OPEN CFTOOL MANUALLY ===========================
% Uncomment to explore fit types interactively before committing above:
% cftool(x_len_all, y_vol_len_all)
% cftool(x_len_all, y_flat_len_all)
% cftool(x_dia_all, y_vol_dia_all)
% cftool(x_dia_all, y_flat_dia_all)
