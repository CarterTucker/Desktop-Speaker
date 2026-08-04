% For each trial
% - loads the time-domain mic data
% - loads the input sweep .wav and resamples to match mic
% - computes the frequency response from deconvolution
% - convert Bode to dB
% - compute average volume (dB) and flatness (RMS deviation from mean)
% over min and max frequencies (adjust according to mic specs)
% - plot:
%  -> Bode plots for each trial with the independent axis being length or
%  diameter (level 1)
%  -> Volume/Flatness vs. Length/Diameter (4 total) w/ error bars (level 2/3) 
%
% Filenames:
% - length_(lenght in mm)_dia_(diameter as 10xmm)_trial(#).csv
%   (45 mm long, 17.5 mm diameter would be length_45_dia_175_trial1.csv)
%
% CSV columns
% Latest: Time (s), Latest: Mic 1 (2 cm), Col 3: Latest: Offset Pressure 1,
% Col 4: Latest: Mic 2 (1 m), Col 5: Latest: Offset Pressure 2

% SET PARAMETERS ---------------------------------------------------------
csv_dir  = 'CSVs';
wav_file = 'sweep_sine.wav'; 

Fs_mic         = 10000;
wav_Fs_nominal = 48000;

lengths_mm   = [0, 45, 90, 135, 180];
diameters_mm = [17.5, 35.0, 52.5, 70.0];

n_trials = 5;

f_min = 126; 
f_max = 500;

smooth_frac = 0.05; 
reg_eps     = 1e-6;

% LOAD THE WAV FILE -------------------------------------------------------
disp('Loading sweep .wav...');
[x_wav, Fs_wav] = audioread(wav_file);

resampled_wav = resample(x_wav, Fs_mic, Fs_wav);
resampled_wav = resampled_wav(:);

fprintf('Sweep loaded: %.1f s at %d Hz, resampled from %d Hz\n', ...
    length(x_wav)/Fs_wav, Fs_mic, Fs_wav);

% HELPER FUNCTIONS --------------------------------------------------------
    function [f, H_dB] = compute_response(y, resampled_wav, Fs, f_min, f_max, reg_eps, smooth_frac)
        N = length(y);
        x = resampled_wav;
        if length(x) < N
            x = [x; zeros(N - length(x), 1)];
        else
            x = x(1:N);
        end

        output_FFT = fft(y, N);
        input_FFT  = fft(x, N);

        H = output_FFT .* conj(input_FFT) ./ (abs(input_FFT).^2 + reg_eps * max(abs(input_FFT).^2));

        f_full = (0:N-1)' * Fs / N;
        H_full = abs(H(1:N));

        idx_ny  = floor(N/2);
        f       = f_full(1:idx_ny);
        H_mag   = H_full(1:idx_ny);

        win      = max(1, round(smooth_frac * idx_ny));
        H_smooth = movmean(H_mag, win);

        H_dB = 20 * log10(H_smooth + eps) + 30;

        band = f >= f_min & f <= f_max;
        f    = f(band);
        H_dB = H_dB(band);
    end

    function [avg_vol, flatness] = compute_metrics(H_dB)
        avg_vol  = mean(H_dB);
        flatness = sqrt(mean((H_dB - avg_vol).^2));
    end

% PROCESS THE CSVs --------------------------------------------------------
fprintf('Processing CSV files...\n');

n_L = length(lengths_mm);
n_D = length(diameters_mm);

vol_close  = nan(n_L, n_D, n_trials);
vol_far    = nan(n_L, n_D, n_trials);
flat_close = nan(n_L, n_D, n_trials);
flat_far   = nan(n_L, n_D, n_trials);

bode_store_close = cell(n_L, n_D);
bode_store_far   = cell(n_L, n_D);
bode_f           = [];
bode_f_store     = cell(n_L, n_D);

vol_close_capped  = nan(1, n_trials);
vol_far_capped    = nan(1, n_trials);
flat_close_capped = nan(1, n_trials);
flat_far_capped   = nan(1, n_trials);

% Process ported configurations
for li = 1:n_L
    for di = 1:n_D
        for ti = 1:n_trials
            fname = sprintf('length_%d_dia_%d_trial%d.csv', ...
                lengths_mm(li), 10*diameters_mm(di), ti);
            fpath = fullfile(csv_dir, fname);

            if ~isfile(fpath)
                warning('File not found: %s', fpath);
                continue;
            end

            try
                T = readtable(fpath, 'NumHeaderLines', 1, ...
                    'VariableNamingRule', 'preserve');
            catch
                warning('Could not read: %s', fpath);
                continue;
            end

            data    = T{:,:}; 
            y_close = data(:,2);
            y_far   = data(:,4);

            [f_band, H_close] = compute_response(y_close, resampled_wav, Fs_mic, f_min, f_max, reg_eps, smooth_frac);
            [~,      H_far]   = compute_response(y_far,   resampled_wav, Fs_mic, f_min, f_max, reg_eps, smooth_frac);

            if isempty(bode_f)
                bode_f = f_band;
            end

            [vol_close(li,di,ti),  flat_close(li,di,ti)]  = compute_metrics(H_close);
            [vol_far(li,di,ti),    flat_far(li,di,ti)]    = compute_metrics(H_far);

            if ti == 1
                bode_store_close{li,di} = H_close;
                bode_store_far{li,di}   = H_far;
                bode_f_store{li,di}     = f_band;
            end

            fprintf('  Processed: %s\n', fname);
        end
    end
end

% Process capped (sealed) condition
for ti = 1:n_trials
    fname = sprintf('capped_trial%d.csv', ti);
    fpath = fullfile(csv_dir, fname);

    if ~isfile(fpath)
        warning('File not found: %s', fpath);
        continue;
    end

    try
        T = readtable(fpath, 'NumHeaderLines', 1, ...
            'VariableNamingRule', 'preserve');
    catch
        warning('Could not read: %s', fpath);
        continue;
    end

    data    = T{:,:};
    y_close = data(:,2);
    y_far   = data(:,4);

    [f_band, H_close] = compute_response(y_close, resampled_wav, Fs_mic, f_min, f_max, reg_eps, smooth_frac);
    [~,      H_far]   = compute_response(y_far,   resampled_wav, Fs_mic, f_min, f_max, reg_eps, smooth_frac);

    if isempty(bode_f)
        bode_f = f_band;
    end

    [vol_close_capped(ti),  flat_close_capped(ti)]  = compute_metrics(H_close);
    [vol_far_capped(ti),    flat_far_capped(ti)]    = compute_metrics(H_far);

    fprintf('  Processed capped trial %d\n', ti);
end

% COMPUTE MEANS AND STANDARD DEVIATIONS -----------------------------------
vol_close_mean  = mean(vol_close,  3, 'omitnan');
vol_far_mean    = mean(vol_far,    3, 'omitnan');
flat_close_mean = mean(flat_close, 3, 'omitnan');
flat_far_mean   = mean(flat_far,   3, 'omitnan');

vol_close_std   = std(vol_close,  0, 3, 'omitnan');
vol_far_std     = std(vol_far,    0, 3, 'omitnan');
flat_close_std  = std(flat_close, 0, 3, 'omitnan');
flat_far_std    = std(flat_far,   0, 3, 'omitnan');

% Capped reference values with uncertainty (95% CI = 2*std/sqrt(n))
cap_vol_far    = mean(vol_far_capped,    'omitnan');
cap_flat_far   = mean(flat_far_capped,   'omitnan');
cap_vol_far_unc   = 2 * std(vol_far_capped,   'omitnan') / sqrt(sum(~isnan(vol_far_capped)));
cap_flat_far_unc  = 2 * std(flat_far_capped,  'omitnan') / sqrt(sum(~isnan(flat_far_capped)));

cap_vol_close  = mean(vol_close_capped,  'omitnan');
cap_flat_close = mean(flat_close_capped, 'omitnan');
cap_vol_close_unc  = 2 * std(vol_close_capped,  'omitnan') / sqrt(sum(~isnan(vol_close_capped)));
cap_flat_close_unc = 2 * std(flat_close_capped, 'omitnan') / sqrt(sum(~isnan(flat_close_capped)));

% COLORS AND LABELS -------------------------------------------------------
red_shades    = [1.00 0.70 0.70; 0.90 0.40 0.40; 0.75 0.15 0.15; 0.55 0.00 0.00];
purple_shades = [0.80 0.70 0.90; 0.65 0.45 0.80; 0.50 0.25 0.70; 0.38 0.10 0.58; 0.25 0.00 0.45];

dia_colors = red_shades;
len_colors = purple_shades;

dia_labels = arrayfun(@(d) sprintf('%.1f mm', d), diameters_mm, 'UniformOutput', false);
len_labels = arrayfun(@(l) sprintf('%d mm', l),   lengths_mm,   'UniformOutput', false);

% REPRESENTATIVE BODE PLOT WITH AVERAGE VOLUME LINE -----------------------
rep_li = 4; % length index  (4 = 135 mm)
rep_di = 2; % diameter index (2 = 35.0 mm)
rep_ti = 1; % trial number

fname = sprintf('length_%d_dia_%d_trial%d.csv', ...
    lengths_mm(rep_li), 10*diameters_mm(rep_di), rep_ti);
fpath = fullfile(csv_dir, fname);

T     = readtable(fpath, 'NumHeaderLines', 1, 'VariableNamingRule', 'preserve');
data  = T{:,:};
y_far = data(:,4);

[f_rep, H_rep] = compute_response(y_far, resampled_wav, Fs_mic, ...
    f_min, f_max, reg_eps, smooth_frac);

avg = mean(H_rep);

figure('Position', [100 100 700 500]);
hold on; grid on;
set(gca, 'XScale', 'log');

plot(f_rep, H_rep, 'Color', [0.55 0.00 0.00], 'LineWidth', 1.5, ...
    'DisplayName', sprintf('L=%d mm, D=%.1f mm', lengths_mm(rep_li), diameters_mm(rep_di)));

yline(avg, '--k', 'LineWidth', 1.5, ...
    'Label', sprintf('\\bar{L} = %.1f dB', avg), ...
    'LabelHorizontalAlignment', 'right', ...
    'FontSize', 12, 'FontName', 'Arial');

xlabel('Frequency (Hz)', 'FontSize', 14, 'FontName', 'Arial');
ylabel('Response (dB)',  'FontSize', 14, 'FontName', 'Arial');
title(sprintf('Frequency Response — L=%d mm, D=%.1f mm', ...
    lengths_mm(rep_li), diameters_mm(rep_di)), ...
    'FontSize', 16, 'FontName', 'Arial');
legend('Location', 'best', 'FontSize', 12);
ax = gca; ax.FontSize = 12; ax.FontName = 'Arial'; ax.LineWidth = 1.2;

% BODE PLOTS — one figure per diameter, lines vary over length ------------
for di = 1:n_D
    figure('Name', sprintf('Bode - Diameter %.1f mm', diameters_mm(di)), ...
        'Position', [100 100 700 500]);
    hold on; grid on;
    set(gca, 'XScale', 'log');
    xlabel('Frequency (Hz)', 'FontSize', 14, 'FontName', 'Arial');
    ylabel('Response (dB)',  'FontSize', 14, 'FontName', 'Arial');
    title(sprintf('Frequency Response — Port Diameter %.1f mm', diameters_mm(di)), ...
        'FontSize', 16, 'FontName', 'Arial');

    for li = 1:n_L
        if ~isempty(bode_store_far{li,di})
            plot(bode_f_store{li,di}, bode_store_far{li,di}, ...
                'Color', len_colors(li,:), ...
                'LineWidth', 1.5, 'DisplayName', len_labels{li});
        end
    end
    legend('Location', 'best');
    ax = gca; ax.FontSize = 12; ax.FontName = 'Arial'; ax.LineWidth = 1.2;
end

% COMPARISON BODE PLOT: SHORTEST vs LONGEST PORT FOR ONE DIAMETER
% Shows average volume as dashed line and RMS deviation as shaded band
% Change rep_di to select which diameter to show

rep_di    = 4;   % diameter index (2 = 35.0 mm)
li_short  = 1;   % shortest length index (1 = 0 mm)  -- change if excluding 0mm
li_long   = 5;   % longest length index  (5 = 180 mm)

% Colors for short and long
col_short = [0.20 0.45 0.70];   % blue for short
col_long  = [0.75 0.15 0.15];   % red for long

% --- Load and compute short port ---
fname_s = sprintf('length_%d_dia_%d_trial%d.csv', ...
    lengths_mm(li_short), 10*diameters_mm(rep_di), 1);
T_s     = readtable(fullfile(csv_dir, fname_s), ...
    'NumHeaderLines', 1, 'VariableNamingRule', 'preserve');
data_s  = T_s{:,:};
[f_s, H_s] = compute_response(data_s(:,4), resampled_wav, Fs_mic, ...
    f_min, f_max, reg_eps, smooth_frac);

avg_s  = mean(H_s);
flat_s = sqrt(mean((H_s - avg_s).^2));

% --- Load and compute long port ---
fname_l = sprintf('length_%d_dia_%d_trial%d.csv', ...
    lengths_mm(li_long), 10*diameters_mm(rep_di), 1);
T_l     = readtable(fullfile(csv_dir, fname_l), ...
    'NumHeaderLines', 1, 'VariableNamingRule', 'preserve');
data_l  = T_l{:,:};
[f_l, H_l] = compute_response(data_l(:,4), resampled_wav, Fs_mic, ...
    f_min, f_max, reg_eps, smooth_frac);

avg_l  = mean(H_l);
flat_l = sqrt(mean((H_l - avg_l).^2));

% --- Plot ---
figure('Position', [100 100 800 550]);
hold on; grid on;
set(gca, 'XScale', 'log');

% RMS deviation shading for short port (behind everything)
x_fill = [f_s', fliplr(f_s')];
y_fill = [(avg_s + flat_s) * ones(1, length(f_s)), ...
           (avg_s - flat_s) * ones(1, length(f_s))];
hfs = fill(x_fill, y_fill, col_short);
hfs.FaceAlpha      = 0.12;
hfs.EdgeColor      = 'none';
hfs.HandleVisibility = 'off';

% RMS deviation shading for long port (behind everything)
x_fill = [f_l', fliplr(f_l')];
y_fill = [(avg_l + flat_l) * ones(1, length(f_l)), ...
           (avg_l - flat_l) * ones(1, length(f_l))];
hfl = fill(x_fill, y_fill, col_long);
hfl.FaceAlpha      = 0.12;
hfl.EdgeColor      = 'none';
hfl.HandleVisibility = 'off';

% Response curves
plot(f_s, H_s, 'Color', col_short, 'LineWidth', 1.8, ...
    'DisplayName', sprintf('L = %d mm', lengths_mm(li_short)));
plot(f_l, H_l, 'Color', col_long,  'LineWidth', 1.8, ...
    'DisplayName', sprintf('L = %d mm', lengths_mm(li_long)));

% Relative volume lines
yline(avg_s, '--', 'Color', col_short, 'LineWidth', 1.5, ...
    'Label', sprintf('Avg. Vol. = %.1f dB,  \\sigma = %.2f dB', avg_s, flat_s), ...
    'LabelHorizontalAlignment', 'right', ...
    'FontSize', 11, 'FontName', 'Arial', ...
    'HandleVisibility', 'off');

yline(avg_l, '--', 'Color', col_long, 'LineWidth', 1.5, ...
    'Label', sprintf('Avg. Vol. = %.1f dB,  \\sigma = %.2f dB', avg_l, flat_l), ...
    'LabelHorizontalAlignment', 'right', ...
    'FontSize', 11, 'FontName', 'Arial', ...
    'HandleVisibility', 'off');

xlabel('Frequency (Hz)', 'FontSize', 14, 'FontName', 'Arial');
ylabel('Response (dB)',  'FontSize', 14, 'FontName', 'Arial');
title(sprintf('Frequency Response — D = %.1f mm, Short vs. Long Port', ...
    diameters_mm(rep_di)), 'FontSize', 16, 'FontName', 'Arial');
legend('Location', 'best', 'FontSize', 12, 'FontName', 'Arial');

ax = gca;
ax.FontSize = 12; ax.FontName = 'Arial'; ax.LineWidth = 1.2;

% =========================================================================
% LEVEL 3 DATA POOLING — PORTED ONLY (length = 0 excluded as sealed proxy)
% Zero-length ports with nonzero diameter are open holes, not true ports,
% so they are excluded from the trend fits along with the fully capped case.
% =========================================================================

% lengths_mm index 1 corresponds to 0 mm — skip it for trend pooling
x_len_all      = [];
y_vol_len_all  = [];
y_flat_len_all = [];

for li = 1:n_L          % start at 2 to skip 0 mm length
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

x_dia_all      = [];
y_vol_dia_all  = [];
y_flat_dia_all = [];

for li = 1:n_L          % start at 2 to skip 0 mm length
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

% =========================================================================
% HELPER: fit and plot with 95% prediction interval
% No R^2 in legend. Equation shown as y = p1*x + p2 with CI on parameters.
% Sealed reference shown as dashed line with uncertainty band.
% =========================================================================
function plot_fit_with_CI(x, y, cap_val, cap_unc, x_label, y_label, fig_title, ...
    fit_type_str, start_point, line_color, x_lim)

    figure('Position', [100 100 800 580]);
    hold on; grid on;

    % Set x axis limits
    xlim(x_lim);

    % Curve fit
    ft = fittype(fit_type_str);
    linear_fits = {'poly1', 'poly2', 'poly3', 'poly4'};
    if any(strcmp(fit_type_str, linear_fits))
        [fit_obj, ~, ~] = fit(x, y, ft);
    else
        [fit_obj, ~, ~] = fit(x, y, ft, 'StartPoint', start_point);
    end

    % Parameter CIs
    param_CI   = confint(fit_obj, 0.95);
    param_vals = coeffvalues(fit_obj);

    significance_check = (param_CI(1,:) < 0 & param_CI(2,:) < 0) | ...
                         (param_CI(1,:) > 0 & param_CI(2,:) > 0);

    % Fine x vector across plot range
    x_ = linspace(x_lim(1), x_lim(2), 200)';

    % 95% prediction interval
    y_predint = predint(fit_obj, x_, 0.95, 'functional', 'off');
    y_upper   = y_predint(:, 2);
    y_lower   = y_predint(:, 1);

    % Sealed reference uncertainty band (behind everything)
    yline_lo = cap_val - cap_unc;
    yline_hi = cap_val + cap_unc;
    x_fill   = [x_lim(1), x_lim(2), x_lim(2), x_lim(1)];
    y_fill   = [yline_lo, yline_lo, yline_hi, yline_hi];
    hcap = fill(x_fill, y_fill, [0.7 0.7 0.7]);
    hcap.FaceAlpha      = 0.30;
    hcap.EdgeColor      = 'none';
    hcap.HandleVisibility = 'off';

    % Sealed reference center line
    yline(cap_val, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.8, ...
        'Label', sprintf('Sealed: %.2f \\pm %.2f dB', cap_val, cap_unc), ...
        'LabelHorizontalAlignment', 'left', ...
        'FontSize', 12, 'FontName', 'Arial');

    % Fit prediction interval shading
    x2         = [x_', fliplr(x_')];
    in_between = [y_upper; flipud(y_lower)]';
    hf = fill(x2, in_between, line_color);
    hf.FaceAlpha = 0.15;
    hf.EdgeColor = 'none';
    hf.HandleVisibility = 'off';

    % Fit line — solid if significant, dashed if not
    if min(significance_check) > 0
        plot(x_, feval(fit_obj, x_), '-', 'Color', line_color, ...
            'LineWidth', 2.5, 'HandleVisibility', 'off');
    else
        plot(x_, feval(fit_obj, x_), '--', 'Color', line_color, ...
            'LineWidth', 2.5, 'HandleVisibility', 'off');
    end

    % Raw data points
    plot(x, y, 'o', 'MarkerFaceColor', line_color, ...
        'MarkerEdgeColor', line_color, 'MarkerSize', 5, ...
        'HandleVisibility', 'off');

    % Axis labels and title
    xlabel(x_label,  'FontSize', 18, 'FontName', 'Arial');
    ylabel(y_label,  'FontSize', 18, 'FontName', 'Arial');
    title(fig_title, 'FontSize', 20, 'FontName', 'Arial');

    ax = gca;
    ax.FontSize = 14; ax.FontName = 'Arial'; ax.LineWidth = 1.2;

    % Annotation box: equation as y = p1*x + p2 with CIs, no R^2
    p1     = param_vals(1);
    p2     = param_vals(2);
    p1_lo  = param_CI(1,1);  p1_hi  = param_CI(2,1);
    p2_lo  = param_CI(1,2);  p2_hi  = param_CI(2,2);

    annotation_str = sprintf(['y = p_{1} x + p_{2}\n' ...
        'p_{1} = %.4g  [%.4g, %.4g]\n' ...
        'p_{2} = %.4g  [%.4g, %.4g]'], ...
        p1, p1_lo, p1_hi, ...
        p2, p2_lo, p2_hi);

    text(0.03, 0.97, annotation_str, ...
        'Units',             'normalized', ...
        'VerticalAlignment', 'top', ...
        'FontSize',          11, ...
        'FontName',          'Arial', ...
        'BackgroundColor',   [1 1 1 0.88], ...
        'EdgeColor',         [0.4 0.4 0.4], ...
        'Margin',            6);
end

% FIT TYPE SETTINGS -------------------------------------------------------
fit_vol_len  = 'poly1';  sp_vol_len  = [1, 1];
fit_flat_len = 'poly1';  sp_flat_len = [1, 1];
fit_vol_dia  = 'poly1';  sp_vol_dia  = [1, 1];
fit_flat_dia = 'poly1';  sp_flat_dia = [1, 1];

col_len = [0.25 0.00 0.45];
col_dia = [0.55 0.00 0.00];

x_lim_len = [-10, 190];
x_lim_dia = [-5,  75];

% GENERATE THE FOUR LEVEL 3 PLOTS -----------------------------------------
plot_fit_with_CI(...
    x_len_all, y_vol_len_all, ...
    cap_vol_far, cap_vol_far_unc, ...
    'Port Length (mm)', 'Relative Output Level (dB)', ...
    'Relative Output Level vs. Port Length', ...
    fit_vol_len, sp_vol_len, col_len, x_lim_len);

plot_fit_with_CI(...
    x_len_all, y_flat_len_all, ...
    cap_flat_far, cap_flat_far_unc, ...
    'Port Length (mm)', 'RMS Deviation (dB)', ...
    'RMS Deviation vs. Port Length', ...
    fit_flat_len, sp_flat_len, col_len, x_lim_len);

plot_fit_with_CI(...
    x_dia_all, y_vol_dia_all, ...
    cap_vol_far, cap_vol_far_unc, ...
    'Port Diameter (mm)', 'Relative Output Level (dB)', ...
    'Relative Output Level vs. Port Diameter', ...
    fit_vol_dia, sp_vol_dia, col_dia, x_lim_dia);

plot_fit_with_CI(...
    x_dia_all, y_flat_dia_all, ...
    cap_flat_far, cap_flat_far_unc, ...
    'Port Diameter (mm)', 'RMS Deviation (dB)', ...
    'RMS Deviation vs. Port Diameter', ...
    fit_flat_dia, sp_flat_dia, col_dia, x_lim_dia);


% =========================================================================
% FLASH TALK SCATTER PLOT — Flatness vs Volume, all configurations
% Color = diameter, shape = length, sealed shown separately
% =========================================================================

figure('Position', [100 100 780 620]);
hold on; grid on;

% --- Color map for the 4 diameters (light → dark red) ---
dia_scatter_colors = [
    1.00 0.70 0.70;   % 17.5 mm — light red
    0.90 0.40 0.40;   % 35.0 mm
    0.75 0.15 0.15;   % 52.5 mm
    0.50 0.00 0.00;   % 70.0 mm — dark red
];

% --- Marker shapes for the 5 lengths ---
len_markers = {'o', 's', '^', 'v', 'd'};  % circle, square, up-tri, down-tri, diamond

marker_size = 10;

% --- Plot ported configurations ---
for di = 1:n_D
    for li = 1:n_L
        v = vol_far_mean(li, di);
        f = flat_far_mean(li, di);
        if isnan(v) || isnan(f)
            continue;
        end
        plot(v, f, ...
            len_markers{li}, ...
            'MarkerFaceColor', dia_scatter_colors(di,:), ...
            'MarkerEdgeColor', dia_scatter_colors(di,:) * 0.6, ...
            'MarkerSize', marker_size, ...
            'HandleVisibility', 'off');
    end
end

% --- Plot sealed reference ---
plot(cap_vol_far, cap_flat_far, ...
    'p', ...                          % pentagram for visual distinction
    'MarkerFaceColor', [0.4 0.4 0.4], ...
    'MarkerEdgeColor', [0.1 0.1 0.1], ...
    'MarkerSize', 14, ...
    'DisplayName', 'Sealed');

% --- Highlight best configuration ---
% Best = highest volume with flatness below median (adjust threshold as needed)
flat_thresh = median(flat_far_mean(:), 'omitnan');
vol_flat    = vol_far_mean;
vol_flat(flat_far_mean > flat_thresh) = NaN;
[~, best_idx] = max(vol_flat(:));
[best_li, best_di] = ind2sub([n_L, n_D], best_idx);

best_v = vol_far_mean(best_li, best_di);
best_f = flat_far_mean(best_li, best_di);

plot(best_v, best_f, 'k*', 'MarkerSize', 18, 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Best: L=%d mm, D=%.1f mm', ...
        lengths_mm(best_li), diameters_mm(best_di)));

% --- Shade "good region" (top-right: high volume, low flatness) ---
% Defined as above-median volume AND below-median flatness
vol_thresh = median(vol_far_mean(:), 'omitnan');
ax = gca;
xl = xlim; yl = ylim;  % get auto limits first
% nudge axis limits so patch doesn't get clipped
xlim([xl(1)-0.3, xl(2)+0.3]);
ylim([yl(1)-0.1, yl(2)+0.1]);
xl = xlim; yl = ylim;

patch([vol_thresh, xl(2), xl(2), vol_thresh], ...
      [yl(1), yl(1), flat_thresh, flat_thresh], ...
      [0.6 0.9 0.6], ...
      'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
      'HandleVisibility', 'off');

text(xl(2) - 0.05*(xl(2)-xl(1)), yl(1) + 0.06*(yl(2)-yl(1)), ...
    'loud + flat', ...
    'HorizontalAlignment', 'right', ...
    'FontSize', 11, 'FontName', 'Arial', ...
    'Color', [0.15 0.5 0.15]);

% --- Axis labels (flip y so "flatter = up") ---
set(gca, 'YDir', 'reverse');   % low RMS = top = better
xlabel('Relative Output Level (dB)', 'FontSize', 14, 'FontName', 'Arial');
ylabel('RMS Deviation (dB)  ←  flatter', 'FontSize', 14, 'FontName', 'Arial');
title('Port Geometry: Loudness–Flatness Tradeoff', 'FontSize', 16, 'FontName', 'Arial');

ax.FontSize = 12; ax.FontName = 'Arial'; ax.LineWidth = 1.2;

% --- Legend: diameter (color) + length (shape) ---
% Build manual legend handles
h_dia = gobjects(n_D, 1);
for di = 1:n_D
    h_dia(di) = plot(NaN, NaN, 'o', ...
        'MarkerFaceColor', dia_scatter_colors(di,:), ...
        'MarkerEdgeColor', dia_scatter_colors(di,:) * 0.6, ...
        'MarkerSize', marker_size, ...
        'DisplayName', sprintf('D = %.1f mm', diameters_mm(di)));
end

h_len = gobjects(n_L, 1);
for li = 1:n_L
    h_len(li) = plot(NaN, NaN, ...
        len_markers{li}, ...
        'MarkerFaceColor', [0.5 0.5 0.5], ...
        'MarkerEdgeColor', [0.3 0.3 0.3], ...
        'MarkerSize', marker_size, ...
        'DisplayName', sprintf('L = %d mm', lengths_mm(li)));
end

% Grab the sealed and best handles already plotted
h_sealed = findobj(gca, 'DisplayName', 'Sealed');
h_best   = findobj(gca, 'DisplayName', ...
    sprintf('Best: L=%d mm, D=%.1f mm', lengths_mm(best_li), diameters_mm(best_di)));

legend([h_dia; h_len; h_sealed(1); h_best(1)], ...
    'Location', 'eastoutside', 'FontSize', 10, 'FontName', 'Arial');