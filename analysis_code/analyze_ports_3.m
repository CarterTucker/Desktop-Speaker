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
% add directories to csv folder and input wav
csv_dir  = 'CSVs';
wav_file = 'sweep_sine.wav'; 

% set sampling rates for the mic and wav file
Fs_mic      = 10000;
wav_Fs_nominal = 48000;

% port geometries
lengths_mm = [0, 45, 90, 135, 180];
diameters_mm = [17.5, 35.0, 52.5, 70.0];

n_trials = 5;

% Analysis frequency range (set to be bass response, not total)
f_min = 126; 
f_max = 700;

% smoothing factor
smooth_frac = 0.05; 

% prevent division by 0 in deconvolution
reg_eps = 1e-6;

% LOAD THE WAV FILE -------------------------------------------------------
disp('Loading sweep .wav...\n');
[x_wav, Fs_wav] = audioread(wav_file);

% resample wav to mic sampling rate
resampled_wav = resample(x_wav, Fs_mic, Fs_wav);
resampled_wav = resampled_wav(:);

fprintf('Sweep loaded: %.1f s at %d Hz, resampled from %d Hz\n', ...
    length(x_wav)/Fs_wav, Fs_mic, Fs_wav);

% Define functions for deconvolution and grabbing parameters
    function [f, H_dB] = compute_response(y, resampled_wav, Fs, f_min, f_max, reg_eps, smooth_frac)

        % make sure the sizes match
        N = length(y);
        x = resampled_wav;
        if length(x) < N
            x = [x; zeros(N - length(x), 1)];
        else
            x = x(1:N);
        end

        output_FFT = fft(y, N);
        input_FFT = fft(x, N);

        % deconvolution
        H = output_FFT .* conj(input_FFT) ./ (abs(input_FFT).^2 + reg_eps * max(abs(input_FFT).^2));

        f_full = (0:N-1)' * Fs / N;
        H_full = abs(H(1:N));

        % keep the response in the specified range
        idx_ny = floor(N/2);
        f    = f_full(1:idx_ny);
        H_mag = H_full(1:idx_ny);

        % smooth noise with SMA
        win = max(1, round(smooth_frac * idx_ny));
        H_smooth = movmean(H_mag, win);

        % convert response to dB for Bode
        H_dB = 20 * log10(H_smooth + eps);

        % only look at the frequency range defined previously
        band = f >= f_min & f <= f_max;
        f    = f(band);
        H_dB = H_dB(band);
    end

    function [avg_vol, flatness] = compute_metrics(H_dB)
        % avg_vol  : mean of H_dB across frequency band (dB)
        % flatness : RMS deviation from mean (lower = flatter)
        avg_vol  = mean(H_dB);
        flatness = sqrt(mean((H_dB - avg_vol).^2));
    end

% PROCESS THE CSVs
fprintf('Processing CSV files...\n');

% find the number of lengths and diameters
n_L = length(lengths_mm);
n_D = length(diameters_mm);

% add NANs to geometries which have less than 5 trials
vol_close   = nan(n_L, n_D, n_trials);
vol_far    = nan(n_L, n_D, n_trials);
flat_close  = nan(n_L, n_D, n_trials);
flat_far   = nan(n_L, n_D, n_trials);

bode_store_close = cell(n_L, n_D);
bode_store_far  = cell(n_L, n_D);
bode_f          = [];

bode_f_store = cell(n_L, n_D);

% capped condition
vol_close_capped  = nan(1, n_trials);
vol_far_capped   = nan(1, n_trials);
flat_close_capped = nan(1, n_trials);
flat_far_capped  = nan(1, n_trials);

% do processing for non-capped port geometries
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

            % try to read the CSV
            try
                T = readtable(fpath, 'NumHeaderLines', 1, ...
                    'VariableNamingRule', 'preserve');
            catch
                warning('Could not read: %s', fpath);
                continue;
            end

            data = T{:,:}; 
            y_close = data(:,2);
            y_far  = data(:,4);

            % use function to get response over frequency
            [f_band, H_close] = compute_response(y_close, resampled_wav, Fs_mic, f_min, f_max, reg_eps, smooth_frac);
            [~, H_far]  = compute_response(y_far,  resampled_wav, Fs_mic, f_min, f_max, reg_eps, smooth_frac);

            if isempty(bode_f)
                bode_f = f_band;
            end

            % use function to get volume and flatness for each mic
            [vol_close(li,di,ti),  flat_close(li,di,ti)]  = compute_metrics(H_close);
            [vol_far(li,di,ti),   flat_far(li,di,ti)]   = compute_metrics(H_far);

            % store an example Bode for the progress report
            if ti == 1
                bode_store_close{li,di} = H_close;
                bode_store_far{li,di}  = H_far;
                bode_f_store{li,di}    = f_band;
            end

            fprintf('  Processed: %s\n', fname);
        end
    end
end

% do processing for the capped condition (to compare results with a
% sealed, non-ported design)
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

    data   = T{:,:};
    y_close = data(:,2);
    y_far  = data(:,4);

    % use function to get response over frequency
    [f_band, H_close] = compute_response(y_close, resampled_wav, Fs_mic, f_min, f_max, reg_eps, smooth_frac);
    [~, H_far]  = compute_response(y_far,  resampled_wav, Fs_mic, f_min, f_max, reg_eps, smooth_frac);

    if isempty(bode_f)
        bode_f = f_band; 
    end
    
    % use function to get volume and flatness for each mic
    [vol_close(li,di,ti),  flat_close(li,di,ti)]  = compute_metrics(H_close);
    [vol_far(li,di,ti),   flat_far(li,di,ti)]   = compute_metrics(H_far);

    disp("Processed");
end

% get the mean volume and flatness without NANs if any were found
vol_close_mean  = mean(vol_close,  3, 'omitnan');
vol_far_mean   = mean(vol_far,   3, 'omitnan');
flat_close_mean = mean(flat_close, 3, 'omitnan');
flat_far_mean  = mean(flat_far,  3, 'omitnan');

% get standard deviation volume and flatness without NANs again
vol_close_std   = std(vol_close,  0, 3, 'omitnan');
vol_far_std    = std(vol_far,   0, 3, 'omitnan');
flat_close_std  = std(flat_close, 0, 3, 'omitnan');
flat_far_std   = std(flat_far,  0, 3, 'omitnan');

% do the same for capped trials
cap_vol_close  = mean(vol_close_capped,  'omitnan');
cap_vol_far   = mean(vol_far_capped,   'omitnan');
cap_flat_close = mean(flat_close_capped, 'omitnan');
cap_flat_far  = mean(flat_far_capped,  'omitnan');

cap_vol_close_std  = std(vol_close_capped,  'omitnan');
cap_vol_far_std   = std(vol_far_capped,   'omitnan');
cap_flat_close_std = std(flat_close_capped, 'omitnan');
cap_flat_far_std  = std(flat_far_capped,  'omitnan');

% GRAPH EVERYTHING
% Red gradient for diameter (light to dark), purple gradient for length (light to dark)
red_shades = [1.00 0.70 0.70; 0.90 0.40 0.40; 0.75 0.15 0.15; 0.55 0.00 0.00];
purple_shades = [0.80 0.70 0.90; 0.65 0.45 0.80; 0.50 0.25 0.70; 0.38 0.10 0.58; 0.25 0.00 0.45];

dia_colors = red_shades;    % diameter lines = red gradient
len_colors = purple_shades; % length lines   = purple gradient

dia_labels = arrayfun(@(d) sprintf('%.1f mm', d), diameters_mm, 'UniformOutput', false);
len_labels = arrayfun(@(l) sprintf('%d mm', l), lengths_mm, 'UniformOutput', false);

% BODE PLOTS — Far Mic only, one figure per diameter, lines vary over length
for di = 1:n_D
    figure('Name', sprintf('Bode - Diameter %.1f mm', diameters_mm(di)), ...
        'Position', [100 100 700 500]);
    hold on; grid on;
    title(sprintf('Far Mic (1 m) — Dia = %.1f mm', diameters_mm(di)));
    xlabel('Frequency (Hz)'); ylabel('Response (dB)');
    set(gca, 'XScale', 'log');

    for li = 1:n_L
        if ~isempty(bode_store_far{li,di})
            plot(bode_f_store{li,di}, bode_store_far{li,di}, ...
                'Color', len_colors(li,:), ...
                'LineWidth', 1.5, 'DisplayName', len_labels{li});
        end
    end
    legend('Location','best');
    title(sprintf('Frequency Response — Port Diameter %.1f mm', diameters_mm(di)));
end


% =========================================================================
% THIELE-SMALL + ENCLOSURE PARAMETERS FOR MODEL PREDICTION
% =========================================================================
Fs_ts  = 53.8;        % Free-air resonance (Hz)
Qts_ts = 0.35;        % Total Q
Vas_ts = 3.1e-3;      % Equivalent compliance volume (m^3)
Vb     = 3.99e-3;     % Box volume (m^3)
c_spd  = 343;         % Speed of sound (m/s)
rho    = 1.21;        % Air density (kg/m^3)

% System parameters with box loading
Fs_sys  = Fs_ts * sqrt(1 + Vas_ts/Vb);
Qts_sys = Qts_ts * sqrt(1 + Vas_ts/Vb);
Cab     = Vb / (rho * c_spd^2);

% Fine frequency vector for smooth predicted curves
f_pred = linspace(f_min, f_max, 2000)';
w_pred = 2 * pi * f_pred;

% -------------------------------------------------------------------------
% HELPER: bass-reflex transfer function for one port geometry
% -------------------------------------------------------------------------
function H_dB = br_response(w, Fs_sys, Qts_sys, Cab, c_spd, rho, L_m, D_m)
    Sp = pi * (D_m/2)^2;
    if L_m == 0
        L_eff = 0.1 * D_m;
    else
        L_eff = L_m + 0.6 * (D_m/2);
    end
    fp = (c_spd/(2*pi)) * sqrt(Sp / (Cab * c_spd^2 * rho * L_eff));
    ws = 2*pi*Fs_sys;
    wp = 2*pi*fp;
    h  = w / ws;
    a1 = (1/Qts_sys) * (wp/ws);
    a2 = (wp/ws)^2 + 1 + 1/Qts_sys^2;
    a3 = (1/Qts_sys) * (1 + (wp/ws)^2);
    H_mag = abs(h.^4 ./ (h.^4 + a3.*h.^3 + a2.*h.^2 + a1.*h + 1));
    H_dB  = 20 * log10(H_mag + eps);
end

% -------------------------------------------------------------------------
% HELPER: compute volume and flatness
% -------------------------------------------------------------------------
function [avg_vol, flatness] = get_metrics(H_dB)
    avg_vol  = mean(H_dB);
    flatness = sqrt(mean((H_dB - avg_vol).^2));
end

% =========================================================================
% COLLAPSE DATA ACROSS OTHER VARIABLE
% =========================================================================
% vol_far and flat_far are (n_L x n_D x n_trials)
% Collapse all trials and all conditions of the OTHER variable into one
% vector per condition of the variable of interest

% --- vs. LENGTH: collapse across diameters and trials ---
% reshape to (n_L x n_D*n_trials), then mean/std across dim 2
vol_far_all  = reshape(vol_far,  n_L, []);   % (n_L x n_D*n_trials)
flat_far_all = reshape(flat_far, n_L, []);

vol_vs_len_mean  = mean(vol_far_all,  2, 'omitnan');  % (n_L x 1)
vol_vs_len_std   = std(vol_far_all,   0, 2, 'omitnan');
flat_vs_len_mean = mean(flat_far_all, 2, 'omitnan');
flat_vs_len_std  = std(flat_far_all,  0, 2, 'omitnan');

% --- vs. DIAMETER: collapse across lengths and trials ---
vol_far_perm  = permute(vol_far,  [2 1 3]);  % (n_D x n_L x n_trials)
flat_far_perm = permute(flat_far, [2 1 3]);
vol_far_Dall  = reshape(vol_far_perm,  n_D, []);
flat_far_Dall = reshape(flat_far_perm, n_D, []);

vol_vs_dia_mean  = mean(vol_far_Dall,  2, 'omitnan');  % (n_D x 1)
vol_vs_dia_std   = std(vol_far_Dall,   0, 2, 'omitnan');
flat_vs_dia_mean = mean(flat_far_Dall, 2, 'omitnan');
flat_vs_dia_std  = std(flat_far_Dall,  0, 2, 'omitnan');

% =========================================================================
% LINEAR REGRESSION HELPER
% =========================================================================
function [p, ci, y_fit, y_upper, y_lower] = linreg_with_ci(x, y, x_fine)
    % x, y: column vectors of data
    % Returns fit line, 95% CI band over x_fine
    n  = length(x);
    p  = polyfit(x, y, 1);
    y_hat = polyval(p, x);
    res   = y - y_hat;
    s2    = sum(res.^2) / (n - 2);
    Sxx   = sum((x - mean(x)).^2);
    t_crit = tinv(0.975, n-2);
    % CI band
    se_band = sqrt(s2 * (1/n + (x_fine - mean(x)).^2 / Sxx));
    y_fit   = polyval(p, x_fine);
    y_upper = y_fit + t_crit * se_band;
    y_lower = y_fit - t_crit * se_band;
    % Parameter CI
    se_slope = sqrt(s2 / Sxx);
    se_int   = sqrt(s2 * (1/n + mean(x)^2/Sxx));
    ci = [p(1) - t_crit*se_slope, p(1) + t_crit*se_slope;
          p(2) - t_crit*se_int,   p(2) + t_crit*se_int];
end

% =========================================================================
% PREDICTED MODEL LINES (collapsed across other variable)
% =========================================================================
% Sweep diameter continuously, average over all lengths
D_fine_mm = linspace(min(diameters_mm), max(diameters_mm), 200);
L_mean_m  = mean(lengths_mm(lengths_mm > 0)) / 1000;  % mean nonzero length

pred_vol_vs_dia  = zeros(size(D_fine_mm));
pred_flat_vs_dia = zeros(size(D_fine_mm));
for k = 1:length(D_fine_mm)
    H = br_response(w_pred, Fs_sys, Qts_sys, Cab, c_spd, rho, L_mean_m, D_fine_mm(k)/1000);
    [pred_vol_vs_dia(k), pred_flat_vs_dia(k)] = get_metrics(H);
end

% Sweep length continuously, average over all diameters
L_fine_mm = linspace(0, max(lengths_mm), 200);
D_mean_m  = mean(diameters_mm) / 1000;

pred_vol_vs_len  = zeros(size(L_fine_mm));
pred_flat_vs_len = zeros(size(L_fine_mm));
for k = 1:length(L_fine_mm)
    L_m = L_fine_mm(k) / 1000;
    if L_m == 0, L_m = 1e-4; end
    H = br_response(w_pred, Fs_sys, Qts_sys, Cab, c_spd, rho, L_m, D_mean_m);
    [pred_vol_vs_len(k), pred_flat_vs_len(k)] = get_metrics(H);
end

% Shift predicted curves to match data mean (shape comparison only)
pred_vol_vs_dia  = pred_vol_vs_dia  - mean(pred_vol_vs_dia)  + mean(vol_vs_dia_mean,  'omitnan');
pred_flat_vs_dia = pred_flat_vs_dia - mean(pred_flat_vs_dia) + mean(flat_vs_dia_mean, 'omitnan');
pred_vol_vs_len  = pred_vol_vs_len  - mean(pred_vol_vs_len)  + mean(vol_vs_len_mean,  'omitnan');
pred_flat_vs_len = pred_flat_vs_len - mean(pred_flat_vs_len) + mean(flat_vs_len_mean, 'omitnan');

% =========================================================================
% ANNOTATION HELPER
% =========================================================================
function annotation_box(ax, p, ci, fp_mean)
    xl = xlim(ax); yl = ylim(ax);
    tx = xl(1) + 0.03*(xl(2)-xl(1));
    ty = yl(2) - 0.03*(yl(2)-yl(1));
    str = sprintf(['Linear fit\n' ...
        'Slope: %.4f  [%.4f, %.4f]\n' ...
        'Intercept: %.3f  [%.3f, %.3f]\n' ...
        'Model f_p (mean): %.1f Hz'], ...
        p(1), ci(1,1), ci(1,2), ...
        p(2), ci(2,1), ci(2,2), ...
        fp_mean);
    text(ax, tx, ty, str, ...
        'VerticalAlignment',   'top', ...
        'HorizontalAlignment', 'left', ...
        'FontSize', 10, ...
        'FontName', 'Arial', ...
        'BackgroundColor', [1 1 1 0.8], ...
        'EdgeColor', [0.4 0.4 0.4], ...
        'Margin', 5);
end

% =========================================================================
% PLOT COLORS (single color per graph)
% =========================================================================
col_dia = [0.55 0.00 0.00];   % dark red  — diameter graphs
col_len = [0.25 0.00 0.45];   % dark purple — length graphs

% =========================================================================
% FINE X VECTORS FOR REGRESSION CI BAND
% =========================================================================
x_dia_fine = linspace(min(diameters_mm), max(diameters_mm), 200)';
x_len_fine = linspace(min(lengths_mm),   max(lengths_mm),   200)';

% =========================================================================
% MEAN HELMHOLTZ FREQUENCY FOR ANNOTATION
% =========================================================================
fp_at_dia = zeros(size(diameters_mm));
for di = 1:n_D
    Sp = pi*(diameters_mm(di)/2000)^2;
    L_eff = L_mean_m + 0.6*(diameters_mm(di)/2000);
    fp_at_dia(di) = (c_spd/(2*pi))*sqrt(Sp/(Cab*c_spd^2*rho*L_eff));
end
fp_mean_dia = mean(fp_at_dia);

fp_at_len = zeros(size(lengths_mm(lengths_mm>0)));
Lnn = lengths_mm(lengths_mm>0);
for li = 1:length(Lnn)
    Sp = pi*(D_mean_m/2)^2;
    L_eff = Lnn(li)/1000 + 0.6*(D_mean_m/2);
    fp_at_len(li) = (c_spd/(2*pi))*sqrt(Sp/(Cab*c_spd^2*rho*L_eff));
end
fp_mean_len = mean(fp_at_len);

% =========================================================================
% FIGURE 1: VOLUME vs. LENGTH
% =========================================================================
[p_vl, ci_vl, yfit_vl, yup_vl, ylo_vl] = linreg_with_ci(...
    lengths_mm', vol_vs_len_mean, x_len_fine);

figure('Name','Volume vs Length (Level 3)','Position',[100 100 700 520]);
ax1 = axes; hold on; grid on;
% Predicted model curve
plot(L_fine_mm, pred_vol_vs_len, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, ...
    'DisplayName', 'Lumped-parameter model');
% CI shading
fill([x_len_fine; flipud(x_len_fine)], [yup_vl; flipud(ylo_vl)], ...
    col_len, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
% Linear fit
plot(x_len_fine, yfit_vl, '-', 'Color', col_len, 'LineWidth', 2, ...
    'DisplayName', 'Linear fit');
% Data points
errorbar(lengths_mm, vol_vs_len_mean, vol_vs_len_std, ...
    'o', 'Color', col_len, 'MarkerFaceColor', col_len, ...
    'LineWidth', 1.5, 'MarkerSize', 7, 'DisplayName', 'Measured (mean \pm std)');
xlabel('Port Length (mm)', 'FontSize', 14, 'FontName', 'Arial');
ylabel('Average Response (dB)', 'FontSize', 14, 'FontName', 'Arial');
title('Volume vs. Port Length — Far Mic (1 m)', 'FontSize', 15, 'FontName', 'Arial');
legend('Location', 'northeast', 'FontSize', 11);
ax1.FontSize = 12; ax1.FontName = 'Arial'; ax1.LineWidth = 1.1;
xlim([min(lengths_mm)-10, max(lengths_mm)+10]);
annotation_box(ax1, p_vl, ci_vl, fp_mean_len);

% =========================================================================
% FIGURE 2: FLATNESS vs. LENGTH
% =========================================================================
[p_fl, ci_fl, yfit_fl, yup_fl, ylo_fl] = linreg_with_ci(...
    lengths_mm', flat_vs_len_mean, x_len_fine);

figure('Name','Flatness vs Length (Level 3)','Position',[100 100 700 520]);
ax2 = axes; hold on; grid on;
plot(L_fine_mm, pred_flat_vs_len, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, ...
    'DisplayName', 'Lumped-parameter model');
fill([x_len_fine; flipud(x_len_fine)], [yup_fl; flipud(ylo_fl)], ...
    col_len, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(x_len_fine, yfit_fl, '-', 'Color', col_len, 'LineWidth', 2, ...
    'DisplayName', 'Linear fit');
errorbar(lengths_mm, flat_vs_len_mean, flat_vs_len_std, ...
    'o', 'Color', col_len, 'MarkerFaceColor', col_len, ...
    'LineWidth', 1.5, 'MarkerSize', 7, 'DisplayName', 'Measured (mean \pm std)');
xlabel('Port Length (mm)', 'FontSize', 14, 'FontName', 'Arial');
ylabel('RMS Deviation (dB)', 'FontSize', 14, 'FontName', 'Arial');
title('RMS Deviation vs. Port Length — Far Mic (1 m)', 'FontSize', 15, 'FontName', 'Arial');
legend('Location', 'northeast', 'FontSize', 11);
ax2.FontSize = 12; ax2.FontName = 'Arial'; ax2.LineWidth = 1.1;
xlim([min(lengths_mm)-10, max(lengths_mm)+10]);
annotation_box(ax2, p_fl, ci_fl, fp_mean_len);

% =========================================================================
% FIGURE 3: VOLUME vs. DIAMETER
% =========================================================================
[p_vd, ci_vd, yfit_vd, yup_vd, ylo_vd] = linreg_with_ci(...
    diameters_mm', vol_vs_dia_mean, x_dia_fine);

figure('Name','Volume vs Diameter (Level 3)','Position',[100 100 700 520]);
ax3 = axes; hold on; grid on;
plot(D_fine_mm, pred_vol_vs_dia, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, ...
    'DisplayName', 'Lumped-parameter model');
fill([x_dia_fine; flipud(x_dia_fine)], [yup_vd; flipud(ylo_vd)], ...
    col_dia, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(x_dia_fine, yfit_vd, '-', 'Color', col_dia, 'LineWidth', 2, ...
    'DisplayName', 'Linear fit');
errorbar(diameters_mm, vol_vs_dia_mean, vol_vs_dia_std, ...
    'o', 'Color', col_dia, 'MarkerFaceColor', col_dia, ...
    'LineWidth', 1.5, 'MarkerSize', 7, 'DisplayName', 'Measured (mean \pm std)');
xlabel('Port Diameter (mm)', 'FontSize', 14, 'FontName', 'Arial');
ylabel('Average Response (dB)', 'FontSize', 14, 'FontName', 'Arial');
title('Volume vs. Port Diameter — Far Mic (1 m)', 'FontSize', 15, 'FontName', 'Arial');
legend('Location', 'northwest', 'FontSize', 11);
ax3.FontSize = 12; ax3.FontName = 'Arial'; ax3.LineWidth = 1.1;
xlim([min(diameters_mm)-5, max(diameters_mm)+5]);
annotation_box(ax3, p_vd, ci_vd, fp_mean_dia);

% =========================================================================
% FIGURE 4: FLATNESS vs. DIAMETER
% =========================================================================
[p_fd, ci_fd, yfit_fd, yup_fd, ylo_fd] = linreg_with_ci(...
    diameters_mm', flat_vs_dia_mean, x_dia_fine);

figure('Name','Flatness vs Diameter (Level 3)','Position',[100 100 700 520]);
ax4 = axes; hold on; grid on;
plot(D_fine_mm, pred_flat_vs_dia, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, ...
    'DisplayName', 'Lumped-parameter model');
fill([x_dia_fine; flipud(x_dia_fine)], [yup_fd; flipud(ylo_fd)], ...
    col_dia, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(x_dia_fine, yfit_fd, '-', 'Color', col_dia, 'LineWidth', 2, ...
    'DisplayName', 'Linear fit');
errorbar(diameters_mm, flat_vs_dia_mean, flat_vs_dia_std, ...
    'o', 'Color', col_dia, 'MarkerFaceColor', col_dia, ...
    'LineWidth', 1.5, 'MarkerSize', 7, 'DisplayName', 'Measured (mean \pm std)');
xlabel('Port Diameter (mm)', 'FontSize', 14, 'FontName', 'Arial');
ylabel('RMS Deviation (dB)', 'FontSize', 14, 'FontName', 'Arial');
title('RMS Deviation vs. Port Diameter — Far Mic (1 m)', 'FontSize', 15, 'FontName', 'Arial');
legend('Location', 'northeast', 'FontSize', 11);
ax4.FontSize = 12; ax4.FontName = 'Arial'; ax4.LineWidth = 1.1;
xlim([min(diameters_mm)-5, max(diameters_mm)+5]);
annotation_box(ax4, p_fd, ci_fd, fp_mean_dia);
