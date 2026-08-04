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

% VOLUME vs. LENGTH — Far Mic only
figure('Name','Volume vs. Length','Position',[200 200 700 500]);
hold on; grid on;
title('Average Volume vs. Port Length — Far Mic (1 m)');
xlabel('Port Length (mm)'); ylabel('Average Response (dB)');

for di = 1:n_D
    errorbar(lengths_mm, vol_far_mean(:,di)', vol_far_std(:,di)', ...
        '-o', 'Color', dia_colors(di,:), 'LineWidth', 1.5, ...
        'MarkerFaceColor', dia_colors(di,:), ...
        'DisplayName', dia_labels{di});
end
yline(cap_vol_far, '--k', 'LineWidth', 1.5, 'DisplayName', 'Capped', ...
    'Label', sprintf('Capped (%.1f dB)', cap_vol_far), 'LabelHorizontalAlignment','left');
legend('Location','best'); xlim([min(lengths_mm)-10, max(lengths_mm)+10]);

% FLATNESS vs. LENGTH — Far Mic only
figure('Name','Flatness vs. Length','Position',[200 200 700 500]);
hold on; grid on;
title('RMS Deviation vs. Port Length — Far Mic (1 m)');
xlabel('Port Length (mm)'); ylabel('RMS Deviation (dB)');

for di = 1:n_D
    errorbar(lengths_mm, flat_far_mean(:,di)', flat_far_std(:,di)', ...
        '-s', 'Color', dia_colors(di,:), 'LineWidth', 1.5, ...
        'MarkerFaceColor', dia_colors(di,:), ...
        'DisplayName', dia_labels{di});
end
yline(cap_flat_far, '--k', 'LineWidth', 1.5, ...
    'Label', sprintf('Capped (%.2f dB)', cap_flat_far), 'LabelHorizontalAlignment','left');
legend('Location','best'); xlim([min(lengths_mm)-10, max(lengths_mm)+10]);

% VOLUME vs. DIAMETER — Far Mic only
figure('Name','Volume vs. Diameter','Position',[200 200 700 500]);
hold on; grid on;
title('Average Volume vs. Port Diameter — Far Mic (1 m)');
xlabel('Port Diameter (mm)'); ylabel('Average Response (dB)');

for li = 1:n_L
    errorbar(diameters_mm, vol_far_mean(li,:), vol_far_std(li,:), ...
        '-o', 'Color', len_colors(li,:), 'LineWidth', 1.5, ...
        'MarkerFaceColor', len_colors(li,:), ...
        'DisplayName', len_labels{li});
end
yline(cap_vol_far, '--k', 'LineWidth', 1.5, ...
    'Label', sprintf('Capped (%.1f dB)', cap_vol_far), 'LabelHorizontalAlignment','left');
legend('Location','best'); xlim([min(diameters_mm)-5, max(diameters_mm)+5]);

% FLATNESS vs. DIAMETER — Far Mic only
figure('Name','Flatness vs. Diameter','Position',[200 200 700 500]);
hold on; grid on;
title('RMS Deviation vs. Port Diameter — Far Mic (1 m)');
xlabel('Port Diameter (mm)'); ylabel('RMS Deviation (dB)');

for li = 1:n_L
    errorbar(diameters_mm, flat_far_mean(li,:), flat_far_std(li,:), ...
        '-s', 'Color', len_colors(li,:), 'LineWidth', 1.5, ...
        'MarkerFaceColor', len_colors(li,:), ...
        'DisplayName', len_labels{li});
end
yline(cap_flat_far, '--k', 'LineWidth', 1.5, ...
    'Label', sprintf('Capped (%.2f dB)', cap_flat_far), 'LabelHorizontalAlignment','left');
legend('Location','best'); xlim([min(diameters_mm)-5, max(diameters_mm)+5]);