clear; close all

% Pick file
[f,p] = uigetfile({'*.smr;*.smrx','Spike2 files (*.smr, *.smrx)'}, 'Select spike2 file to analyze');
if isequal(f,0), return; end
smrFile = fullfile(p,f);

% Point to your CEDS64ML folder
cedPath = 'C:\CEDMATLAB\CEDS64ML';

% Extract DigMark ch32
out = extractDigMark32_CEDS64(smrFile, cedPath);

% Decode to trials
decoded = decodeEyeblinkDigMark(out);
trialT = struct2table(decoded.trials);
% Neuropixels start and stop times
% NPX recording segments (start / end)
npx_start = [out.time_sec(2); out.time_sec(4)];
npx_end   = [out.time_sec(3); out.time_sec(end)];

spike2npxT = table(npx_start, npx_end, ...
    'VariableNames', {'npx_start_s','npx_end_s'});
spike2npxT.recording = ["spontaneous"; "task"];

% Spike2 file name without extension
[~, baseName, ~] = fileparts(smrFile);

% Output Excel file (same folder)
excelFile = fullfile(p, [baseName '_eyeblink_timestamps.xlsx']);

writetable(trialT,      excelFile, 'Sheet', 'Trials');
writetable(spike2npxT,  excelFile, 'Sheet', 'Neuropixels');

% Subtract spike2npxT start time (row 1) from all timestamp columns in trialT
% to align trial timestamps to the NPX recording clock
trialT_npx = trialT;
timestampCols = trialT.Properties.VariableNames(1:8);
for c = 1:8
    trialT_npx.(timestampCols{c}) = trialT.(timestampCols{c}) - spike2npxT.npx_start_s(1);
end
writetable(trialT_npx, excelFile, 'Sheet', 'Trials_NPX');

%% Extract Neuropixels TTL edge times (GLX clock)
npxTimes = npxTimestamps;

%% Match Spike2 and Neuropixels time

% 1. Calculate Durations
% t1_s2, t2_s2: Times of start/end pulses in Spike2 (seconds).
% t1_npx, t2_npx: Times of start/end pulses in NPX (seconds).
dur_s2(1)  = npx_end(1) - npx_start(1); %t2_s2 - t1_s2;
dur_s2(2)  = npx_end(2) - npx_start(2);

dur_npx(1) = npxTimes.fourEdgeTimesS(2) - npxTimes.fourEdgeTimesS(1); %t2_npx - t1_npx
dur_npx(2) = npxTimes.fourEdgeTimesS(4) - npxTimes.fourEdgeTimesS(3);

%% Import unit spike timestamps from workspace

% 1. Select the workspace file
[file, pth] = uigetfile('*.mat', 'Select spike timestamp workspace (.mat)');
if isequal(file, 0), return; end
filePath = fullfile(pth, file);

% 2. Read sampling rate from ap.meta file
[metaFile, metaPath] = uigetfile('*ap.meta', 'Select ap.meta file');
if isequal(metaFile, 0), return; end
metaFilePath = fullfile(metaPath, metaFile);

fid = fopen(metaFilePath, 'r');
metaText = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);
metaLines = metaText{1};
fs_npx = NaN;
for k = 1:numel(metaLines)
    if startsWith(metaLines{k}, 'imSampRate')
        parts = strsplit(metaLines{k}, '=');
        fs_npx = str2double(strtrim(parts{2}));
        break;
    end
end
if isnan(fs_npx)
    error('Could not find imSampRate in the selected .meta file.');
end
fprintf('Sampling rate read from ap.meta: %.1f Hz\n', fs_npx);

% 3. Load workspace and collect all variables into a cell array
% Each variable is expected to be a vector of sample indices for one unit
wsVars = load(filePath);
varNames = fieldnames(wsVars);
numUnits = numel(varNames);
spikeSamples = cell(numUnits, 1);
for i = 1:numUnits
    spikeSamples{i} = double(wsVars.(varNames{i}));
end
fprintf('Loaded %d units from workspace.\n', numUnits);

% 4. Segment spike timestamps using edge times converted to AP sample indices
% fourEdgeTimesS is clock-agnostic; multiply by fs_npx (AP clock) to get AP samples
% fs_nidq pulled directly from npxTimes rather than requiring a separate meta file read
fs_nidq = npxTimes.fsHz;
spont_start_samp = round(npxTimes.fourEdgeTimesS(1) * fs_npx);
spont_end_samp   = round(npxTimes.fourEdgeTimesS(2) * fs_npx);
task_start_samp  = round(npxTimes.fourEdgeTimesS(3) * fs_npx);
task_end_samp    = round(npxTimes.fourEdgeTimesS(4) * fs_npx);

spikeSamples_spont = cell(numUnits, 1);
spikeSamples_task  = cell(numUnits, 1);

for i = 1:numUnits
    ts = spikeSamples{i};
    spikeSamples_spont{i} = ts(ts >= spont_start_samp & ts <= spont_end_samp);
    spikeSamples_task{i}  = ts(ts >= task_start_samp  & ts <= task_end_samp);
end
fprintf('Segmentation complete. Units: %d | Spont spikes (unit 1): %d | Task spikes (unit 1): %d\n', ...
    numUnits, numel(spikeSamples_spont{1}), numel(spikeSamples_task{1}));

% 5. Subtract spont_start_samp from all timestamps, then convert to seconds
spikeTimestamps_spont = cell(numUnits, 1);
spikeTimestamps_task  = cell(numUnits, 1);

for i = 1:numUnits
    spikeTimestamps_spont{i} = (spikeSamples_spont{i} - spont_start_samp) / fs_npx;
    spikeTimestamps_task{i}  = (spikeSamples_task{i}  - spont_start_samp) / fs_npx;
end
fprintf('Timestamps zero-referenced to fourEdgeTimesS(1) and converted to seconds.\n');

%% Group task spike timestamps around trial events
nTrials = height(trialT_npx);

baseAPs    = cell(nTrials, numUnits);
toneAPs    = cell(nTrials, numUnits);
traceAPs   = cell(nTrials, numUnits);
airPuffAPs = cell(nTrials, numUnits);

for u = 1:numUnits
    ts = spikeTimestamps_task{u};
    for t = 1:nTrials
        camera_on  = trialT_npx.camera_on(t);
        tone_on    = trialT_npx.tone_on(t);
        tone_off   = trialT_npx.tone_off(t);
        trace_start = trialT_npx.trace_start(t);
        airpuff_on  = trialT_npx.airpuff_on(t);
        airpuff_off = trialT_npx.airpuff_off(t);

        baseAPs{t,u}    = ts(ts >= camera_on   & ts < tone_on);
        toneAPs{t,u}    = ts(ts >= tone_on     & ts < tone_off);
        traceAPs{t,u}   = ts(ts >= trace_start & ts < airpuff_on);
        airPuffAPs{t,u} = ts(ts >= airpuff_on  & ts <= airpuff_off + 0.1);
    end
end
fprintf('Event-grouped spike times computed for %d trials x %d units.\n', nTrials, numUnits);

%% Raster plots — one figure per unit, full trial aligned to camera_on = 0
% for u = 1:numUnits
%     figure('Name', sprintf('Unit %d Raster', u), 'NumberTitle', 'off');
%     hold on;
% 
%     for t = 1:nTrials
%         t0        = trialT_npx.camera_on(t);
%         trial_end = trialT_npx.airpuff_off(t) + 0.1;
%         ts        = spikeTimestamps_task{u};
%         spk       = ts(ts >= t0 & ts <= trial_end) - t0;
%         plot(spk, t * ones(size(spk)), '|', 'Color', 'k', 'MarkerSize', 6);
%     end
% 
%     % Shade epoch regions using mean event times relative to camera_on
%     mean_t0       = mean(trialT_npx.camera_on);
%     tone_start    = mean(trialT_npx.tone_on)     - mean_t0;
%     tone_end      = mean(trialT_npx.tone_off)    - mean_t0;
%     trace_start   = mean(trialT_npx.trace_start) - mean_t0;
%     airpuff_start = mean(trialT_npx.airpuff_on)  - mean_t0;
%     airpuff_end   = mean(trialT_npx.airpuff_off) - mean_t0 + 0.1;
% 
%     ylims = [0.5, nTrials + 0.5];
%     patch([tone_start    tone_end      tone_end      tone_start],    [ylims(1) ylims(1) ylims(2) ylims(2)], [0.2 0.6 0.2], 'FaceAlpha', 0.08, 'EdgeColor', 'none');
%     patch([trace_start   airpuff_start airpuff_start trace_start],   [ylims(1) ylims(1) ylims(2) ylims(2)], [0.2 0.2 0.8], 'FaceAlpha', 0.08, 'EdgeColor', 'none');
%     patch([airpuff_start airpuff_end   airpuff_end   airpuff_start], [ylims(1) ylims(1) ylims(2) ylims(2)], [0.8 0.2 0.2], 'FaceAlpha', 0.08, 'EdgeColor', 'none');
% 
%     xlabel('Time relative to camera on (s)');
%     ylabel('Trial');
%     title(sprintf('Unit %d', u));
%     ylim(ylims);
%     xlim([0, airpuff_end]);
%     set(gca, 'YDir', 'reverse');
%     hold off;
% end


%% Common parameters and epoch definitions
binWidth = 0.01;  % 5ms bins

mean_cam     = mean(trialT_npx.camera_on);
mean_toneon  = mean(trialT_npx.tone_on)     - mean_cam;
mean_toneoff = mean(trialT_npx.tone_off)    - mean_cam;
mean_trace   = mean(trialT_npx.trace_start) - mean_cam;
mean_puffon  = mean(trialT_npx.airpuff_on)  - mean_cam;
mean_puffoff = mean(trialT_npx.airpuff_off) - mean_cam + 0.1;

bin_edges   = 0 : binWidth : mean_puffoff;
bin_centers = bin_edges(1:end-1) + binWidth / 2;
nBins       = numel(bin_centers);

epochs(1).name = 'Baseline'; epochs(1).trigger = 'camera_on';   epochs(1).t0 = 0;           epochs(1).win_start = 0;                  epochs(1).win_end = mean_toneon;  epochs(1).color = [0.4 0.4 0.4];
epochs(2).name = 'Tone';     epochs(2).trigger = 'tone_on';     epochs(2).t0 = mean_toneon; epochs(2).win_start = mean_toneon - 0.05; epochs(2).win_end = mean_toneoff; epochs(2).color = [0.2 0.6 0.2];
epochs(3).name = 'Trace';    epochs(3).trigger = 'trace_start'; epochs(3).t0 = mean_trace;  epochs(3).win_start = mean_trace - 0.05;  epochs(3).win_end = mean_puffon;  epochs(3).color = [0.2 0.2 0.8];
epochs(4).name = 'Airpuff';  epochs(4).trigger = 'airpuff_on';  epochs(4).t0 = mean_puffon; epochs(4).win_start = mean_puffon - 0.05; epochs(4).win_end = mean_puffoff; epochs(4).color = [0.8 0.2 0.2];

%% Compute PSTH for all units
smoothFR = zeros(numUnits, nBins);
for u = 1:numUnits
    counts = zeros(1, nBins);
    for t = 1:nTrials
        cam_t = trialT_npx.camera_on(t);
        ts    = spikeTimestamps_task{u};
        spk   = ts(ts >= cam_t & ts <= cam_t + mean_puffoff) - cam_t;
        counts = counts + histcounts(spk, bin_edges);
    end
    smoothFR(u,:) = counts / nTrials / binWidth;
end

% %% Peristimulus histograms — one figure per unit, shared x-axis (0 = camera_on)
% for u = 1:numUnits
%     figure('Name', sprintf('Unit %d PSTH', u), 'NumberTitle', 'off');
%     for e = 1:4
%         subplot(4, 1, e);
%         hold on;
%         mask = bin_centers >= epochs(e).win_start & bin_centers <= epochs(e).win_end;
%         bar(bin_centers(mask), smoothFR(u, mask), 1, 'FaceColor', epochs(e).color, 'EdgeColor', 'none');
%         xline(epochs(e).t0, '-k', epochs(e).name, 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
%         ylabel('Firing rate (Hz)');
%         xlim([0, mean_puffoff]);
%         if e == 4, xlabel('Time relative to camera on (s)'); end
%         hold off;
%     end
%     all_ax = findobj(gcf, 'Type', 'Axes');
%     max_ylim = max(arrayfun(@(ax) ax.YLim(2), all_ax));
%     set(all_ax, 'YLim', [0, max_ylim]);
%     sgtitle(sprintf('Unit %d PSTH', u));
% end

%% Raster + PSTH per unit, 4 epoch panels
% for u = 1:numUnits
%     figure('Name', sprintf('Unit %d Raster+PSTH', u), 'NumberTitle', 'off');
%     for e = 1:4
%         % --- Raster ---
%         subplot(4, 2, (e-1)*2 + 1);
%         hold on;
%         for t = 1:nTrials
%             cam_t  = trialT_npx.camera_on(t);
%             trig_t = trialT_npx.(epochs(e).trigger)(t);
%             ts     = spikeTimestamps_task{u};
%             spk    = ts(ts >= trig_t + (epochs(e).win_start - epochs(e).t0) & ...
%                         ts <  trig_t + (epochs(e).win_end   - epochs(e).t0)) - cam_t;
%             plot(spk, t * ones(size(spk)), '|', 'Color', epochs(e).color, 'MarkerSize', 4);
%         end
%         xline(epochs(e).t0, '-k', 'LineWidth', 1);
%         xlim([epochs(e).win_start, epochs(e).win_end]);
%         ylim([0.5, nTrials + 0.5]);
%         set(gca, 'YDir', 'reverse');
%         ylabel('Trial');
%         title(epochs(e).name);
%         if e == 4, xlabel('Time rel. camera on (s)'); end
%         hold off;
% 
%         % --- PSTH ---
%         subplot(4, 2, (e-1)*2 + 2);
%         hold on;
%         mask = bin_centers >= epochs(e).win_start & bin_centers <= epochs(e).win_end;
%         bar(bin_centers(mask), smoothFR(u, mask), 1, 'FaceColor', epochs(e).color, 'EdgeColor', 'none');
%         xline(epochs(e).t0, '-k', 'LineWidth', 1);
%         xlim([epochs(e).win_start, epochs(e).win_end]);
%         ylabel('Firing rate (Hz)');
%         if e == 4, xlabel('Time rel. camera on (s)'); end
%         hold off;
%     end
%     sgtitle(sprintf('Unit %d', u));
% end

%% Full trial raster + PSTH per unit
bin_edges_full   = 0 : binWidth : mean_puffoff;
bin_centers_full = bin_edges_full(1:end-1) + binWidth / 2;

for u = 1:numUnits
    figure('Name', sprintf('Unit %d Full Trial', u), 'NumberTitle', 'off');

    subplot(2, 1, 1);
    hold on;
    for t = 1:nTrials
        cam_t = trialT_npx.camera_on(t);
        ts    = spikeTimestamps_task{u};
        spk   = ts(ts >= cam_t & ts <= cam_t + mean_puffoff) - cam_t;
        plot(spk, t * ones(size(spk)), '|', 'Color', 'k', 'MarkerSize', 4);
    end
    ylims = [0.5, nTrials + 0.5];
    patch([mean_toneon  mean_toneoff mean_toneoff mean_toneon],  [ylims(1) ylims(1) ylims(2) ylims(2)], [0.2 0.6 0.2], 'FaceAlpha', 0.08, 'EdgeColor', 'none');
    patch([mean_trace   mean_puffon  mean_puffon  mean_trace],   [ylims(1) ylims(1) ylims(2) ylims(2)], [0.2 0.2 0.8], 'FaceAlpha', 0.08, 'EdgeColor', 'none');
    patch([mean_puffon  mean_puffoff mean_puffoff mean_puffon],  [ylims(1) ylims(1) ylims(2) ylims(2)], [0.8 0.2 0.2], 'FaceAlpha', 0.08, 'EdgeColor', 'none');
    set(gca, 'YDir', 'reverse');
    xlim([0, mean_puffoff]); ylim(ylims);
    ylabel('Trial'); title(sprintf('Unit %d — Full Trial', u));
    hold off;

    subplot(2, 1, 2);
    hold on;
    counts_full = zeros(1, numel(bin_centers_full));
    for t = 1:nTrials
        cam_t = trialT_npx.camera_on(t);
        ts    = spikeTimestamps_task{u};
        spk   = ts(ts >= cam_t & ts <= cam_t + mean_puffoff) - cam_t;
        counts_full = counts_full + histcounts(spk, bin_edges_full);
    end
    fr_full = counts_full / nTrials / binWidth;
    bar(bin_centers_full, fr_full, 1, 'FaceColor', 'k', 'EdgeColor', 'none');
    yl = ylim;
    patch([mean_toneon  mean_toneoff mean_toneoff mean_toneon],  [yl(1) yl(1) yl(2) yl(2)], [0.2 0.6 0.2], 'FaceAlpha', 0.08, 'EdgeColor', 'none');
    patch([mean_trace   mean_puffon  mean_puffon  mean_trace],   [yl(1) yl(1) yl(2) yl(2)], [0.2 0.2 0.8], 'FaceAlpha', 0.08, 'EdgeColor', 'none');
    patch([mean_puffon  mean_puffoff mean_puffoff mean_puffon],  [yl(1) yl(1) yl(2) yl(2)], [0.8 0.2 0.2], 'FaceAlpha', 0.08, 'EdgeColor', 'none');
    xlim([0, mean_puffoff]);
    xlabel('Time relative to camera on (s)'); ylabel('Firing rate (Hz)');
    hold off;
end

%% Response characterization — classify units and compute response latency
alpha = 0.05;
unitClass   = struct();
respLatency = nan(numUnits, 4);
epoch_labels = {'Baseline','Tone','Trace','Airpuff'};

% Compute per-trial firing rates for each epoch
trialRates = zeros(nTrials, 4, numUnits);
for u = 1:numUnits
    for t = 1:nTrials
        cam_t = trialT_npx.camera_on(t);
        ts    = spikeTimestamps_task{u};
        for e = 1:4
            trig_t  = trialT_npx.(epochs(e).trigger)(t);
            win_dur = epochs(e).win_end - epochs(e).t0;  % from trigger onset to epoch end
            spk     = ts(ts >= trig_t & ts < trig_t + win_dur);
            trialRates(t,e,u) = numel(spk) / win_dur;
        end
    end
end

for u = 1:numUnits
    baseline_rates = trialRates(:,1,u);
    unitClass(u).unit = u;

    % Test normality of baseline with Lilliefors test
    [h_norm, ~] = lillietest(baseline_rates);

    for e = 2:4
        epoch_rates = trialRates(:,e,u);

        if h_norm == 0
            % Normal — use one-tailed t-tests
            [~, p_exc] = ttest2(baseline_rates, epoch_rates, 'Tail', 'left');   % epoch > baseline
            [~, p_inh] = ttest2(baseline_rates, epoch_rates, 'Tail', 'right');  % epoch < baseline
            test_used = 't-test';
        else
            % Non-normal — use one-tailed Wilcoxon rank-sum
            p_exc = ranksum(baseline_rates, epoch_rates, 'Tail', 'left');   % epoch > baseline
            p_inh = ranksum(baseline_rates, epoch_rates, 'Tail', 'right');  % epoch < baseline
            test_used = 'Wilcoxon';
        end

        if p_exc < alpha
            unitClass(u).(epoch_labels{e}) = 'excited';
            unitClass(u).([epoch_labels{e} '_p']) = p_exc;
        elseif p_inh < alpha
            unitClass(u).(epoch_labels{e}) = 'inhibited';
            unitClass(u).([epoch_labels{e} '_p']) = p_inh;
        else
            unitClass(u).(epoch_labels{e}) = 'unresponsive';
            unitClass(u).([epoch_labels{e} '_p']) = min(p_exc, p_inh);
        end
        unitClass(u).([epoch_labels{e} '_test']) = test_used;
    end

    % Response latency — first bin exceeding mean+3SD of baseline PSTH
    base_edges = epochs(1).win_start : binWidth : epochs(1).win_end;
    base_cents = base_edges(1:end-1) + binWidth/2;
    base_counts = zeros(1, numel(base_cents));
    for t = 1:nTrials
        cam_t = trialT_npx.camera_on(t);
        ts    = spikeTimestamps_task{u};
        spk   = ts(ts >= cam_t + epochs(1).win_start & ts < cam_t + epochs(1).win_end) - cam_t;
        base_counts = base_counts + histcounts(spk, base_edges);
    end
    base_fr   = base_counts / nTrials / binWidth;
    threshold = mean(base_fr) + 3 * std(base_fr);

    for e = 2:4
        edges = epochs(e).win_start : binWidth : epochs(e).win_end;
        cents = edges(1:end-1) + binWidth/2;
        counts = zeros(1, numel(cents));
        for t = 1:nTrials
            cam_t  = trialT_npx.camera_on(t);
            trig_t = trialT_npx.(epochs(e).trigger)(t);
            offset = trig_t - cam_t - epochs(e).t0;
            ts     = spikeTimestamps_task{u};
            spk    = ts(ts >= cam_t + epochs(e).win_start + offset & ...
                        ts <  cam_t + epochs(e).win_end   + offset) - cam_t - offset;
            counts = counts + histcounts(spk, edges);
        end
        fr = counts / nTrials / binWidth;

        upper_threshold = mean(base_fr) + 3 * std(base_fr);
        lower_threshold = mean(base_fr) - 3 * std(base_fr);

        % Only search for latency in bins at or after epoch onset
        epoch_mask    = cents >= epochs(e).t0;
        fr_masked     = fr;
        fr_masked(~epoch_mask) = NaN;

        if strcmp(unitClass(u).(epoch_labels{e}), 'excited')
            latency_idx = find(fr_masked > upper_threshold, 1, 'first');
        elseif strcmp(unitClass(u).(epoch_labels{e}), 'inhibited')
            latency_idx = find(fr_masked < lower_threshold, 1, 'first');
        else
            latency_idx = [];
        end

        if ~isempty(latency_idx)
            respLatency(u,e) = cents(latency_idx) - epochs(e).t0;
        end
    end
end

% Print summary
fprintf('\n--- Unit Classification ---\n');
for u = 1:numUnits
    fprintf('Unit %d | Tone: %-12s (p=%.3f, %s) | Trace: %-12s (p=%.3f, %s) | Airpuff: %-12s (p=%.3f, %s)\n', ...
        u, ...
        unitClass(u).Tone,    unitClass(u).Tone_p,    unitClass(u).Tone_test, ...
        unitClass(u).Trace,   unitClass(u).Trace_p,   unitClass(u).Trace_test, ...
        unitClass(u).Airpuff, unitClass(u).Airpuff_p, unitClass(u).Airpuff_test);
end
fprintf('\n--- Response Latencies (s relative to epoch onset) ---\n');
fprintf('%-8s %-12s %-12s %-12s\n', 'Unit', 'Tone', 'Trace', 'Airpuff');
for u = 1:numUnits
    fprintf('%-8d %-12s %-12s %-12s\n', u, ...
        num2str(respLatency(u,2), '%.3f'), ...
        num2str(respLatency(u,3), '%.3f'), ...
        num2str(respLatency(u,4), '%.3f'));
end

%% ISI distributions — spontaneous and task separately
% figure('Name', 'ISI Distributions', 'NumberTitle', 'off');
% nCols = ceil(sqrt(numUnits));
% nRows = ceil(numUnits / nCols);
% 
% for u = 1:numUnits
%     subplot(nRows, nCols, u);
%     hold on;
%     if numel(spikeTimestamps_spont{u}) > 1
%         isi_spont = diff(sort(spikeTimestamps_spont{u})) * 1000;
%         histogram(isi_spont, 'BinWidth', 5, 'Normalization', 'probability', ...
%             'FaceColor', [0.4 0.4 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
%     end
%     if numel(spikeTimestamps_task{u}) > 1
%         isi_task = diff(sort(spikeTimestamps_task{u})) * 1000;
%         histogram(isi_task, 'BinWidth', 5, 'Normalization', 'probability', ...
%             'FaceColor', [0.8 0.4 0.4], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
%     end
%     xlabel('ISI (ms)'); ylabel('Probability');
%     title(sprintf('Unit %d', u));
%     xlim([0, 500]);
%     legend('Spont', 'Task', 'Location', 'northeast');
%     hold off;
% end
% sgtitle('ISI Distributions');

%% ISI distributions and CV analysis
cv_spont  = nan(numUnits, 1);
cv_task   = nan(numUnits, 1);
cv2_spont = nan(numUnits, 1);
cv2_task  = nan(numUnits, 1);

for u = 1:numUnits
    if numel(spikeTimestamps_spont{u}) > 2
        isi_spont = diff(sort(spikeTimestamps_spont{u}));
        cv_spont(u) = std(isi_spont) / mean(isi_spont);
    end
    if numel(spikeTimestamps_spont{u}) > 3
        isi_spont = diff(sort(spikeTimestamps_spont{u}));
        cv2_vals  = 2 * abs(diff(isi_spont)) ./ (isi_spont(1:end-1) + isi_spont(2:end));
        cv2_spont(u) = mean(cv2_vals);
    end
    if numel(spikeTimestamps_task{u}) > 2
        isi_task = diff(sort(spikeTimestamps_task{u}));
        cv_task(u) = std(isi_task) / mean(isi_task);
    end
    if numel(spikeTimestamps_task{u}) > 3
        isi_task = diff(sort(spikeTimestamps_task{u}));
        cv2_vals = 2 * abs(diff(isi_task)) ./ (isi_task(1:end-1) + isi_task(2:end));
        cv2_task(u) = mean(cv2_vals);
    end
end

% Plot ISI histograms with CV and CV2 annotations
figure('Name', 'ISI Distributions', 'NumberTitle', 'off');
nCols = ceil(sqrt(numUnits));
nRows = ceil(numUnits / nCols);

for u = 1:numUnits
    subplot(nRows, nCols, u);
    hold on;
    if numel(spikeTimestamps_spont{u}) > 2
        isi_spont_ms = diff(sort(spikeTimestamps_spont{u})) * 1000;
        histogram(isi_spont_ms, 'BinWidth', 5, 'Normalization', 'probability', ...
            'FaceColor', [0.4 0.4 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    end
    if numel(spikeTimestamps_task{u}) > 2
        isi_task_ms = diff(sort(spikeTimestamps_task{u})) * 1000;
        histogram(isi_task_ms, 'BinWidth', 5, 'Normalization', 'probability', ...
            'FaceColor', [0.8 0.4 0.4], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    end
    xlabel('ISI (ms)'); ylabel('Probability');
    title(sprintf('Unit %d\nCV=%.2f/%.2f  CV2=%.2f/%.2f', u, ...
        cv_spont(u), cv_task(u), cv2_spont(u), cv2_task(u)));
    xlim([0 500]);
    if u == 1
        legend('Spont', 'Task', 'Location', 'northeast');
    end
    hold off;
end
sgtitle('ISI Distributions (spont=blue, task=red)');

% Summary table
fprintf('\n--- ISI CV Summary ---\n');
fprintf('%-8s %-12s %-12s %-12s %-12s %-25s\n', ...
    'Unit', 'CV Spont', 'CV Task', 'CV2 Spont', 'CV2 Task', 'Putative Type');
for u = 1:numUnits
    cv_use  = cv_task(u);
    cv2_use = cv2_task(u);
    if isnan(cv_use),  cv_use  = cv_spont(u);  end
    if isnan(cv2_use), cv2_use = cv2_spont(u); end

    if isnan(cv_use)
        cell_type = 'insufficient data';
    elseif cv_use < 0.5
        cell_type = 'regular (interneuron?)';
    elseif cv_use > 1.0 && cv2_use > 0.8
        cell_type = 'bursty (pyramidal?)';
    elseif cv_use > 1.0 && cv2_use <= 0.8
        cell_type = 'non-stationary rate';
    else
        cell_type = 'irregular (pyramidal?)';
    end

    fprintf('%-8d %-12s %-12s %-12s %-12s %-25s\n', u, ...
        num2str(cv_spont(u),  '%.2f'), ...
        num2str(cv_task(u),   '%.2f'), ...
        num2str(cv2_spont(u), '%.2f'), ...
        num2str(cv2_task(u),  '%.2f'), ...
        cell_type);
end
%% ADD SPIKE WIDTH FOR CELL TYPE CLASSIFICATION


%% Burst analysis and spreadsheet export
% Max interval method parameters (standard cortical thresholds)
burst_isi_start  = 0.010;  % 10ms — max ISI to start/continue a burst
burst_isi_end    = 0.100;  % 100ms — ISI that terminates a burst
burst_min_spikes = 2;      % minimum spikes to count as a burst

% ISI feature and burst metrics storage
isi_min        = nan(numUnits, 2);  % cols: spont, task
isi_peak       = nan(numUnits, 2);
isi_p75        = nan(numUnits, 2);
isi_p95        = nan(numUnits, 2);
burst_rate     = nan(numUnits, 2);  % bursts per second
burst_mean_dur = nan(numUnits, 2);  % mean burst duration (s)
burst_mean_spk = nan(numUnits, 2);  % mean spikes per burst
burst_mean_ifr = nan(numUnits, 2);  % mean intra-burst firing rate (Hz)
burst_mean_ibi = nan(numUnits, 2);  % mean inter-burst interval (s)

for u = 1:numUnits
    spkTrains = {spikeTimestamps_spont{u}, spikeTimestamps_task{u}};
    recDurs   = [(spont_end_samp - spont_start_samp) / fs_npx, ...
                 (task_end_samp  - spont_start_samp) / fs_npx];

    for s = 1:2
        spk = sort(spkTrains{s});
        if numel(spk) < 2, continue; end

        isi = diff(spk);

        % ISI features
        isi_min(u,s)  = prctile(isi, 1) * 1000;   % in ms
        isi_peak(u,s) = mode(round(isi * 1000));   % in ms, rounded to nearest ms
        isi_p75(u,s)  = prctile(isi, 75) * 1000;  % in ms
        isi_p95(u,s)  = prctile(isi, 95) * 1000;  % in ms

        % Burst detection — max interval method
        burst_start_idx = [];
        burst_end_idx   = [];
        in_burst = false;
        b_start  = 1;

        for k = 1:numel(isi)
            if isi(k) <= burst_isi_start
                if ~in_burst
                    in_burst = true;
                    b_start  = k;
                end
            else
                if in_burst
                    if isi(k) > burst_isi_end
                        % Burst ends at spike k
                        b_end = k;
                        if (b_end - b_start + 1) >= burst_min_spikes
                            burst_start_idx(end+1) = b_start;
                            burst_end_idx(end+1)   = b_end;
                        end
                        in_burst = false;
                    end
                end
            end
        end
        % Close any open burst at end of train
        if in_burst
            b_end = numel(spk);
            if (b_end - b_start + 1) >= burst_min_spikes
                burst_start_idx(end+1) = b_start;
                burst_end_idx(end+1)   = b_end;
            end
        end

        nBursts = numel(burst_start_idx);
        if nBursts == 0, continue; end

        % Burst metrics
        dur  = zeros(nBursts, 1);
        nspk = zeros(nBursts, 1);
        ifr  = zeros(nBursts, 1);
        for b = 1:nBursts
            idx      = burst_start_idx(b) : burst_end_idx(b);
            dur(b)   = spk(burst_end_idx(b)) - spk(burst_start_idx(b));
            nspk(b)  = numel(idx) + 1;
            ifr(b)   = (nspk(b) - 1) / dur(b);
        end

        burst_rate(u,s)     = nBursts / recDurs(s);
        burst_mean_dur(u,s) = mean(dur)  * 1000;  % ms
        burst_mean_spk(u,s) = mean(nspk);
        burst_mean_ifr(u,s) = mean(ifr);

        % Inter-burst intervals
        if nBursts > 1
            ibi = zeros(nBursts-1, 1);
            for b = 1:nBursts-1
                ibi(b) = spk(burst_start_idx(b+1)) - spk(burst_end_idx(b));
            end
            burst_mean_ibi(u,s) = mean(ibi) * 1000;  % ms
        end
    end
end

%% Assemble and export spreadsheet
unitNums = (1:numUnits)';

% Classification and latency
tone_class    = {unitClass.Tone}';
trace_class   = {unitClass.Trace}';
airpuff_class = {unitClass.Airpuff}';
tone_p        = [unitClass.Tone_p]';
trace_p       = [unitClass.Trace_p]';
airpuff_p     = [unitClass.Airpuff_p]';
tone_test     = {unitClass.Tone_test}';
trace_test    = {unitClass.Trace_test}';
airpuff_test  = {unitClass.Airpuff_test}';
lat_tone      = respLatency(:,2) * 1000;   % ms
lat_trace     = respLatency(:,3) * 1000;   % ms
lat_airpuff   = respLatency(:,4) * 1000;   % ms

% CV metrics
outT = table(unitNums, ...
    ... % Classification
    tone_class, tone_p, tone_test, ...
    trace_class, trace_p, trace_test, ...
    airpuff_class, airpuff_p, airpuff_test, ...
    ... % Response latency (ms)
    lat_tone, lat_trace, lat_airpuff, ...
    ... % CV metrics
    cv_spont, cv_task, cv2_spont, cv2_task, ...
    ... % ISI features — spontaneous
    isi_min(:,1),  isi_peak(:,1), isi_p75(:,1), isi_p95(:,1), ...
    ... % ISI features — task
    isi_min(:,2),  isi_peak(:,2), isi_p75(:,2), isi_p95(:,2), ...
    ... % Burst metrics — spontaneous
    burst_rate(:,1), burst_mean_dur(:,1), burst_mean_spk(:,1), ...
    burst_mean_ifr(:,1), burst_mean_ibi(:,1), ...
    ... % Burst metrics — task
    burst_rate(:,2), burst_mean_dur(:,2), burst_mean_spk(:,2), ...
    burst_mean_ifr(:,2), burst_mean_ibi(:,2), ...
    'VariableNames', { ...
    'Unit', ...
    'Tone_Class', 'Tone_p', 'Tone_Test', ...
    'Trace_Class', 'Trace_p', 'Trace_Test', ...
    'Airpuff_Class', 'Airpuff_p', 'Airpuff_Test', ...
    'Latency_Tone_ms', 'Latency_Trace_ms', 'Latency_Airpuff_ms', ...
    'CV_Spont', 'CV_Task', 'CV2_Spont', 'CV2_Task', ...
    'ISI_Min_Spont_ms', 'ISI_Peak_Spont_ms', 'ISI_P75_Spont_ms', 'ISI_P95_Spont_ms', ...
    'ISI_Min_Task_ms',  'ISI_Peak_Task_ms',  'ISI_P75_Task_ms',  'ISI_P95_Task_ms', ...
    'BurstRate_Spont',    'BurstDur_Spont_ms',    'BurstSpkCount_Spont', ...
    'BurstIFR_Spont_Hz',  'BurstIBI_Spont_ms', ...
    'BurstRate_Task',     'BurstDur_Task_ms',     'BurstSpkCount_Task', ...
    'BurstIFR_Task_Hz',   'BurstIBI_Task_ms'});

unitStatsFile = fullfile(pth, [baseName '_unit_stats.xlsx']);
writetable(outT, unitStatsFile, 'Sheet', 'UnitStats');
fprintf('Unit stats saved to %s\n', unitStatsFile);

%% Save workspace
save(fullfile(pth, [baseName '_workspace.mat']));
fprintf('Workspace saved to %s\n', fullfile(pth, [baseName '_workspace.mat']));