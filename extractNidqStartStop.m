function edges = extractNidqStartStop(nidqBinPath, varargin)
% npxTimestamps  Extract 4 TTL transition times (R1,F1,R2,F2) from a Neuropixels SpikeGLX *.nidq.bin file.
%
% USAGE
%   edges = npxTimestamps("C:\...\something_g0_t0.nidq.bin");
%   edges = npxTimestamps("C:\...\something_g0_t0.nidq.bin", "Channel", 6);
%
% WHAT IT DOES
%   - Reads SpikeGLX nidq.bin (int16 interleaved) via memmap.
%   - Finds TTL edges on a selected channel (or auto-selects).
%   - Returns all rise/fall edges AND the first two pulses => 4 edges: R1,F1,R2,F2.
%   - Plots overview of all channels + TTL channel with edge markers.
%
% KEY OPTIONS
%   "Channel"          : TTL channel number (MATLAB 1-based). If empty => auto-pick.
%   "AutoPickChannel"  : true/false (default true if Channel not provided)
%   "Threshold"        : manual threshold in int16 counts (default auto from percentiles)
%   "MinStateMs"       : debounce duration in ms (default 5 ms)
%   "PlotAllChannels"  : true/false (default true)
%   "PlotTTL"          : true/false (default true)
%   "OverviewHz"       : downsample target rate for plots (default 1000 Hz)
%
% OUTPUT STRUCT edges
%   edges.fsHz
%   edges.nSavedChans
%   edges.nSamples
%   edges.chanNumber
%   edges.thresholdUsed
%   edges.riseSamples, edges.fallSamples
%   edges.riseTimesS,  edges.fallTimesS
%   edges.fourEdgeSamples [R1 F1 R2 F2]
%   edges.fourEdgeTimesS  [R1 F1 R2 F2]
%   edges.metaPath

% -----------------------------
% Parse inputs
% -----------------------------
p = inputParser;
p.addRequired('nidqBinPath', @(x) ischar(x) || isstring(x));

p.addParameter('Channel', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>=1));
p.addParameter('AutoPickChannel', [], @(x) isempty(x) || islogical(x) || (isnumeric(x)&&isscalar(x)));
p.addParameter('Threshold', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('MinStateMs', 5, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('PlotAllChannels', true, @(x) islogical(x) || (isnumeric(x)&&isscalar(x)));
p.addParameter('PlotTTL', true, @(x) islogical(x) || (isnumeric(x)&&isscalar(x)));
p.addParameter('OverviewHz', 1000, @(x) isnumeric(x) && isscalar(x) && x>0);

p.parse(nidqBinPath, varargin{:});

nidqBinPath = char(p.Results.nidqBinPath);
chanNumber  = p.Results.Channel;
thrManual   = p.Results.Threshold;
minStateMs  = double(p.Results.MinStateMs);
plotAll     = logical(p.Results.PlotAllChannels);
plotTTL     = logical(p.Results.PlotTTL);
overviewHz  = double(p.Results.OverviewHz);

autoPick = p.Results.AutoPickChannel;
if isempty(autoPick)
    autoPick = isempty(chanNumber);
else
    autoPick = logical(autoPick);
end

assert(exist(nidqBinPath,'file')==2, 'nidq.bin not found: %s', nidqBinPath);

% Meta path
[folder, base, ~] = fileparts(nidqBinPath);
metaPath = fullfile(folder, [base '.meta']);
assert(exist(metaPath,'file')==2, 'Missing .meta file next to nidq.bin: %s', metaPath);

% -----------------------------
% Read meta: fs and nSavedChans
% -----------------------------
meta = readSpikeGLXMeta(metaPath);

assert(isfield(meta,'niSampRate'), 'Meta missing niSampRate.');
fs = str2double(meta.niSampRate);

nSaved = [];
if isfield(meta,'nSavedChans')
    nSaved = str2double(meta.nSavedChans);
elseif isfield(meta,'snsSaveChanSubset')
    nSaved = countChanSubset(meta.snsSaveChanSubset);
else
    error('Could not determine number of saved channels (nSavedChans or snsSaveChanSubset).');
end
assert(~isnan(nSaved) && nSaved>=1, 'Invalid nSavedChans.');

% File length -> samples
d = dir(nidqBinPath);
nInt16 = d.bytes/2;
assert(mod(nInt16, nSaved)==0, 'File length not divisible by nSavedChans. Check meta.');
nSamples = nInt16 / nSaved;

% -----------------------------
% Memmap and reshape
% -----------------------------
raw = memmapfile(nidqBinPath, 'Format', 'int16');
X = reshape(raw.Data, [nSaved, nSamples]);   % channels x samples

% -----------------------------
% Auto-pick TTL channel if needed
% -----------------------------
if autoPick
    [bestChan, scoreTable] = autoPickTtlChannel(X, fs, overviewHz);
    fprintf('\nAuto-pick TTL channel => %d\n', bestChan);
    disp(scoreTable);
    chanNumber = bestChan;
else
    assert(~isempty(chanNumber), 'Channel must be provided if AutoPickChannel is false.');
end

assert(chanNumber>=1 && chanNumber<=nSaved, 'Channel out of range. nSaved=%d', nSaved);

% Extract channel trace
x = double(X(chanNumber, :));

% Debounce in samples
minStateSamples = max(1, round((minStateMs/1000) * fs));

% -----------------------------
% Detect ALL edges (robust)
% -----------------------------
[riseSamp, fallSamp, thrUsed, debug] = detectAllEdgesTTL(x, thrManual, minStateSamples);

% Pair first two pulses => 4 time points
[fourSamp, fourTime] = pickFirstTwoPulses(riseSamp, fallSamp, fs);

% -----------------------------
% Build output struct
% -----------------------------
edges = struct();
edges.fsHz = fs;
edges.nSavedChans = nSaved;
edges.nSamples = nSamples;
edges.chanNumber = chanNumber;
edges.thresholdUsed = thrUsed;
edges.minStateSamples = minStateSamples;

edges.riseSamples = riseSamp;
edges.fallSamples = fallSamp;

edges.riseTimesS = (riseSamp-1)/fs;
edges.fallTimesS = (fallSamp-1)/fs;

edges.fourEdgeSamples = fourSamp;     % [R1 F1 R2 F2]
edges.fourEdgeTimesS  = fourTime;     % seconds
edges.metaPath = metaPath;

edges.debug = debug;

% -----------------------------
% Plots
% -----------------------------
ds = max(1, round(fs/overviewHz)); % downsample for plotting

if plotAll
    plotAllChannelsOverview(X, fs, ds);
end

if plotTTL
    plotTtlWithEdges(x, fs, ds, thrUsed, edges.riseTimesS, edges.fallTimesS, edges.fourEdgeTimesS, chanNumber);
end

end % end main function


% =========================================================================
% Helpers
% =========================================================================

function meta = readSpikeGLXMeta(metaPath)
% Read key=value pairs from SpikeGLX .meta into a struct
meta = struct();
fid = fopen(metaPath,'r');
assert(fid>0, 'Failed to open meta: %s', metaPath);
cleanup = onCleanup(@() fclose(fid));

while true
    t = fgetl(fid);
    if ~ischar(t), break; end
    t = strtrim(t);
    if isempty(t) || startsWith(t,'#'), continue; end
    eq = strfind(t,'=');
    if isempty(eq), continue; end
    k = strtrim(t(1:eq(1)-1));
    v = strtrim(t(eq(1)+1:end));
    k = matlab.lang.makeValidName(k);
    meta.(k) = v;
end
end

function n = countChanSubset(subsetStr)
% Parse snsSaveChanSubset like "0:7,9,11"
parts = split(string(subsetStr), ",");
count = 0;
for i = 1:numel(parts)
    p = strtrim(parts(i));
    if strlength(p)==0, continue; end
    if contains(p, ":")
        ab = split(p, ":");
        a = str2double(ab(1)); b = str2double(ab(2));
        count = count + (b - a + 1);
    else
        count = count + 1;
    end
end
n = count;
end

function [bestChan, scoreTable] = autoPickTtlChannel(X, fs, overviewHz)
% Heuristic: TTL-like channel usually has huge two-level range (p99 - p1).
% X is [nSaved x nSamples].
[nSaved, nSamples] = size(X);
ds = max(1, round(fs/overviewHz));
idx = 1:ds:nSamples;

p1s = zeros(nSaved,1);
p99s = zeros(nSaved,1);
ranges = zeros(nSaved,1);

for ch = 1:nSaved
    x = double(X(ch, idx));
    p1 = prctile(x, 1);
    p99 = prctile(x, 99);
    p1s(ch) = p1;
    p99s(ch) = p99;
    ranges(ch) = p99 - p1;
end

[~, bestChan] = max(ranges);

scoreTable = table((1:nSaved)', p1s, p99s, ranges, ...
    'VariableNames', {'Chan','P1','P99','Range'});
scoreTable = sortrows(scoreTable, 'Range', 'descend');
end

function [riseSamp, fallSamp, thrUsed, debug] = detectAllEdgesTTL(x, thrManual, minStateSamples)
% Robust edge detection for two-level (TTL-like) analog trace:
% - Auto threshold = midpoint between p1 and p99
% - Debounce via "stable for N samples" using convolution

x = double(x(:))'; % ensure row

p1  = prctile(x, 1);
p99 = prctile(x, 99);

if isempty(thrManual)
    thrUsed = (p1 + p99)/2;
else
    thrUsed = thrManual;
end

isHigh = x > thrUsed;

% Debounce: stable high/low for minStateSamples
kernel = ones(1, minStateSamples);
highStable = conv(double(isHigh), kernel, 'same') >= minStateSamples;
lowStable  = conv(double(~isHigh), kernel, 'same') >= minStateSamples;

state = isHigh;
state(highStable) = 1;
state(lowStable)  = 0;

riseSamp = find(state(1:end-1)==0 & state(2:end)==1) + 1;
fallSamp = find(state(1:end-1)==1 & state(2:end)==0) + 1;

debug = struct();
debug.min = min(x);
debug.max = max(x);
debug.p1 = p1;
debug.p99 = p99;
debug.highFrac = mean(isHigh);

fprintf('\nTTL debug: min=%.0f max=%.0f p1=%.1f p99=%.1f thr=%.1f highFrac=%.3f minStateSamples=%d\n', ...
    debug.min, debug.max, debug.p1, debug.p99, thrUsed, debug.highFrac, minStateSamples);
end

function [fourSamp, fourTime] = pickFirstTwoPulses(rises, falls, fs)
% Pair each rise with the first subsequent fall.
% Return first two pairs => [R1 F1 R2 F2]
fourSamp = nan(1,4);
fourTime = nan(1,4);

if isempty(rises) || isempty(falls)
    return
end

pairs = zeros(0,2);
for i = 1:numel(rises)
    f = falls(find(falls > rises(i), 1, 'first'));
    if ~isempty(f)
        pairs(end+1,:) = [rises(i), f]; %#ok<AGROW>
    end
    if size(pairs,1) >= 2
        break
    end
end

if size(pairs,1) >= 2
    fourSamp = [pairs(1,1), pairs(1,2), pairs(2,1), pairs(2,2)];
    fourTime = (fourSamp-1)/fs;
end
end

function plotAllChannelsOverview(X, fs, ds)
% Plot all channels downsampled (one tile per channel)
[nSaved, nSamples] = size(X);
idx = 1:ds:nSamples;
t = (idx-1)/fs;

figure('Name','NIDQ overview (downsampled)','Color','w');
tiledlayout(nSaved,1,'TileSpacing','compact','Padding','compact');

for ch = 1:nSaved
    nexttile;
    x = double(X(ch, idx));
    plot(t, x, 'k');
    ylabel(sprintf('Ch %d', ch));
    if ch ~= nSaved
        set(gca, 'XTickLabel', []);
    end
    if ch == 1
        title(sprintf('All NIDQ channels (downsampled to ~%.0f Hz)', fs/ds));
    end
end
xlabel('Time (s)');
end

function plotTtlWithEdges(x, fs, ds, thr, riseTimes, fallTimes, fourTimes, chanNumber)
% Plot TTL channel downsampled with threshold and edge markers
nSamples = numel(x);
idx = 1:ds:nSamples;
t = (idx-1)/fs;
x_ds = x(idx);

figure('Name',sprintf('TTL channel %d', chanNumber),'Color','w');
plot(t, x_ds, 'k'); grid on; hold on;
yline(thr, 'r--', 'Threshold');

% Mark rises/falls (at top of plot for visibility)
yTop = max(x_ds);
plot(riseTimes, yTop*ones(size(riseTimes)), 'g^', 'MarkerFaceColor','g', 'MarkerSize',6);
plot(fallTimes, yTop*0.98*ones(size(fallTimes)), 'mv', 'MarkerFaceColor','m', 'MarkerSize',6);

% Emphasize the 4 key edges if present
if ~isempty(fourTimes) && all(~isnan(fourTimes))
    for k = 1:4
        xline(fourTimes(k), 'b', sprintf('E%d', k), 'LineWidth', 1.5);
    end
end

xlabel('Time (s)');
ylabel('int16 counts');
title(sprintf('TTL channel %d (downsampled to ~%.0f Hz)', chanNumber, fs/ds));
ylim([min(x_ds)-500, max(x_ds)+500]);
legend({'signal','threshold','rises','falls'}, 'Location','best');
end
