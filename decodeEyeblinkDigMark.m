function decoded = decodeEyeblinkDigMark(out)
%DECODEEYEBLINKDIGMARK Decode DigMark codes for your trace eyeblink experiment.
%
% Updated code meanings:
%   01 -> spontaneous activity Neuropixels recording START
%   02 -> (within spontaneous recording; ignored here)
%   03 -> spontaneous activity Neuropixels recording END
%   04 -> camera ON (trial starts camera opening)
%   05 -> tone ON (after a delay)
%   06 -> tone OFF + trace period START
%   07 -> air puff ON
%   08 -> air puff OFF
%   09 -> camera OFF (trial ends)
%
% Inputs:
%   out.time_sec (preferred) or out.time_ticks
%   out.code
%
% Outputs:
%   decoded.spontaneous_segments : Nx2 [start end] from 01->03 pairs
%   decoded.events_table         : table(time, code, label)
%   decoded.trials               : struct array; trials segmented by 04->09
%   decoded.qc                   : summary flags

    % ----- choose timebase
    if isfield(out,'time_sec') && ~isempty(out.time_sec)
        t = out.time_sec(:);
        tUnits = 's';
    else
        t = double(out.time_ticks(:));
        tUnits = 'ticks';
        warning('out.time_sec is empty. Using ticks instead of seconds.');
    end

    code = double(out.code(:));

    % ----- sort by time
    [t, idx] = sort(t);
    code = code(idx);

    % ----- label map
    labels = strings(size(code));
    labels(code==1) = "spontaneous_npx_start";
    labels(code==2) = "spontaneous_npx_mid";
    labels(code==3) = "spontaneous_npx_end";
    labels(code==4) = "camera_on";
    labels(code==5) = "tone_on";
    labels(code==6) = "tone_off_trace_start";
    labels(code==7) = "airpuff_on";
    labels(code==8) = "airpuff_off";
    labels(code==9) = "camera_off";
    labels(labels=="") = "unknown";

    decoded = struct();
    decoded.time_units = tUnits;

    decoded.events_table = table(t, code, labels, ...
        'VariableNames', {'time','code','label'});

    % ----- spontaneous segments: pair each 01 with the next 03
    spStart = t(code==1);
    spEnd   = t(code==3);
    decoded.spontaneous_segments = pairStartEnd(spStart, spEnd);

    % ----- trials: pair each 04 with the next 09
    camOn  = t(code==4);
    camOff = t(code==9);
    [trialSeg, segNotes] = pairStartEndWithNotes(camOn, camOff);

    trials = struct( ...
        'trial_idx', {}, ...
        'camera_on', {}, 'camera_off', {}, ...
        'tone_on', {}, 'tone_off', {}, 'trace_start', {}, ...
        'airpuff_on', {}, 'airpuff_off', {}, ...
        'cam_dur', {}, 'tone_dur', {}, 'trace_dur', {}, 'puff_dur', {}, ...
        'cam_to_tone', {}, 'tone_to_puff', {}, ...
        'notes', {} );

    for i = 1:size(trialSeg,1)
        tOn  = trialSeg(i,1);
        tOff = trialSeg(i,2);

        inWin = (t >= tOn) & (t <= tOff);
        tw = t(inWin);
        cw = code(inWin);

        this = struct();
        this.trial_idx  = i;
        this.camera_on  = tOn;
        this.camera_off = tOff;

        this.tone_on     = firstOrNaN(tw(cw==5));
        this.tone_off    = firstOrNaN(tw(cw==6));
        this.trace_start = this.tone_off; % alias
        this.airpuff_on  = firstOrNaN(tw(cw==7));
        this.airpuff_off = firstOrNaN(tw(cw==8));

        % Durations
        this.cam_dur  = safeDiff(this.camera_off, this.camera_on);

        this.tone_dur = NaN;
        if ~isnan(this.tone_on) && ~isnan(this.tone_off) && this.tone_off >= this.tone_on
            this.tone_dur = this.tone_off - this.tone_on;
        end

        this.trace_dur = NaN;
        if ~isnan(this.tone_off) && ~isnan(this.airpuff_on) && this.airpuff_on >= this.tone_off
            this.trace_dur = this.airpuff_on - this.tone_off;
        end

        this.puff_dur = NaN;
        if ~isnan(this.airpuff_on) && ~isnan(this.airpuff_off) && this.airpuff_off >= this.airpuff_on
            this.puff_dur = this.airpuff_off - this.airpuff_on;
        end

        % Key intervals
        this.cam_to_tone = NaN;
        if ~isnan(this.tone_on)
            this.cam_to_tone = this.tone_on - this.camera_on;
        end

        this.tone_to_puff = NaN;
        if ~isnan(this.tone_on) && ~isnan(this.airpuff_on)
            this.tone_to_puff = this.airpuff_on - this.tone_on;
        end

        % Notes / QC
        notes = strings(0,1);
        if strlength(segNotes(i)) > 0, notes(end+1) = segNotes(i); end

        if isnan(this.tone_on),  notes(end+1) = "missing tone_on (05)"; end
        if isnan(this.tone_off), notes(end+1) = "missing tone_off/trace_start (06)"; end
        if ~isnan(this.tone_on) && ~isnan(this.tone_off) && this.tone_off < this.tone_on
            notes(end+1) = "tone_off (06) before tone_on (05)";
        end

        if isnan(this.airpuff_on),  notes(end+1) = "missing airpuff_on (07)"; end
        if isnan(this.airpuff_off), notes(end+1) = "missing airpuff_off (08)"; end
        if ~isnan(this.airpuff_on) && ~isnan(this.airpuff_off) && this.airpuff_off < this.airpuff_on
            notes(end+1) = "airpuff_off (08) before airpuff_on (07)";
        end

        if ~isnan(this.tone_off) && ~isnan(this.airpuff_on) && this.airpuff_on < this.tone_off
            notes(end+1) = "airpuff_on (07) before trace_start (06)";
        end

        this.notes = strjoin(notes, "; ");
        trials(end+1) = this; %#ok<AGROW>
    end

    decoded.trials = trials;

    % ----- QC summary
    qc = struct();
    qc.n_events = numel(code);
    qc.n_trials = numel(trials);
    qc.n_spontaneous_segments = size(decoded.spontaneous_segments,1);

    qc.overlapping_trials = false;
    if size(trialSeg,1) >= 2
        qc.overlapping_trials = any(trialSeg(2:end,1) < trialSeg(1:end-1,2));
    end

    decoded.qc = qc;
end

% ---------- helpers ----------
function seg = pairStartEnd(starts, ends)
    starts = starts(:);
    ends   = ends(:);
    seg = zeros(0,2);
    if isempty(starts) || isempty(ends), return; end

    j = 1;
    for i = 1:numel(starts)
        while j <= numel(ends) && ends(j) <= starts(i)
            j = j + 1;
        end
        if j <= numel(ends)
            seg(end+1,:) = [starts(i), ends(j)]; %#ok<AGROW>
            j = j + 1;
        end
    end
end

function [seg, notes] = pairStartEndWithNotes(starts, ends)
    starts = starts(:);
    ends   = ends(:);

    seg = zeros(numel(starts),2);
    notes = strings(numel(starts),1);

    j = 1;
    for i = 1:numel(starts)
        seg(i,1) = starts(i);

        while j <= numel(ends) && ends(j) <= starts(i)
            j = j + 1;
        end

        if j <= numel(ends)
            seg(i,2) = ends(j);
            j = j + 1;
        else
            seg(i,2) = starts(i); % degenerate
            notes(i) = "missing camera_off (09) after camera_on (04)";
        end
    end
end

function v = firstOrNaN(x)
    if isempty(x), v = NaN; else, v = x(1); end
end

function d = safeDiff(a,b)
    if isnan(a) || isnan(b), d = NaN; else, d = a - b; end
end
