function out = extractDigMark32_CEDS64(smrFile, cedPath)
% out = extractDigMark32_CEDS64(smrFile, cedPath)
%
% Reads DigMark/Marker channel 32 from a Spike2 .smr file using CEDS64ML.
%
% Inputs
%   smrFile : e.g. 'fullTest013026.smr'
%   cedPath : folder containing CEDS64ML m-files (the download from CED)
%
% Output struct fields
%   out.time_ticks : int64 marker times in ticks
%   out.time_sec   : marker times in seconds (if timebase available)
%   out.code       : marker code (0..7 expected) extracted from code1 by default
%   out.codes4     : Nx4 marker bytes [code1 code2 code3 code4]

    if nargin < 2 || isempty(cedPath)
        error('Provide cedPath (folder containing CEDS64ML .m files, including CEDS64Open.m, CEDS64LoadLib.m).');
    end

    addpath(genpath(cedPath));

    % ---- Load the CEDS64ML DLL (ceds64int)
    % This is the standard entry-point recommended by CED / many users. :contentReference[oaicite:0]{index=0}
    if exist('CEDS64LoadLib','file') ~= 2
        error('CEDS64LoadLib.m not found. Check cedPath.');
    end
    CEDS64LoadLib(cedPath);

    % ---- Open file
    if exist('CEDS64Open','file') ~= 2
        error('CEDS64Open.m not found. Check cedPath.');
    end
    fhand = CEDS64Open(smrFile, 1); % 1 = read-only in most CED examples
    if fhand < 0
        error('CEDS64Open failed with code %d', fhand);
    end

    chan = 32;

    % ---- Read markers
    % CEDS64ReadMarkers(fhand, iChan, iN, i64From, i64To, maskh)
    % i64To = -1 means "to end of channel" per the file header comment you uploaded.
    iN      = 1000000;   % max markers to read in one call (bump up if needed)
    i64From = int64(0);
    i64To   = int64(-1);

    [nRead, markers] = CEDS64ReadMarkers(fhand, chan, iN, i64From, i64To);
    if nRead <= 0 || isempty(markers)
        CEDS64Close(fhand);
        error('No markers read from channel %d (nRead=%d). Check that chan 32 is a marker/DigMark channel.', chan, nRead);
    end

    % ---- Extract time + the 4 marker bytes
    % CEDS64ReadMarkers builds CEDMarker objects and sets:
    %   SetTime(m_Time), SetCode(1..4, m_Code1..4)
    n = numel(markers);

    time_ticks = zeros(n,1,'int64');
    codes4     = zeros(n,4,'uint8');

    for i = 1:n
        mk = markers(i);

        % Time
        if ismethod(mk,'GetTime')
            time_ticks(i) = int64(mk.GetTime());
        else
            % Fallbacks (depending on class version)
            try
                time_ticks(i) = int64(mk.m_Time);
            catch
                error('Could not read marker time. Your CEDMarker class lacks GetTime and m_Time.');
            end
        end

        % Codes 1..4
        for k = 1:4
            if ismethod(mk,'GetCode')
                codes4(i,k) = uint8(mk.GetCode(k));
            else
                % Fallback guesses if class exposes properties
                try
                    codes4(i,1) = uint8(mk.m_Code1);
                    codes4(i,2) = uint8(mk.m_Code2);
                    codes4(i,3) = uint8(mk.m_Code3);
                    codes4(i,4) = uint8(mk.m_Code4);
                    break;
                catch
                    error('Could not read marker codes. Your CEDMarker class lacks GetCode and m_Code1..4.');
                end
            end
        end
    end

    % The DigMark code is *usually* in Code1
    code = double(codes4(:,1));

    % ---- Convert ticks -> seconds (if we can find timebase function)
    time_sec = [];
    if exist('CEDS64TimeBase','file') == 2
        % many CEDS64ML installs include this helper
        tb = CEDS64TimeBase(fhand);   % seconds per tick (typical)
        time_sec = double(time_ticks) * double(tb);
    else
        % leave in ticks if timebase helper not present
        time_sec = [];
    end

    % ---- Close file
    if exist('CEDS64Close','file') == 2
        CEDS64Close(fhand);
    else
        % Some versions name it CEDS64CloseFile; try both
        if exist('CEDS64CloseFile','file') == 2
            CEDS64CloseFile(fhand);
        end
    end

    out = struct();
    out.time_ticks = time_ticks;
    out.time_sec   = time_sec;
    out.code       = code;     % expected 0..7
    out.codes4     = codes4;

    % Helpful warning if code isn't in byte 1
    if any(out.code > 9)
        warning(['Some codes > 9 found in codes4(:,1). ' ...
                 'Check which column contains [0..7] by running unique(out.codes4(:,k)).']);
    end
end
