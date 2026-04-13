function edges = npxTimestamps

nidqPath = uigetfile('*.bin', 'Select nidq.bin file');

% Let it auto-pick the TTL channel (recommended first time)
edges = extractNidqStartStop(nidqPath, "Channel", 6);