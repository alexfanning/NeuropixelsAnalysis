#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Import the necessary libraries
import spikeinterface.full as si
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import logging


# Set up the base and SpikeGLX folder paths
base_folder = Path(r'C:\Alex\project\ilnStim\data\test3\test3diIRightHemiStraight121825_g0')
base_folder.mkdir(parents=True, exist_ok=True)

spikeglx_folder = base_folder / 'test3diIRightHemiStraight121825_g0_imec0'

# Logging setup (with force=True to override existing logging)
log_file = base_folder / 'postprocessing_log.txt'
logging.basicConfig(
    filename=str(log_file),
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    force=True
)

logging.info("---- Logging started ----")

# Get the available streams in the SpikeGLX folder
stream_names, stream_ids = si.get_neo_streams('spikeglx', spikeglx_folder)
logging.info(f"Available stream names: {stream_names}")

# Load the specific AP stream (we are not loading the sync channel)
raw_rec = si.read_spikeglx(spikeglx_folder, stream_name='imec0.ap', load_sync_channel=False)
logging.info("Raw recording loaded successfully.")

# Print a summary of the recording
logging.info(f"Raw recording summary: {raw_rec}")

# Display the probe information as a DataFrame
logging.info(f"Probe information: \n{raw_rec.get_probe().to_dataframe()}")
logging.info("Section 1 completed successfully.")


# In[2]:


# Section 2: Preprocessing (Filtering, Removing Bad Channels, etc.)
logging.info("---- Section 2: Preprocessing ----")

# Apply high-pass filter
rec1 = si.highpass_filter(raw_rec, freq_min=300.)
logging.info("High-pass filter applied.")

# Detect and remove bad channels
bad_channel_ids, _ = si.detect_bad_channels(rec1)
rec2 = rec1.remove_channels(bad_channel_ids)
logging.info(f"Removed bad channels with IDs: {bad_channel_ids}")

# Further preprocessing: Phase shift and common reference
rec3 = si.phase_shift(rec2)
rec4 = si.common_reference(rec3, operator="median", reference="global")
logging.info("Phase shift and common reference applied.")

# Save the preprocessed data
job_kwargs = dict(n_jobs=40, chunk_duration='1s', progress_bar=True)
rec4.save(folder=base_folder / 'preprocess', format='binary', **job_kwargs)
logging.info("Preprocessed data saved to 'preprocess' folder.")

logging.info("Section 2 completed successfully.")

# Section 3: Visualizations (Probe Map, Signal Traces)
logging.info("---- Section 3: Visualizations ----")

# Plot the probe map with channel IDs
fig, ax = plt.subplots(figsize=(15, 10))
si.plot_probe_map(raw_rec, ax=ax, with_channel_ids=True)
ax.set_ylim(-100, 100)
ax.set_title("Probe Map with Channel IDs")
ax.set_xlabel("X-axis")
ax.set_ylabel("Y-axis")
plt.show()
logging.info("Probe map visualization completed.")

# Plot signal traces using ipywidgets and matplotlib backends
get_ipython().run_line_magic('matplotlib', 'widget')
si.plot_traces({'cmr': rec4}, backend='ipywidgets')

# Static plot with matplotlib
fig, axs = plt.subplots(ncols=3, figsize=(20, 10))
si.plot_traces(rec1, backend='matplotlib', clim=(-50, 50), ax=axs[0])
si.plot_traces(rec2, backend='matplotlib', clim=(-50, 50), ax=axs[1])
si.plot_traces(rec4, backend='matplotlib', clim=(-50, 50), ax=axs[2])
for i, label in enumerate(('filter', 'cmr', 'final')):
    axs[i].set_title(label)
logging.info("Signal traces visualization completed.")

logging.info("Section 3 completed successfully.")

# Section 4: Noise Estimation and Peak Detection
logging.info("---- Section 4: Noise Estimation and Peak Detection ----")

# Estimate noise levels (scaled and raw)
noise_levels_microV = si.get_noise_levels(rec4, return_scaled=True)
noise_levels_int16 = si.get_noise_levels(rec4, return_scaled=False)

# Plot noise levels
fig, ax = plt.subplots()
ax.hist(noise_levels_microV, bins=np.arange(5, 30, 2.5))
ax.set_xlabel('Noise [microV]')
logging.info("Noise estimation and histogram plot completed.")

# Detect peaks
from spikeinterface.sortingcomponents.peak_detection import detect_peaks
peaks = detect_peaks(rec4, method='locally_exclusive', noise_levels=noise_levels_int16,
                     detect_threshold=5, radius_um=50., **job_kwargs)
logging.info("Peak detection completed.")

# Localize peaks
from spikeinterface.sortingcomponents.peak_localization import localize_peaks
peak_locations = localize_peaks(rec4, peaks, method='center_of_mass', radius_um=50., **job_kwargs)
logging.info("Peak localization completed.")

# Visualize drift
fs = rec4.sampling_frequency
fig, ax = plt.subplots(figsize=(10, 8))
ax.scatter(peaks['sample_index'] / fs, peak_locations['y'], color='k', marker='.', alpha=0.002)
logging.info("Drift visualization completed.")

# Plot probe map with peak locations
fig, ax = plt.subplots(figsize=(15, 10))
si.plot_probe_map(rec4, ax=ax, with_channel_ids=True)
ax.set_ylim(-100, 150)
ax.scatter(peak_locations['x'], peak_locations['y'], color='purple', alpha=0.002)
logging.info("Probe map with peak locations visualization completed.")

logging.info("Section 4 completed successfully.")


# In[3]:


# Section 5: Kilosort4 Sorting
params_kilosort4 = {
    "do_correction": True,
    "highpass_cutoff": 300,
    "Th_universal": 5,
    "Th_learned": 4,
    "Th_single_ch": 4,
    "nearest_chans": 16,
    "batch_size": 20000,
    "nblocks": 5
}

sorting = si.run_sorter('kilosort4', rec4,
                        folder=base_folder / 'kilosort4_output',
                        remove_existing_folder=True,
                        verbose=True,
                        **params_kilosort4)
logging.info("Kilosort4 sorting completed.")

# Read the sorted output back
sorting = si.read_sorter_folder(base_folder / 'kilosort4_output')
logging.info("Sorted output loaded.")

logging.info("Section 5 completed successfully.")


# In[4]:


# Section 6: Sorting Analyzer and Quality Metrics
logging.info("---- Section 6: Sorting Analyzer and Quality Metrics ----")

# Set global job kwargs for analysis
si.set_global_job_kwargs(n_jobs=4)

# Create sorting analyzer
analyzer = si.create_sorting_analyzer(sorting, rec4, sparse=True, format="memory")
analyzer.compute("random_spikes", method="uniform", max_spikes_per_unit=500)
analyzer.compute("waveforms", ms_before=1.5, ms_after=2.)
analyzer.compute("templates", operators=["average", "median", "std"])
analyzer.compute("noise_levels")
analyzer.compute("correlograms")
analyzer.compute("unit_locations")
analyzer.compute("spike_amplitudes")
analyzer.compute("template_similarity")
analyzer.compute("principal_components")
logging.info("Sorting analyzer computations completed.")

# Save the analyzer results
analyzer_saved = analyzer.save_as(folder=base_folder / "analyzer", format="binary_folder")
logging.info("Analyzer results saved.")

# Define the quality metrics
metric_names = [
    "firing_rate", "presence_ratio", "snr", "isi_violation", "amplitude_cutoff", "isolation_distance", "l_ratio", "d_prime", "nn_isolation", 
]
metrics = si.compute_quality_metrics(analyzer, metric_names=metric_names)
logging.info("Quality metrics computed for units.")
logging.info(f"Computed metrics columns: {metrics.columns.tolist()}")


# In[5]:


# Section: tiered quality filtering
def filter_units_by_quality(metrics_df, tier="medium"):
    thresholds = {
        "strict": {
            "snr": 6.0,
            "presence_ratio": 0.95,
            "isi_violations_ratio": 0.20,
            "amplitude_cutoff": 0.05,
            "isolation_distance": 30,
            "l_ratio": 0.05,
            "d_prime": 3.0,
            "nn_isolation": 0.95
        },
        "medium": {
            "snr": 4.5,
            "presence_ratio": 0.80,
            "isi_violations_ratio": 0.30,
            "amplitude_cutoff": 0.10,
            "isolation_distance": 15,
            "l_ratio": 0.15,
            "d_prime": 1.5,
            "nn_isolation": 0.80
        },
        "lenient": {
            "snr": 3.0,
            "presence_ratio": 0.70,
            "isi_violations_ratio": 0.50,
            "amplitude_cutoff": 0.20,
            "isolation_distance": 5,
            "l_ratio": 0.30,
            "d_prime": 1.0,
            "nn_isolation": 0.60
        }
    }
    if tier not in thresholds:
        raise ValueError(f"Unknown tier '{tier}'")
    th = thresholds[tier]
    mask = (
        (metrics_df["snr"] > th["snr"]) &
        (metrics_df["presence_ratio"] > th["presence_ratio"]) &
        (metrics_df["isi_violations_ratio"] < th["isi_violations_ratio"]) &
        (metrics_df["amplitude_cutoff"] < th["amplitude_cutoff"]) &
        (metrics_df["isolation_distance"] > th["isolation_distance"]) &
        (metrics_df["l_ratio"] < th["l_ratio"]) &
        (metrics_df["d_prime"] > th["d_prime"]) &
        (metrics_df["nn_isolation"] > th["nn_isolation"])
    )
    return metrics_df[mask].copy()

# apply for all tiers
metrics_strict  = filter_units_by_quality(metrics, tier="strict")
metrics_medium  = filter_units_by_quality(metrics, tier="medium")
metrics_lenient = filter_units_by_quality(metrics, tier="lenient")

report_folder = base_folder / "metrics_tier_report"
report_folder.mkdir(parents=True, exist_ok=True)

# save unit IDs of each tier to CSV
(metrics_strict.index.to_series()
     .to_frame(name="unit_id")
     .to_csv(base_folder / "metrics_tier_report" / "units_strict.csv", index=False))
(metrics_medium.index.to_series()
     .to_frame(name="unit_id")
     .to_csv(base_folder / "metrics_tier_report" / "units_medium.csv", index=False))
(metrics_lenient.index.to_series()
     .to_frame(name="unit_id")
     .to_csv(base_folder / "metrics_tier_report" / "units_lenient.csv", index=False))

logging.info(f"Unit counts per tier: strict {len(metrics_strict)}, medium {len(metrics_medium)}, lenient {len(metrics_lenient)}")

# Create a folder to save lists & reports per tier
tier_report_folder = base_folder / "metrics_tier_report"
tier_report_folder.mkdir(parents=True, exist_ok=True)

for tier in ("strict", "medium", "lenient"):
    df_t = filter_units_by_quality(metrics, tier=tier)
    unit_ids = df_t.index.values
    logging.info(f"Tier '{tier}' — keeping {len(unit_ids)} units")
    # Save unit-IDs list
    df_t.index.to_series().to_frame(name="unit_id").to_csv(
        tier_report_folder / f"units_{tier}.csv", index=False
    )
    # Build cleaned analyzer for that tier
    out_analyzer_folder = base_folder / f"analyzer_{tier}"
    clean_an = analyzer.select_units(unit_ids,
                                     folder=out_analyzer_folder,
                                     format="binary_folder")
    logging.info(f"Saved cleaned analyzer for tier '{tier}' at {out_analyzer_folder}")
    # --- NEW CHECK: skip report if no units ---
    if len(unit_ids) == 0:
        logging.warning(f"Tier '{tier}' has 0 units — skipping export_report()")
        continue
    # Export report for cleaned analyzer
    report_out = base_folder / f"report_{tier}"
    si.export_report(clean_an, report_out, format="png")
    logging.info(f"Exported report for tier '{tier}' at {report_out}")

logging.info("Tiered QC + cleaned analyzers done for strict/medium/lenient sets.")


# In[8]:


analyzer_medium = si.load_sorting_analyzer(base_folder / 'analyzer_medium')
si.plot_sorting_summary(analyzer_medium, curation=True, backend='spikeinterface_gui')
logging.info("Sorting summary plot completed.")


# In[ ]:




