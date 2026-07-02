# NeuropixelsAnalysis

Python and MATLAB tools for analyzing Neuropixels high-density electrophysiology data, including single-unit spike analysis, LFP extraction, digital event decoding, and synchronization across recording systems.

## Overview

This repository contains analysis tools for Neuropixels probe recordings in awake-behaving mice. Neuropixels probes record from hundreds of channels simultaneously, enabling large-scale single-unit and LFP characterization across cortical and subcortical circuits. The pipeline covers post-spike-sorting analysis of sorted units, LFP extraction and processing, digital event decoding for behavioral synchronization, and timestamp alignment across NIDQ and Spike2 acquisition systems.

## Contents

- `npxAnalysis.py` — Python pipeline for post-sorting Neuropixels single-unit analysis; loads Kilosort/Phy-sorted spike data, computes per-unit firing rates across behavioral epochs (baseline, tone, trace, airpuff), generates ISI histograms and raster plots, and maps units to recording channels via CSV parsing
- `AlexLFPextract.ipynb` — Jupyter notebook for LFP band extraction and spectral analysis from Neuropixels AP and LFP band recordings
- `apAnalysis.m` — MATLAB analysis script for action potential band data; extracts spike timestamps and computes epoch-specific firing rate statistics
- `decodeEyeblinkDigMark.m` — Decodes digital marker channels to extract eye-blink conditioning event timestamps (CS onset, US onset, inter-trial intervals) from Spike2 recordings
- `extractDigMark32_CEDS64.m` — Extracts 32-bit digital marker values from CED Spike2 files using the CEDS64 library
- `extractNidqStartStop.m` — Extracts NIDQ acquisition start and stop timestamps for aligning Neuropixels recording epochs to behavioral events
- `npxTimestamps.m` — Timestamp conversion and alignment between Neuropixels sample clocks and behavioral event markers

## Requirements

- Python 3.7+ (numpy, scipy, pandas, matplotlib, pathlib)
- Jupyter Notebook
- MATLAB R2018a or later
- CED CEDS64 library (for Spike2 file access)
- Kilosort/Phy-sorted data as input

## Usage

1. Run spike sorting on raw Neuropixels AP band data using Kilosort and manually curate in Phy
2. Use `npxAnalysis.py` to load sorted units and compute behavioral epoch firing rates
3. Use `decodeEyeblinkDigMark.m` and `extractNidqStartStop.m` to extract and align behavioral event timestamps
4. Use `AlexLFPextract.ipynb` for LFP band spectral analysis
