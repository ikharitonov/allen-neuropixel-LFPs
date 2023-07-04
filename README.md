# Analysis of LFP data from Allen Visual Coding Neuropixels dataset

__tutorial.ipynb__ : following the tutorial https://allensdk.readthedocs.io/en/latest/_static/examples/nb/ecephys_session.html

__session_715093703.ipynb__ : exploring LFP data, current source density (CSD), flash stimuli presentations, spectrograms, power spectral density plots in one session

__sessions_and_probes.py__ : selecting session and probe IDs which correspond to VISpm and VISp areas from the database

__download_lfp_from_dandi_manifest.py__ : downloading NWB files with LFP data using dandi https://dandiarchive.org/dandiset/000021/draft

__rename_nwb_files.py__ : moving and renaming donwloaded NWB files according to directory format of AllenSDK (EcephysProjectCache)

__spectral_analysis_pipepline.ipynb__ : main pipeline for accessing LFP data from NWB files, selecting maximum signal channels, spectrograms, PSD, running speed analysis (initial) and stimulus-averaged LFP trace visualisation (initial)

__generalised_phase.ipynb__ : initial implementation of https://github.com/mullerlab/generalized-phase in Python
