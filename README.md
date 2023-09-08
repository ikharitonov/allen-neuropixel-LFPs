# Analysis of LFP data from Allen Visual Coding Neuropixels dataset

__tutorial.ipynb__ : following the tutorial https://allensdk.readthedocs.io/en/latest/_static/examples/nb/ecephys_session.html

__session_715093703.ipynb__ : exploring LFP data, current source density (CSD), flash stimuli presentations, spectrograms, power spectral density plots in one session

__sessions_and_probes.py__ : selecting session and probe IDs which correspond to VISpm and VISp areas from the database

__download_lfp_from_dandi_manifest.py__ : downloading NWB files with LFP data using dandi https://dandiarchive.org/dandiset/000021/draft

__rename_nwb_files.py__ : moving and renaming donwloaded NWB files according to directory format of AllenSDK (EcephysProjectCache)

__spectral_analysis_pipepline.ipynb__ : pipeline for accessing LFP data from NWB files, selecting maximum signal channels, spectrograms, PSD, running speed analysis (initial), stimulus-averaged LFP trace visualisation (initial), testing normality of data, propagating uncertainty and calculating mean power spectrum

__generalised_phase.ipynb__ : initial implementation of https://github.com/mullerlab/generalized-phase in Python

__lfp_grand_average/1-velocity_threshold.ipynb__ : fitting a bimodal curve to running velocity distribution in each session to determine threshold

__lfp_grand_average/2-stimulus_averaged_lfp.py__ : saving LFP traces aligned to flash stimuli presentations over 8 conditions (toWhite/toBlack and running/nonrunning) 

__lfp_grand_average/2-stimulus_averaged_power_spectrum.py__ : saving power spectra (Welch's method) calculated from LFP traces aligned to flash stimuli presentations.

__lfp_grand_average/parse_filenames.py__ : returns a Pandas DataFrame containing experimental condition metadata after parsing specified directory.

__lfp_grand_average/3-grand_average.ipynb__ : calculates median + IQR for each session, weighted average across sessions

__lfp_grand_average/3-power_spectrum_grand_average.ipynb__ : same approach as in 3-grand_average.ipynb applied to power spectra.