import numpy as np
import pandas as pd
import pickle
import time
from pathlib import Path
from allensdk.brain_observatory.ecephys.ecephys_project_cache import EcephysProjectCache


"""

This script is for computing 8 averaged LFP traces for toWhite/toBlack flash presentation conditions and running/nonrunning states.
It requires velocity thresholds (chosen manually or automatically in 1-velocity_threshold.ipynb) provided in chosen_velocity_thresholds.pkl.
As well as paths to directories with all relevant ecephys session/LFP data and reference of session IDs and VISpm/VISp probe IDs (VISpm_VISp_probes.pkl).

1) Available session IDs and corresponding velocity thresholds are read from a file.
2) Cache is initialised and probe IDs are read.
3) Session data is loaded and probe IDs are assigned with corresponding areas.
4) LFP files are loaded into memory.
5) LFP traces are quality controlled, clipped to the period of flash stimuli presentations.
6) Channels with maximum amplitude LFP traces are selected within VISpm/VISp.
7) Velocity time series of the mouse are loaded from session data.
8) Each flash presentation is assigned with a class (toWhite OR toBlack AND running OR nonrunning) based on the flash data itself and chosen velocity threshold.
9) Electrode depth of the maximum amplitude VISpm/VISp channel is read from metadata.
10) VISpm/VISp LFP traces are aligned to each class of flash presentation according to chosen time window range (before and after the start time of the flash).
11) Aligned LFP traces are averaged over presentations and saved into .npy files.

"""


def align_LFPs(stim_presentations, lfp_slice, channel_to_select, window_range=[-0.5,0.5], sampling_frequency=1250):
    presentation_times = stim_presentations.start_time.values
    presentation_ids = stim_presentations.index.values
    
    window_length = window_range[1]-window_range[0] # in seconds
    # trial_window = np.arange(-0.5, 3.5, 1/1250) # np.arange(start, stop, step)
    trial_window = np.linspace(window_range[0], window_range[1], int(window_length*sampling_frequency))
    # Trial window array (including 0.5s before and 0.5s after stimulus) is added to each timepoint linked to flash stimulus presentation
    # Then all of these stim-presentation-related trial windows are concatenated into a single array
    time_selection = np.concatenate([trial_window + t for t in presentation_times])
    
    # Important to note that 1/500 increment here corresponds to downsampling to 500Hz
    
    # This repeatedly assigns basic trial window values (-0.5s:0.5s) to every stimulus presentation id, under the column label of 'time_from_presentation_onset'
    inds = pd.MultiIndex.from_product((presentation_ids, trial_window), 
                                      names=('presentation_id', 'time_from_presentation_onset'))
    
    # Selecting LFP data at trial window constrained flash presentation times by the nearest value available
    # Saving this into a newly created 'aligned_lfp' inside of xarray dataset
    dataset = lfp_slice.sel(channel=channel_to_select).sel(time = time_selection, method='nearest').to_dataset(name = 'aligned_lfp')
    # Assigning new variabes ('presentation_id' and 'time_from_presentation_onset' from inds) to the xarray dataset
    # Removing time variable
    # Now, for each presentation_id, there is a corresponding time_from_presentation_onset-long array of LFP data for all the channels
    dataset = dataset.assign(time=inds).unstack('time')
    return dataset['aligned_lfp']



def get_every_area_channel_id_dict(merged_df):
    probe_every_area_channel_ids = {}
    for area in merged_df.ecephys_structure_acronym.unique():
        if pd.isnull(area): # to select NaN area key
            probe_every_area_channel_ids['NaN'] = merged_df[pd.isna(merged_df.ecephys_structure_acronym)].index.to_list()  # works for selecting channels from the NaN area
        else:
            probe_every_area_channel_ids[area] = merged_df[merged_df.ecephys_structure_acronym==area].index.to_list()
    return probe_every_area_channel_ids



def get_max_channel_id(data):
    max_amplitudes = np.max(np.abs(data), axis=0)
    return int(np.argmax(max_amplitudes.data))



def get_max_channel_id_list(data):
    max_amplitudes = np.max(np.abs(data), axis=0)
    return np.argsort(-max_amplitudes).data



def load_and_process_lfps(session, selected_VISpm_probe_id, selected_VISp_probe_id, channel_df):
    t1 = time.time()
    VISpm_lfp = session.get_lfp(selected_VISpm_probe_id)
    print(f'LFP data loaded in {time.time()-t1} seconds.')
    t1 = time.time()
    VISp_lfp = session.get_lfp(selected_VISp_probe_id)
    print(f'LFP data loaded in {time.time()-t1} seconds.')
    
    stim_epochs = session.get_stimulus_epochs()
    flash_index = stim_epochs[stim_epochs.stimulus_name=='flashes'].index
    print(f'{len(flash_index)} flash block/s present: choosing the first with duration = {stim_epochs[stim_epochs.index==flash_index[0]].duration.item()} seconds.')
    flash_index = flash_index[0] # [0] instead of using ".item()" to only grab the first in case there are multiple flash blocks
    global_start_time = stim_epochs[stim_epochs.index==flash_index-1].start_time.item()
    global_end_time = stim_epochs[stim_epochs.index==flash_index+1].stop_time.item()
    
    # Quality check for channels that are in the LFP NWB file but not in the Allen cache
    for ch in VISpm_lfp.channel.values:
        if ch not in channel_df.index.to_list():
            VISpm_lfp = VISpm_lfp.drop_sel(channel=ch)
    for ch in VISp_lfp.channel.values:
        if ch not in channel_df.index.to_list():
            VISp_lfp = VISp_lfp.drop_sel(channel=ch)
    
    # Organise all channels into a dictionary with keys as area names and values as lists with related LFP recording channels
    VISpm_cache_channel_df_merged = channel_df.loc[VISpm_lfp.channel.values]
    VISp_cache_channel_df_merged = channel_df.loc[VISp_lfp.channel.values]

    VISpm_probe_every_area_channel_ids = get_every_area_channel_id_dict(VISpm_cache_channel_df_merged)
    VISp_probe_every_area_channel_ids = get_every_area_channel_id_dict(VISp_cache_channel_df_merged)

    # VISpm_channel_ids = [ch for ch in VISpm_lfp.channel.values if channel_df[channel_df.index==ch]['ecephys_structure_acronym'].item()=='VISpm']
    # VISp_channel_ids = [ch for ch in VISp_lfp.channel.values if channel_df[channel_df.index==ch]['ecephys_structure_acronym'].item()=='VISp']

    VISpm_channel_ids = VISpm_probe_every_area_channel_ids['VISpm']
    VISp_channel_ids = VISp_probe_every_area_channel_ids['VISp']

    print(f'VISpm channels = {len(VISpm_channel_ids)}')
    print(f'VISp channels = {len(VISp_channel_ids)}')
    
    VISpm_lfp_slice = VISpm_lfp.sel(time=slice(global_start_time, global_end_time))
    VISp_lfp_slice = VISp_lfp.sel(time=slice(global_start_time, global_end_time))
    
    
    # For each area in the two probes, determine max amplitude channel/s
    VISpm_probe_every_area_max_channel = {}
    VISpm_probe_every_area_max_channel_list = {}
    VISp_probe_every_area_max_channel = {}
    VISp_probe_every_area_max_channel_list = {}

    for area, ch_ids in VISpm_probe_every_area_channel_ids.items():
        VISpm_probe_every_area_max_channel[area] = ch_ids[get_max_channel_id(VISpm_lfp_slice.sel(channel=slice(ch_ids[0], ch_ids[-1])))]
        VISpm_probe_every_area_max_channel_list[area] = get_max_channel_id_list(VISpm_lfp_slice.sel(channel=slice(ch_ids[0], ch_ids[-1])))
    for area, ch_ids in VISp_probe_every_area_channel_ids.items():
        VISp_probe_every_area_max_channel[area] = ch_ids[get_max_channel_id(VISp_lfp_slice.sel(channel=slice(ch_ids[0], ch_ids[-1])))]
        VISp_probe_every_area_max_channel_list[area] = get_max_channel_id_list(VISp_lfp_slice.sel(channel=slice(ch_ids[0], ch_ids[-1])))

    VISpm_max_channel = VISpm_probe_every_area_max_channel['VISpm']
    VISpm_max_channel_list = VISpm_probe_every_area_max_channel_list['VISpm']
    VISp_max_channel = VISp_probe_every_area_max_channel['VISp']
    VISp_max_channel_list = VISp_probe_every_area_max_channel_list['VISp']

    print(f'VISpm_max_channel = {VISpm_max_channel}')
    print(f'VISp_max_channel = {VISp_max_channel}')
    
    return (VISpm_lfp_slice, VISp_lfp_slice, VISpm_max_channel, VISp_max_channel)



def get_eight_flash_conditions(session, VELOCITY_THRESHOLD):
    
    running_speed_midpoints = session.running_speed["start_time"] + (session.running_speed["end_time"] - session.running_speed["start_time"]) / 2
    
    time_velocity_df = pd.concat([running_speed_midpoints, session.running_speed["velocity"]],axis=1,keys=['running_speed_midpoint','velocity'])
    
    stim_flashes = session.get_stimulus_table(['flashes'])
    
    # Classifying each flash as running or non-running according to associated running velocity and VELOCITY_THRESHOLD computed above

    flash_cat_inds = []
    flash_cat_isRunning = []

    for stim_id,stim_pres in stim_flashes.iterrows():
        # From running velocity dataframe, select all recorded velocity values during presentation of a single stimulus (between start and end of a flash)
        t1 = time_velocity_df[time_velocity_df.running_speed_midpoint < stim_pres.start_time].index[-1] # getting the closest timepoint index just before the start of the flash as the lower bound
        t2 = time_velocity_df[time_velocity_df.running_speed_midpoint > stim_pres.stop_time].index[0] # getting the closest timepoint index just before the start of the flash as the lower bound
        filtered_df = time_velocity_df[(time_velocity_df.running_speed_midpoint >= stim_pres.start_time) & (time_velocity_df.running_speed_midpoint <= stim_pres.stop_time)]
        # Calculate mean velocity during this presentation
        vel_mean = filtered_df.velocity.mean()
        # Assign boolean classification and store in a Series
        isRunning = False
        if vel_mean >= VELOCITY_THRESHOLD: isRunning = True
        flash_cat_inds.append(stim_id)
        flash_cat_isRunning.append(isRunning)
    flash_cat_series = pd.Series(flash_cat_isRunning, index=flash_cat_inds)
    
    
    stim_flashes['isRunning'] = flash_cat_series
    
    
    print(f'{stim_flashes[stim_flashes.isRunning==True].shape[0]} presentations where animal is running')
    print(f'{stim_flashes[stim_flashes.isRunning==False].shape[0]} presentations where animal is not running')
    
    
    toWhite_running_flashes = stim_flashes[(stim_flashes.stimulus_condition_id==245) & (stim_flashes.isRunning)]
    print('toWhite_running_flashes', toWhite_running_flashes.shape[0], 'presentations')

    toBlack_running_flashes = stim_flashes[(stim_flashes.stimulus_condition_id==244) & (stim_flashes.isRunning)]
    print('toBlack_running_flashes', toBlack_running_flashes.shape[0], 'presentations')

    toWhite_nonrunning_flashes = stim_flashes[(stim_flashes.stimulus_condition_id==245) & (stim_flashes.isRunning==False)]
    print('toWhite_nonrunning_flashes', toWhite_nonrunning_flashes.shape[0], 'presentations')

    toBlack_nonrunning_flashes = stim_flashes[(stim_flashes.stimulus_condition_id==244) & (stim_flashes.isRunning==False)]
    print('toBlack_nonrunning_flashes', toBlack_nonrunning_flashes.shape[0], 'presentations')

    return (toWhite_running_flashes, toBlack_running_flashes, toWhite_nonrunning_flashes, toBlack_nonrunning_flashes)



def run(cache, probe_list, session_id, VELOCITY_THRESHOLD, window_range, sf, output_folder):

    channel_df = cache.get_channels()
    probes_df = cache.get_probes()
    
    probe_dict = {e[0]: [e[1], e[2]] for e in probe_list}    
    
    # Loading the selected session and determining which probe is VISpm/VISp
    selected_session_list = probe_dict[session_id]

    session = cache.get_session_data(session_id)

    if 'VISpm' in probes_df[probes_df.index==selected_session_list[0]].ecephys_structure_acronyms.item():
        selected_VISpm_probe_id = selected_session_list[0]
        selected_VISp_probe_id = selected_session_list[1]
    else:
        selected_VISpm_probe_id = selected_session_list[1]
        selected_VISp_probe_id = selected_session_list[0]
    
    VISpm_lfp_slice, VISp_lfp_slice, VISpm_max_channel, VISp_max_channel = load_and_process_lfps(session, selected_VISpm_probe_id, selected_VISp_probe_id, channel_df)
    
    
    toWhite_running_flashes, toBlack_running_flashes, toWhite_nonrunning_flashes, toBlack_nonrunning_flashes = get_eight_flash_conditions(session, VELOCITY_THRESHOLD)
    
    try:
        VISpm_electrode_depth = int(channel_df[channel_df.index==VISpm_max_channel].dorsal_ventral_ccf_coordinate.item()) # in microns
        VISp_electrode_depth = int(channel_df[channel_df.index==VISp_max_channel].dorsal_ventral_ccf_coordinate.item())
    except:
        VISpm_electrode_depth = 'NaN'
        VISp_electrode_depth = 'NaN'


    # VISpm

    aligned_lfp = align_LFPs(toWhite_running_flashes, VISpm_lfp_slice, VISpm_max_channel, window_range=window_range, sampling_frequency=sf)
    filename = f'sessionID_{session_id}_area_VISpm_condition_toWhite_running_True_flashesAveragedOver_{toWhite_running_flashes.shape[0]}_micronsElectrodeDepth_{VISpm_electrode_depth}.npy'
    np.save(output_folder/filename, [aligned_lfp.time_from_presentation_onset, aligned_lfp.mean(dim='presentation_id')])

    aligned_lfp = align_LFPs(toBlack_running_flashes, VISpm_lfp_slice, VISpm_max_channel, window_range=window_range, sampling_frequency=sf)
    filename = f'sessionID_{session_id}_area_VISpm_condition_toBlack_running_True_flashesAveragedOver_{toBlack_running_flashes.shape[0]}_micronsElectrodeDepth_{VISpm_electrode_depth}.npy'
    np.save(output_folder/filename, [aligned_lfp.time_from_presentation_onset, aligned_lfp.mean(dim='presentation_id')])

    aligned_lfp = align_LFPs(toWhite_nonrunning_flashes, VISpm_lfp_slice, VISpm_max_channel, window_range=window_range, sampling_frequency=sf)
    filename = f'sessionID_{session_id}_area_VISpm_condition_toWhite_running_False_flashesAveragedOver_{toWhite_nonrunning_flashes.shape[0]}_micronsElectrodeDepth_{VISpm_electrode_depth}.npy'
    np.save(output_folder/filename, [aligned_lfp.time_from_presentation_onset, aligned_lfp.mean(dim='presentation_id')])

    aligned_lfp = align_LFPs(toBlack_nonrunning_flashes, VISpm_lfp_slice, VISpm_max_channel, window_range=window_range, sampling_frequency=sf)
    filename = f'sessionID_{session_id}_area_VISpm_condition_toBlack_running_False_flashesAveragedOver_{toBlack_nonrunning_flashes.shape[0]}_micronsElectrodeDepth_{VISpm_electrode_depth}.npy'
    np.save(output_folder/filename, [aligned_lfp.time_from_presentation_onset, aligned_lfp.mean(dim='presentation_id')])

    # VISp

    aligned_lfp = align_LFPs(toWhite_running_flashes, VISp_lfp_slice, VISp_max_channel, window_range=window_range, sampling_frequency=sf)
    filename = f'sessionID_{session_id}_area_VISp_condition_toWhite_running_True_flashesAveragedOver_{toWhite_running_flashes.shape[0]}_micronsElectrodeDepth_{VISp_electrode_depth}.npy'
    np.save(output_folder/filename, [aligned_lfp.time_from_presentation_onset, aligned_lfp.mean(dim='presentation_id')])

    aligned_lfp = align_LFPs(toBlack_running_flashes, VISp_lfp_slice, VISp_max_channel, window_range=window_range, sampling_frequency=sf)
    filename = f'sessionID_{session_id}_area_VISp_condition_toBlack_running_True_flashesAveragedOver_{toBlack_running_flashes.shape[0]}_micronsElectrodeDepth_{VISp_electrode_depth}.npy'
    np.save(output_folder/filename, [aligned_lfp.time_from_presentation_onset, aligned_lfp.mean(dim='presentation_id')])

    aligned_lfp = align_LFPs(toWhite_nonrunning_flashes, VISp_lfp_slice, VISp_max_channel, window_range=window_range, sampling_frequency=sf)
    filename = f'sessionID_{session_id}_area_VISp_condition_toWhite_running_False_flashesAveragedOver_{toWhite_nonrunning_flashes.shape[0]}_micronsElectrodeDepth_{VISp_electrode_depth}.npy'
    np.save(output_folder/filename, [aligned_lfp.time_from_presentation_onset, aligned_lfp.mean(dim='presentation_id')])

    aligned_lfp = align_LFPs(toBlack_nonrunning_flashes, VISp_lfp_slice, VISp_max_channel, window_range=window_range, sampling_frequency=sf)
    filename = f'sessionID_{session_id}_area_VISp_condition_toBlack_running_False_flashesAveragedOver_{toBlack_nonrunning_flashes.shape[0]}_micronsElectrodeDepth_{VISp_electrode_depth}.npy'
    np.save(output_folder/filename, [aligned_lfp.time_from_presentation_onset, aligned_lfp.mean(dim='presentation_id')])

    print(f'Session {session_id}: Averaged LFP traces saved for 8 flash/running conditions.\n')

if __name__ == '__main__':
    output_folder = Path.home() / 'Desktop' / 'disk2' / 'grand_average_lfps'
    
    # Open file with velocity thresholds chosen in 1-velocity_threshold.ipynb
    with open('chosen_velocity_thresholds.pkl', 'rb') as f:
        chosen_velocity_thresholds = pickle.load(f)
      
    # Initialise cache
    output_dir = Path.home() / 'Desktop' / 'disk2' / 'ecephys_data'
    cache = EcephysProjectCache.from_warehouse(manifest=output_dir / 'manifest.json')
    
    with open('VISpm_VISp_probes.pkl', 'rb') as f:
        probe_list = pickle.load(f)
    
    window = [-1, 2.25] # averaging time window with respect to start_time of flash presentation
    sf = 1250 # LFP sampling frequency
     
    # run(cache, probe_list, 755434585, 7, window, sf, output_folder)
    for session_id, thresh in chosen_velocity_thresholds.items():
        try: 
            run(cache, probe_list, session_id, thresh, window, sf, output_folder)
        except:
            print(f'ERROR: SESSION ID = {session_id}.\n')