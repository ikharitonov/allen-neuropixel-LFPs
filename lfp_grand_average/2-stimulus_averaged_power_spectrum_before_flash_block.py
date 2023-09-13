import numpy as np
import pandas as pd
import pickle
import time
from pathlib import Path
from scipy import signal
from allensdk.brain_observatory.ecephys.ecephys_project_cache import EcephysProjectCache


"""

This script is similar to 2-stimulus_averaged_power_spectrum.py, except that it separates spontaneous block (no visual stimuli presented) preceding flashes into running / non-running segments and calculates power spectrum for those.
Velocity is low-passed at 5Hz to smoothen it.
This results into 4 conditions (VISpm_running, VISpm_nonrunning, VISp_running and VISp_nonrunning) instead of 8.

"""



def segment_LFPs_compute_ps(time_segments, lfp_slice, channel_to_select, fs=1250):
    
    desired_frequency_resolution = 1 # hardcoded 1Hz frequency resolution

    lfp_slice = lfp_slice.sel(channel=channel_to_select)

    i = 0
    num_segments_collected = 0
    ps_collected = []

    for segment in time_segments:
        # Only select segments which have >= 1250 samples to get at least 1Hz frequency resolution in Welch's method
        data_selected = lfp_slice.sel(time=slice(segment[0], segment[-1])).data
        if data_selected.shape[0] >= fs // desired_frequency_resolution:
            num_segments_collected += 1
            f,s = signal.welch(data_selected, fs, scaling='spectrum', nperseg=fs//desired_frequency_resolution)
            ps_collected.append(s)
        i+=1
    if len(ps_collected) == 0: return [], [], 0
    else: return f, ps_collected, num_segments_collected



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



def separate_segments(arr):
    if type(arr)==tuple: arr = arr[0]
    
    segments = []
    current_segment = []

    for num in arr:
        if not current_segment or num == current_segment[-1] + 1:
            current_segment.append(num)
        else:
            segments.append(np.array(current_segment))
            current_segment = [num]

    if current_segment:
        segments.append(np.array(current_segment))

    return segments



def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = signal.butter(order, cutoff, fs=fs, btype='low', analog=False)
    y = signal.lfilter(b, a, data)
    return y



def get_two_running_conditions(session, VELOCITY_THRESHOLD):

    # This function retrieves velocity time series (dataframe) of the spontaneous block, separates it into running and non-running time segments
    # Returns a tuple of (running_time_segments, nonrunning_time_segments)

    running_speed_midpoints = session.running_speed["start_time"] + (session.running_speed["end_time"] - session.running_speed["start_time"]) / 2
    time_velocity_df = pd.concat([running_speed_midpoints, session.running_speed["velocity"]],axis=1,keys=['running_speed_midpoint','velocity'])

    stim_epochs = session.get_stimulus_epochs()
    flash_index = stim_epochs[stim_epochs.stimulus_name=='flashes'].index[0]
    spontaneous_start_time = stim_epochs[stim_epochs.index==flash_index-1].start_time.item()
    sponteneuous_stop_time = stim_epochs[stim_epochs.index==flash_index-1].stop_time.item()

    inds = np.where((time_velocity_df.running_speed_midpoint.to_numpy()>spontaneous_start_time) & (time_velocity_df.running_speed_midpoint.to_numpy()<sponteneuous_stop_time))
    spontaneous_t = time_velocity_df.running_speed_midpoint.to_numpy()[inds]
    spontaneous_v = time_velocity_df.velocity.to_numpy()[inds]

    # Filtering at 5Hz to smoothen velocity
    spontaneous_v = butter_lowpass_filter(spontaneous_v, 5, 1250)

    running_inds = np.where(spontaneous_v >= VELOCITY_THRESHOLD)
    nonrunning_inds = np.where(spontaneous_v < VELOCITY_THRESHOLD)

    running_segments = separate_segments(running_inds)
    nonrunning_segments = separate_segments(nonrunning_inds)

    running_time_segments = []
    nonrunning_time_segments = []
    for segment in running_segments:
        segment = (segment,)
        running_time_segments.append(spontaneous_t[segment])
    for segment in nonrunning_segments:
        segment = (segment,)
        nonrunning_time_segments.append(spontaneous_t[segment])

    return (running_time_segments, nonrunning_time_segments)



def run(cache, probe_list, session_id, VELOCITY_THRESHOLD, fs, output_folder):

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
    
    
    # toWhite_running_flashes, toBlack_running_flashes, toWhite_nonrunning_flashes, toBlack_nonrunning_flashes = get_eight_flash_conditions(session, VELOCITY_THRESHOLD)
    running_time_segments, nonrunning_time_segments = get_two_running_conditions(session, VELOCITY_THRESHOLD)

    try:
        VISpm_electrode_depth = int(channel_df[channel_df.index==VISpm_max_channel].dorsal_ventral_ccf_coordinate.item()) # in microns
        VISp_electrode_depth = int(channel_df[channel_df.index==VISp_max_channel].dorsal_ventral_ccf_coordinate.item())
    except:
        VISpm_electrode_depth = 'NaN'
        VISp_electrode_depth = 'NaN'


    # VISpm

    f, s_collected, num_segments = segment_LFPs_compute_ps(running_time_segments, VISpm_lfp_slice, VISpm_max_channel)
    print(f'num_segments collected = {num_segments}')
    filename = f'sessionID_{session_id}_area_VISpm_running_True_segmentsCollected_{num_segments}_micronsElectrodeDepth_{VISpm_electrode_depth}.npy'
    np.save(output_folder/filename, [f,s_collected,num_segments])

    f, s_collected, num_segments = segment_LFPs_compute_ps(nonrunning_time_segments, VISpm_lfp_slice, VISpm_max_channel)
    print(f'num_segments collected = {num_segments}')
    filename = f'sessionID_{session_id}_area_VISpm_running_False_segmentsCollected_{num_segments}_micronsElectrodeDepth_{VISpm_electrode_depth}.npy'
    np.save(output_folder/filename, [f,s_collected,num_segments])

    # VISp

    f, s_collected, num_segments = segment_LFPs_compute_ps(running_time_segments, VISp_lfp_slice, VISp_max_channel)
    print(f'num_segments collected = {num_segments}')
    filename = f'sessionID_{session_id}_area_VISp_running_True_segmentsCollected_{num_segments}_micronsElectrodeDepth_{VISp_electrode_depth}.npy'
    np.save(output_folder/filename, [f,s_collected,num_segments])

    f, s_collected, num_segments = segment_LFPs_compute_ps(nonrunning_time_segments, VISp_lfp_slice, VISp_max_channel)
    print(f'num_segments collected = {num_segments}')
    filename = f'sessionID_{session_id}_area_VISp_running_False_segmentsCollected_{num_segments}_micronsElectrodeDepth_{VISp_electrode_depth}.npy'
    np.save(output_folder/filename, [f,s_collected,num_segments])

    print(f'Session {session_id}: Power spectra of LFP traces saved for 4 areas/running conditions.\n')

if __name__ == '__main__':
    output_folder = Path.home() / 'Desktop' / 'disk2' / 'lfp_power_spectra_dump_before_flashes'
    
    # Open file with velocity thresholds chosen in 1-velocity_threshold.ipynb
    with open('chosen_velocity_thresholds.pkl', 'rb') as f:
        chosen_velocity_thresholds = pickle.load(f)
      
    # Initialise cache
    output_dir = Path.home() / 'Desktop' / 'disk2' / 'ecephys_data'
    cache = EcephysProjectCache.from_warehouse(manifest=output_dir / 'manifest.json')
    
    with open('VISpm_VISp_probes.pkl', 'rb') as f:
        probe_list = pickle.load(f)
    
    fs = 1250 # LFP sampling frequency
     
    # run(cache, probe_list, 755434585, 7, window, sf, output_folder)
    for session_id, thresh in chosen_velocity_thresholds.items():
        try:
            run(cache, probe_list, session_id, thresh, fs, output_folder)
        except Exception as e:
            print(f'ERROR: SESSION ID = {session_id}.')
            print(e)
            print('\n')