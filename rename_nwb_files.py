import os
from pathlib import Path
import shutil
from allensdk.brain_observatory.ecephys.ecephys_project_cache import EcephysProjectCache

"""
This script renames NWB files downloaded with download_lfp_from_dandi_manifest.py according to AllenSDK nomenclature and user specified original_path and resulting_path.
"""

def mkdir(path):
    if not os.path.exists(path): os.makedirs(path)
    return path

def get_session_id_from_probe_file(file):
    for test_probe in probes_df.index.unique():
        if str(test_probe) in file:
            return test_probe, probes_df[probes_df.index==test_probe].ecephys_session_id.item()

def get_session_id_from_session_file(file):
    for test_session in probes_df.ecephys_session_id.unique():
        if str(test_session) in file:
            return test_session
        
original_path = Path.home()/'Desktop'/'disk2'/'dandi_lfp'
resulting_path = Path.home()/'Desktop'/'disk2'/'ecephys_data'


cache = EcephysProjectCache.from_warehouse(manifest = Path.home() / 'Desktop' / 'disk2' / 'ecephys_data' / 'manifest.json')
probes_df = cache.get_probes()


files = os.listdir(original_path)

for f in files:
    if 'probe' not in f: # Session file
        session_id = get_session_id_from_session_file(f)
        shutil.move(original_path/f, mkdir(resulting_path/f'session_{session_id}')/f'session_{session_id}.nwb')
    else: # Probe file
        probe_id, session_id = get_session_id_from_probe_file(f)
        shutil.move(original_path/f, mkdir(resulting_path/f'session_{session_id}')/f'probe_{probe_id}_lfp.nwb')
print('All files moved:')
print('From',original_path)
print('To',resulting_path)