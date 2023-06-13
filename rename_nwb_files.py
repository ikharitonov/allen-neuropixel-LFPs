import os
from pathlib import Path
import shutil

"""
This script renames NWB files downloaded with download_lfp_from_dandi_manifest.py according to AllenSDK nomenclature and user specified original_path and resulting_path.
"""

def mkdir(path):
    if not os.path.exists(path): os.makedirs(path)
    return path

original_path = Path.home()/'Desktop'/'disk1'/'dandi_lfp'
resulting_path = Path.home()/'Desktop'/'disk1'/'ecephys_data'

files = os.listdir(original_path)

for f in files:
    if 'probe' not in f: # Session file
        session_id = int(f.split('_')[1].split('-')[1].split('.')[0])
        shutil.move(original_path/f, mkdir(resulting_path/f'session_{session_id}')/f'session_{session_id}.nwb')
    else: # Probe file
        session_id = int(f.split('_')[1].split('-')[1])
        probe_id = int(f.split('_')[2].split('-')[1])
        shutil.move(original_path/f, mkdir(resulting_path/f'session_{session_id}')/f'probe_{probe_id}_lfp.nwb')
print('All files moved:')
print('From',original_path)
print('To',resulting_path)