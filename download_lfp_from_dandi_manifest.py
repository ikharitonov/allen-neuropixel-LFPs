import numpy as np
import json
import pickle
import subprocess
from pathlib import Path

"""
Using dandi_manifest.json, created manually from https://api.dandiarchive.org/api/dandisets/000021/versions/draft/assets/ , and VISpm_VISp_probes.pkl (from session_and_probes.py), this file automatically downloads session and LFP NWB files using Dandi.
Dataset link: https://dandiarchive.org/dandiset/000021/

Download folder is specified by OUTPUT_DIR. If a file has already been downloaded it will be skipped.
"""

OUTPUT_DIR = Path.home()/'Desktop'/'dandi_lfp'

manifest_path = 'dandi_manifest.json'
with open(manifest_path, 'r') as f: data = json.load(f)
with open('VISpm_VISp_probes.pkl', 'rb') as f:
    probe_list = pickle.load(f)

session_list = np.array(probe_list)[:,0].flatten()
probe_list = np.array(probe_list)[:,1:].flatten()

data = data['results']

for f in data:
    subject_name, probe_file = f['path'].split('/')
    if 'probe' in probe_file:
        probe_id = int(probe_file.split('_')[2].split('-')[1])
        # dandi download the probe files from the list
        if probe_id in probe_list:
            asset_id = f['asset_id']
            result = subprocess.run(['dandi','download','-o',OUTPUT_DIR,'-e','skip',f'https://api.dandiarchive.org/api/dandisets/000021/versions/draft/assets/{asset_id}/download/'], capture_output=True, text=True)
            print(result.stdout)
            print(result.stderr)
    else:
        session_id = int(probe_file.split('-')[-1].split('.')[0])
        # Download the session file
        if session_id in session_list:
            asset_id = f['asset_id']
            result = subprocess.run(['dandi','download','-o',OUTPUT_DIR,'-e','skip',f'https://api.dandiarchive.org/api/dandisets/000021/versions/draft/assets/{asset_id}/download/'], capture_output=True, text=True)
            print(result.stdout)
            print(result.stderr)