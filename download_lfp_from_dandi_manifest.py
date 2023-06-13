import numpy as np
import json
import pickle
import subprocess

manifest_path = 'dandi_manifest.json'
with open(manifest_path, 'r') as f: data = json.load(f)
with open('VISpm_VISp_probes.pkl', 'rb') as f:
    probe_list = pickle.load(f)

probe_list = np.array(probe_list)[:,1:].flatten()

data = data['results']

for f in data:
    session_name, probe_file = f['path'].split('/')
    if 'probe' in probe_file:
        probe_id = int(probe_file.split('_')[2].split('-')[1])
        # dandi download the session file + the probe files from the list
        if probe_id in probe_list:
            asset_id = f['asset_id']
            result = subprocess.run(['dandi','download',f'https://api.dandiarchive.org/api/dandisets/000021/versions/draft/assets/{asset_id}/download/'], capture_output=True, text=True)
            # result = subprocess.run(['dandi', '--version'], capture_output=True, text=True)
            print(result.stdout)
            print(result.stderr)