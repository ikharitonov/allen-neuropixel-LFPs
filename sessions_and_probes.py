import numpy as np
import pandas as pd
from pathlib import Path
from allensdk.brain_observatory.ecephys.ecephys_project_cache import EcephysProjectCache
import pickle

output_dir = Path.home() / 'Desktop' / 'disk1' / 'ecephys_data'

cache = EcephysProjectCache.from_warehouse(manifest=output_dir / 'manifest.json')
probe_df = cache.get_probes()

output_list = []
for session_id in probe_df.ecephys_session_id.unique():
    ind_list = probe_df[(probe_df.ecephys_session_id==session_id) & (probe_df['ecephys_structure_acronyms'].apply(lambda x: 'VISp' in x or 'VISpm' in x)) & (probe_df['has_lfp_data'] == True)].index.to_list()
    if len(ind_list) == 2:
        ind_list = [session_id] + ind_list
        output_list.append(ind_list)

with open('VISpm_VISp_probes.pkl', 'wb') as f:
    pickle.dump(output_list, f)