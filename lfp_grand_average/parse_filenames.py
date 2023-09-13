import os
import numpy as np
import pandas as pd

def tryconvert(value, default, *types):
        for t in types:
            try:
                return t(value)
            except (ValueError, TypeError):
                continue
        return default

def parse_filenames(path, conversion_dict=None):
    file_list = os.listdir(path)
    
    get_splits = lambda l: [x.split('_') for x in l] # split strings
    get_k_v = lambda l: {l[i]:l[i+1] for i in range(0, len(l), 2)} # create dictionaries from successive items
    
    parsed_dict = [ get_k_v(l) for l in get_splits(file_list) ]
    
    for i in range(len(parsed_dict)): parsed_dict[i]['filename'] = file_list[i]
    
    if not conversion_dict:
        conversion_dict = {'sessionID': lambda x: int(x), 'area': lambda x: x, 'condition': lambda x: x, 'running': lambda x: x, 'flashesAveragedOver': lambda x: int(x), 'micronsElectrodeDepth': lambda x: tryconvert(x.split('.')[0], np.NaN, int), 'filename': lambda x: x} # conversion to appropriate datatypes
    apply_conversion = lambda x: {k: func(x[k]) for k, func in conversion_dict.items()}
    
    parsed_df = pd.DataFrame([apply_conversion(e) for e in parsed_dict])
    parsed_df = parsed_df.replace({'running': {'True': True, 'False': False}})
    
    return parsed_df