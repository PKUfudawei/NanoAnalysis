import os
import json


def genWeight_reweighting(genWegiht, JSON_dir):  # .json files stored statistics like n_events
    n_raw_events = 0
    for (current_path, dirs, files) in os.walk(JSON_dir):
        for f in files:
            if f.endswith('.json'):
                with open(f, 'r', encoding ='utf-8') as f:
                    stats = json.load(f)
                n_raw_events += list(stats.values())[0]['n_events']
                
    return genWegiht / n_raw_events
