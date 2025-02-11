import os
for year in ['2016', '2017', '2018', 'Run2']:
    for mass in os.listdir(f'datacard/{year}'):
        prefix = f'datacard/{year}/{mass}'
        for file in set(os.listdir(prefix)):
            name = file.replace('.txt', '')
            #os.system(f"cd {prefix}; combineCards.py {name}={file} > {file}")
            if file.endswith('.txt') and '1' in file:
                os.system(f"cd {prefix}; combineCards.py {name}={file} {name.replace('1', '2')}={file.replace('1', '2')} > {file.replace('1', '')}")