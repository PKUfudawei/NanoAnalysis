import os
for i in os.listdir('.'):
    if '123456' in i:
        os.system(f"mv {i} {i.replace('.123456', '')}")
    if '42' in i:
        os.system(f"mv {i} {i.replace('.42', '')}")
