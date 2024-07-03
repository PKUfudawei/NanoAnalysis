# Instruction about how to run `Analyze.ipynb` on lxplus

To get various plots at the post-process step, we just need to run `Analyze.ipynb` cell by cell. The detailed comments and instructions are also recorded in that file. 

Now let's learn how to prepare the jupyter environment.

## Prerequisite 

1. It is recommended to use [vscode](https://code.visualstudio.com/), please download it.

2. Install the necessary extensions/plugins like *ssh*, *python*, and *jupyter* on vscode. Please follow the tutorial "[Install an extension](https://code.visualstudio.com/docs/editor/extension-marketplace#_install-an-extension)".

3. At your home diretory `~` on lxplus: run `cmsrel CMSSW_14_0_0` to create `CMSSW_14_0_0/` directory.

4. Add the line `cd ~/CMSSW_14_0_0/src; cmsenv; cd ~` to the file `~/.bashrc`.

Now we have prepared the python3 and jupyter environment. Every time we connect to lxplus via vscode, we automatically enter CMSSW_14_0_0 environment with the needed python3 and jupyter.

## How to post-process and get plots

1. Connect to lxplus via ssh on vscode, please follow the tutorial "[Remote Development using SSH](https://code.visualstudio.com/docs/remote/ssh)".

2. Open `NanoAnalysis/postprocess/Analyze.ipynb` and run it cell by cell (keyboard shortcut: Shift+Enter).

