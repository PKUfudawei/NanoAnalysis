import json
import os
import os.path
import re
import subprocess
import sys
import random
import yaml


def get_file_info(ds):
    """
    The files might be from
    RUCIO: rucio list-file-replicas --pfn cms:/DATASET
    EOS: gsiftp://eosuserftp.cern.ch/eos/user/s/somename
    Local: direct available
    DAS: use dasgoclient
    """
    # print("===> get_file_info", "sample:", sample)
    if ds.startswith("gsiftp://"):  # on EOS
        outinfo = subprocess.run(args=f"gfal-ls {ds}", shell=True, stdout=subprocess.PIPE, encoding='utf-8')
        file_list = outinfo.stdout.split("\n")[:-1]
        # file_list = sorted(file_list)  # sort the files
        file_dict = {}
        for ifname in file_list:
            file_path = f"{ds}/{ifname}"
            try:
                proc = subprocess.Popen([f"gfal-sum {file_path} ADLER32"], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE, encoding='utf-8')
                out, err = proc.communicate(input=None)
                local_sum = out.strip("\n").split(" ")[-1]
                file_dict[ifname] = {
                    'checksum': local_sum,
                    'filelist': [file_path],
                }
            except:
                print(f"xxx [Error]: could not run gfal-sum, please check, e.g., unset cmsenv")
                exit(1)
        use_gfal = True
    elif ds.startswith("local:"):  # on local
        outinfo = subprocess.run(args=f"ls {ds.lstrip('local:')}", shell=True,   stdout=subprocess.PIPE,encoding='utf-8')
        file_list = outinfo.stdout.split("\n")[:-1]
        # file_list = sorted(file_list)  # sort the files
        file_dict = {}
        for ifname in file_list:
            file_path = f"{ds.lstrip('local:')}/{ifname}"
            try:
                proc = subprocess.Popen([f"gfal-sum {file_path} ADLER32"], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE, encoding='utf-8')
                out, err = proc.communicate(input=None)
                local_sum = out.strip("\n").split(" ")[-1]
                file_dict[ifname] = {
                    'checksum': local_sum,
                    'filelist': [file_path],
                }
            except:
                print(f"xxx [Error]: could not run gfal-sum, please check, e.g., unset cmsenv")
                exit(1)
        use_gfal = False
    else:  # from DAS
        #########################################################################
        ############################ Use dasgoclient ############################
        #########################################################################
        query_str = f"\"file dataset={ds} system=rucio\""
        print(f"/cvmfs/cms.cern.ch/common/dasgoclient -query={query_str} -json")
        proc = subprocess.Popen([f"/cvmfs/cms.cern.ch/common/dasgoclient -query={query_str} -json"], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE, encoding='utf-8')
        out, err = proc.communicate(input=None)
        out_str = json.loads(out)

        file_dict = {}
        for ifile_block in out_str:
            fname_tmp = ifile_block['file'][0]['name'].split('/')[-1]
            file_dict[fname_tmp] = {
                'checksum': ifile_block['file'][0]['adler32'],
                'filelist': list(ifile_block['file'][0]['pfns'].keys()),
            }
            # remove files from bad site: cms.sscc.uos.ac.kr
            file_dict[fname_tmp]['filelist'] = [
                ifile_path for ifile_path in file_dict[fname_tmp]['filelist'] if not "cms.sscc.uos.ac.kr" in ifile_path
            ]
        use_gfal = True

    print("<=== get_file_info END :)")
    return use_gfal, file_dict



def get_sub(job_path, sample, job_flavour="testmatch",idx=0,long_queue=None):
    if long_queue is None:
        long_queue = []
    # print("===> get_sub", "job_path", job_path, "sample:", sample, "idx:", idx, "long_queue:", len(long_queue) != 0)
    if len(long_queue) == 0:
        job_name = f"{sample}_{idx}"
        file_content = f"executable = {job_path}/{job_name}.sh\n"
        file_content += "universe = vanilla\n"
        file_content += f"output = {job_path}/{job_name}.out\n"
        file_content += f"error = {job_path}/{job_name}.err\n"
        file_content += f"log = {job_path}/{job_name}.log\n"
        # file_content += "request_cpus = 1\n"
        # file_content += "request_memory = 1024\n"
        # file_content += "request_disk = 1024\n"
        # file_content += "requirements = (machine == \"atlas.phy.pku.edu.cn\") || (machine == \"farm.phy.pku.edu.cn\") || (machine == \"node01.phy.pku.edu.cn\") || (machine == \"node02.phy.pku.edu.cn\") || (machine == \"node03.phy.pku.edu.cn\") || (machine == \"node04.phy.pku.edu.cn\") || (machine == \"node05.phy.pku.edu.cn\") || (machine == \"node06.phy.pku.edu.cn\")\n"
        file_content += f"+JobFlavour = \"{job_flavour}\"\n"
        file_content += "queue\n"
        tmp = open(f"{job_path}/{job_name}.sub", "w")
        tmp.write(file_content)
    else:
        file_content = f"executable = {job_path}/$(JNAME).sh\n"
        file_content += "arguments = $(JNAME) $(INFILE)\n"
        file_content += "universe = vanilla\n"
        file_content += f"output = {job_path}/$(JNAME).out\n"
        file_content += f"error = {job_path}/$(JNAME).err\n"
        file_content += f"log = {job_path}/$(JNAME).log\n"
        # file_content += "request_cpus = 1\n"
        # file_content += "request_memory = 1024\n"
        # file_content += "request_disk = 1024\n"
        # file_content += "requirements = (machine == \"atlas.phy.pku.edu.cn\") || (machine == \"farm.phy.pku.edu.cn\") || (machine == \"node01.phy.pku.edu.cn\") || (machine == \"node03.phy.pku.edu.cn\")\n"
        file_content += f"+JobFlavour = \"{job_flavour}\"\n"
        file_content += "queue JNAME in (\n"
        for i in long_queue:
            file_content += f"{i}\n"
        file_content += ")\n"

        # print(file_content)
        proc = subprocess.Popen(["condor_submit"], shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
        out, err = proc.communicate(input=file_content)
        if proc.returncode != 0:
            sys.stderr.write(err)
            raise RuntimeError('Job submission failed.')
        print(out.strip())

        matches = re.match('.*submitted to cluster ([0-9]*).', out.split('\n')[-2])
        if not matches:
            sys.stderr.write('Failed to retrieve the job id. Job submission may have failed.\n')
            for i in long_queue:
                jidFile = f"{job_path}/{i}.jid"
                open(jidFile, 'w').close()
        else:
            clusterId = matches.group(1)
            # now write the jid files
            proc = subprocess.Popen(['condor_q', clusterId, '-l', '-attr', 'ProcId,Cmd', '-json'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
            out, err = proc.communicate(input=None)
            try:
                qlist = json.loads(out.strip())
            except:
                sys.stderr.write('Failed to retrieve job info. Job submission may have failed.\n')
                for jName in long_queue:
                    jidFile = f"{job_path}/{jName}.jid"
                    open(jidFile, 'w').close()
            else:
                for qdict in qlist:
                    with open(qdict['Cmd'].replace('.sh', '.jid'), 'w') as out:
                        out.write('%s.%d\n' % (clusterId, qdict['ProcId']))
    # print("<=== get_sub END :)")
    return


def get_sh(idx, use_gfal, fname, fsum, files, sample, job_path, out_path, sleep):
    HOME_PATH = os.environ['HOME']
    # print("===> get_sh")
    file_content = "#!/bin/bash\n"
    file_content += f"export X509_USER_PROXY={HOME_PATH}/.proxy\n"
    file_content += "voms-proxy-info\n"
    file_content += f"TO_SLEEP={round(random.uniform(0, float(sleep)), 2)}m\n"
    file_content += "echo \"===> Sleep How Long: ${TO_SLEEP}in\"\n"
    file_content += "sleep ${TO_SLEEP}\n"
    # file_content += "source /cvmfs/cms.cern.ch/cmsset_default.sh\n"
    file_content += "# Load data\n"
    # file_content += "echo $TMPDIR\n"
    # file_content += f"cd {work_dir}\n"
    file_content += "echo '===> PWD: '$PWD\n"
    file_content += "export WORK_PATH=$PWD\n"
    # checksum
    file_content += f"export FILECHECKSUM={fsum}\n"
    if use_gfal:
        if len(files) > 1:  # Should from DAS, several replicas are available, try each until successfully download one
            for_str = "for ifile in"
            for i, ifile in enumerate(files):
                file_content += f"FILE{i}=\"{ifile}\"\n"
                for_str += f" $FILE{i}"
            file_content += f"{for_str}; do\n"
            file_content += "    sleep 2s\n"
            file_content += f"    export X509_USER_PROXY={HOME_PATH}/.proxy\n"
            file_content += f"    gfal-copy -t 86400 -T 86400 -f $ifile {out_path}\n"
            file_content += "    if [ $? -eq 0 ]; then break; fi\n"
            file_content += "done\n"
            file_content += "\n"
        else:  # Maybe from eos, try 5 times
            file_content += "for itime in 1 2 3 4 5; do\n"
            file_content += f"    export X509_USER_PROXY={HOME_PATH}/.proxy\n"
            file_content += f"    gfal-copy -t 86400 -T 86400 -f {files[0]} {out_path}\n"
            file_content += "    if [ $? -eq 0 ]; then break; fi\n"
            file_content += "done\n"
            file_content += "\n"
    else:
        file_content += f"cp {files[0]} {out_path}\n"
        file_content += "\n"
    file_content += "# Checksum\n"
    file_content += f"check_file=`gfal-sum {out_path}/{fname} ADLER32`\n"
    file_content += "[ ${check_file#* } = ${FILECHECKSUM} ]\n"
    file_content += f"[ $? -eq 0 ] && mv {job_path}/{sample}_{idx}.jid {job_path}/{sample}_{idx}.done\n"
    tmp = open(f"{job_path}/{sample}_{idx}.sh", "w")
    # print(f"{job_path}/{sample}_{idx}.sh")
    tmp.write(file_content)
    os.system(f"chmod +x {job_path}/{sample}_{idx}.sh")
    # print("<=== get_sh END :)")
    return


def submit_command(cfg, sample=None):
    print("===> submit_command", "year:", cfg['year'], "sample:", sample)
    use_gfal, file_dict = get_file_info(cfg['ds_yml'][sample]['dataset'])
    print("---", "sample:", sample, ", use_fal:", use_gfal, ", #file:", len(file_dict))

    job_path = f"{cfg['job_dir']}/{cfg['sample_type']}/{cfg['year']}/{cfg['channel']}/{cfg['nano_ver']}/{sample}/"
    # sample dataset name
    sname = cfg['ds_yml'][sample]['dataset']
    if sname.startswith("gsiftp://") or sname.startswith("local:"):
        if sname.endswith("/"):
            sname = sname.split("/")[-2]
        else:
            sname = sname.split("/")[-1]
        out_path = f"{cfg['out_dir']}/{cfg['sample_type']}/{cfg['year']}/{cfg['channel']}/private/{cfg['nano_ver']}/{sname}/"
    else:
        out_path = f"{cfg['out_dir']}/{cfg['sample_type']}/{cfg['year']}/{cfg['channel']}/{cfg['nano_ver']}/{sname}/"
    
    base_out_path_list = cfg['out_dir'].split("/")
    base_out_path_list = [[i for i in base_out_path_list if not i==""]]
    out_path_list = out_path.split("/")
    out_path_list = [[i for i in out_path_list if not i==""]]

    if not os.path.exists(job_path):
        os.makedirs(job_path)
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    # make the files can be changed by others
    
    for i in range(len(base_out_path_list),len(out_path_list)):
        os.system(f"chmod g+w {'/'.join(out_path_list[:i])}")
    
    _long_queue = []
    for idx,ifname in enumerate(file_dict):
        out_file_name = f"{out_path}/{ifname}"
        to_do_submit = True
        if os.path.exists(out_file_name):
            if cfg['checkfile']:
                try:
                    proc = subprocess.Popen([f"gfal-sum {out_file_name} ADLER32"], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE, encoding='utf-8')
                    out, err = proc.communicate(input=None)
                    local_sum = out.strip("\n").split(" ")[-1]
                    # print("--- [Test]: checksum:",file_dict[ifname]['checksum'],"local_sum",local_sum)
                    if file_dict[ifname]['checksum'] == local_sum:
                        print("--- The job done already, check:", out_file_name)
                        to_do_submit = False
                    else:
                        print(f"xxx Remove bad file: {out_file_name}")
                        os.remove(out_file_name)
                        to_do_submit = True
                except:
                    print(f"xxx [Error]: could not run gfal-sum, please check, e.g., unset cmsenv")
                    exit(1)
            else:
                print("--- The job done already, check:", out_file_name)
                to_do_submit = False
                pass
        else:
            pass
        if to_do_submit:
            # We write the SUB file for documentation / resubmission, but initial submission will be done in one go below
            get_sub(job_path, sample, cfg['jflavour'],idx)
            get_sh(idx, use_gfal, ifname, file_dict[ifname]['checksum'], file_dict[ifname]['filelist'], sample, job_path, out_path, cfg['sleep'])
            _long_queue.append(f"{sample}_{idx}")
    if cfg['dryrun']:
        pass
    else:
        get_sub(job_path, sample, cfg['jflavour'], long_queue=_long_queue)  # dryrun will just generate submit files, but not run
    print("<=== submit_command END :)")
    return True


def status_command(cfg, sample):
    print("===> status_command", "year:", cfg['year'], "sample:", sample)
    use_gfal, file_dict = get_file_info(cfg['ds_yml'][sample]['dataset'])
    print("---", "sample:", sample, ", use_fal:", use_gfal, ", #file:", len(file_dict))

    # print(file_dict)
    sname = cfg['ds_yml'][sample]['dataset']
    if sname.startswith("gsiftp://") or sname.startswith("local:"):
        if sname.endswith("/"):
            sname = sname.split("/")[-2]
        else:
            sname = sname.split("/")[-1]
        out_path = f"{cfg['out_dir']}/{cfg['nano_ver']}/private/{sname}/"
    else:
        out_path = f"{cfg['out_dir']}/{cfg['nano_ver']}/{sname}/"

    done_file = []
    not_done_file = []
    for idx,ifname in enumerate(file_dict):
        out_file_name = f"{out_path}/{ifname}"
        if os.path.exists(out_file_name):
            if cfg['checkfile']:
                try:
                    proc = subprocess.Popen([f"gfal-sum {out_file_name} ADLER32"], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE, encoding='utf-8')
                    out, err = proc.communicate(input=None)
                    local_sum = out.strip("\n").split(" ")[-1]
                    if file_dict[ifname]['checksum'] == local_sum:
                        print("--- The job done already, check:", out_file_name)
                        done_file.append(out_file_name)
                    else:
                        print(f"xxx Remove bad file: {out_file_name}")
                        os.remove(out_file_name)
                except:
                    print(f"xxx [Error]: could not run gfal-sum, please check, e.g., unset cmsenv")
                    exit(1)
            else:
                print("--- The job done already, check:", out_file_name)
                done_file.append(out_file_name)
                pass
        else:
            not_done_file.append(out_file_name)

    print("--- sample:", sample, ", total files:", len(file_dict), ", done:", len(done_file), "not done:", len(not_done_file), ", they are following jobs:")
    for i in not_done_file:
        print("\t-", i)
    print("<=== submit_command END :)")
    return True

