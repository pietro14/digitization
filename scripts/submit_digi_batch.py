#!/bin/env python
# USAGE: python3.8 scripts/submit_batch.py $PWD --outdir outdir --inputdir inputdir  --dry-run
# remove --dry-run to submit for real (otherwise only the scripts are created and commands are printed)

jobstring  = '''#!/bin/bash
#PBS -j oe  
source /nfs/cygno/users/$USER/CYGNO/setup_cygno_login.sh

cd CYGNOBASE
DIGISTRING
'''

import os, sys, re

if __name__ == "__main__":
    
    from optparse import OptionParser
    parser = OptionParser(usage='%prog workdir runs [options] ')
    parser.add_option(        '--dry-run',  dest='dryRun' ,  action='store_true', default=False, help='Do not run the job, only print the command');
    parser.add_option(        '--outdir',   dest='outdir',   type="string",       default=None,  help='outdirectory');
    parser.add_option(        '--inputdir', dest='inputdir', type="string",       default=None,  help='inputdirectory');
    parser.add_option(        '--config',   dest='config',     type="string",       default=None,  help='config file');
    parser.add_option(        '--tag',      dest='tag',      type="string",       default=None,  help='tag');
    #parser.add_option(        '--nthreads', dest='nthreads', type="string", default=24, help='number of threads / job');
    (options, args) = parser.parse_args()

    if len(args)<1:
        parser.print_help()
        exit(1)

    abswpath  = os.path.abspath(args[0]) 
    if not os.path.isdir(abswpath):
        raise RuntimeError('ERROR: {p} is not a valid directory. This is the base dir where the jobs run'.format(p=abswpath))

    if not options.outdir:
        raise RuntimeError('ERROR: give at least an output directory. This is where the jobs configs and logs are stored')
    else:
        absopath  = os.path.abspath(options.outdir)
        if not os.path.isdir(absopath):
            print ('making a directory ',absopath,' and running in it')
            os.system('mkdir -p {od}'.format(od=absopath))

    jobdir = absopath+'/jobs/'+options.tag+'/'
    if not os.path.isdir(jobdir):
        os.system('mkdir {od}'.format(od=jobdir))
    os.system('cp {cfg} {od}/config_{tag}.txt'.format(cfg=options.config,od=jobdir,tag=options.tag))
    logdir = absopath+'/logs/'+options.tag+'/'
    if not os.path.isdir(logdir):
        os.system('mkdir {od}'.format(od=logdir))
    outdir = absopath+'/out/'+options.tag+'/'
    if not os.path.isdir(logdir):
        os.system('mkdir {od}'.format(od=outdir))
    

    commands = []
    #for run in runs:
    
    job_file_name = jobdir+'/job_{tag}.sh'.format(tag=options.tag)
    log_file_name = logdir+'/job_{tag}.log'.format(tag=options.tag)
    tmp_file = open(job_file_name, 'w')

    tmp_filecont = jobstring
    cmd = 'python3.8 MC_data_new.py '+options.config+' -I '+options.inputdir+' -O '+outdir
    tmp_filecont = tmp_filecont.replace('DIGISTRING',cmd)
    tmp_filecont = tmp_filecont.replace('CYGNOBASE',abswpath+'/')
    tmp_file.write(tmp_filecont)
    tmp_file.close()
    
    RAM=4000
    #RAM=16000

    #sub_cmd = 'qsub -q cygno-custom -d {dpath} -l mem={ram}mb -o localhost:{logf} {jobf}'.format(dpath=abswpath,ram=RAM,logf=log_file_name,jobf=job_file_name)
    sub_cmd = 'qsub -q cygno -l mem={ram}mb -d {dpath} -o localhost:{logf} {jobf}'.format(dpath=abswpath,ram=RAM,logf=log_file_name,jobf=job_file_name)

    #sub_cmd = 'qsub -q cygno-custom -l nodes=1:disk5 -d {dpath} -o localhost:{logf} {jobf}'.format(dpath=abswpath,logf=log_file_name,jobf=job_file_name)
    commands.append(sub_cmd)

    if options.dryRun:
        for c in commands:
            print (c)
    else:
        for c in commands:
            os.system(c)

    print ("DONE")


        
