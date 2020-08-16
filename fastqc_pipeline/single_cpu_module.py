import os, sys,  time, shlex, subprocess, module_helper, helper

def initialize_single_cpu(param):
    #split module list and add module load to the list
    helper.checkParameter(param,key='modules',dType=str)
    param['modules']=[('module load '+mod.strip()) for mod in param['modules'].split(",")]

        
def run_single_job(index, param, job_id, py_file, cores):
    #fetch modules that need to be loaded and add the python command for a single sample
    #initialize_single_cpu(param)
    #command_list = param['modules'][:]
    #command_list.append('python2.7 '+ param['scripts_dir']+ py_file +' -i ' + str(index) + ' -n $NSLOTS')
    command_list=('python2.7 '+ param['scripts_dir']+ py_file +' -i ' + str(index) + ' -n $NSLOTS')

    print(command_list)
    #call qsub script
    for cmd in command_list:
        args = shlex.split(cmd)
        output,error = subprocess.Popen(args,stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()

