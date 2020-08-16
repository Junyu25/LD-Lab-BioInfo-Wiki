import os, sys,  time, shlex, subprocess, module_helper, helper, random

def initialize_qsub(param):
    #split module list and add module load to the list
    helper.checkParameter(param,key='modules_python2.6',dType=str)
    param['modules_python2.6']=[('module load '+mod.strip()) for mod in param['modules_python2.6'].split(",")]
    helper.checkParameter(param,key='modules_python2.7',dType=str)
    param['modules_python2.7']=[('module load '+mod.strip()) for mod in param['modules_python2.7'].split(",")]
    helper.checkParameter(param,key='qsub_email',dType=str)
    helper.checkParameter(param,key='qsub_send_email',dType=bool)
    helper.checkParameter(param,key='qsub_memory',dType=str)
    helper.checkParameter(param,key='qsub_suffix',dType=str)
    helper.checkParameter(param,key='qsub_PROJECT',dType=str)
    helper.checkParameter(param,key='qsub_MACHINE',dType=str)
    helper.checkParameter(param,key='qsub_RUNTIME_LIMIT',dType=str)
    helper.checkParameter(param,key='qsub_wait_time',dType=int)
    helper.checkParameter(param,key='qsub_num_processors',dType=str)


def wait_for_qsub(param, job_id):
    user,error = subprocess.Popen(args = shlex.split('whoami'),stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()
    user = user.rstrip()
    
    #wait until there are no jobs with the current job name in the queue anymore
    qjobs = param['num_samples']
    helper.writeLog('Waiting for single '+param['current_flag']+' jobs to finish.... \n',param)
    while qjobs > 0 :
        time.sleep(param['qsub_wait_time'])
        # check the number of jobs
        cmd = 'qstat -u $USER'
        args = shlex.split(cmd)
        output,error = subprocess.Popen(args,stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()

        #add grabbing the job IDs
        if len(output) == 0:
            qjobs = 0
        else:    
            for current_line in output.split('\n'):
                if ' '+user+' ' in current_line and ' '+job_id+' ' in current_line:
                    qjobs+=1  
        
def submit_jobs(index, param, py_file, job_id, cores, environment):
    #fetch modules that need to be loaded and add the python command for a single sample
    command_list = param[environment][:]
    #if we load the python2.7 module we want to use python2.7 rather than python
    if environment == 'modules_python2.7':
        cmd='python2.7 '
    else:      
        #if we need to run python2.6 just unload python2.7, otherwise it results in errors
        command_list.append('module unload python2.7/Python-2.7.3_gnu446')
        cmd='python '
    
    cmd=cmd+param['scripts_dir']+ py_file +' -i ' + str(index) + ' -n $NSLOTS' + ' -d ' +  param['working_dir']
    command_list.append(cmd)


    #make a directory for all the qsub commands
    param['qsub_dir']=param['working_dir']+'results/qsub/'
    if not os.path.exists(param['qsub_dir']):
        os.makedirs(param['qsub_dir'])
    
    qsub_dir=param['working_dir']+'results/qsub/'
    qsub_filename=qsub_dir+(param['stub'])[index]+'.qsub'
    
    outhandle=PBS_File_class(qsub_filename)
    outhandle.qsub_dir=qsub_dir
    outhandle.job_id = job_id 
    outhandle.email = param['qsub_email']
    outhandle.send_email = param['qsub_send_email']
    outhandle.memory = param['qsub_memory']
    outhandle.suffix = param['qsub_suffix']
    outhandle.PROJECT = param['qsub_PROJECT']
    outhandle.MACHINE = param['qsub_MACHINE']
    outhandle.RUNTIME_LIMIT = param['qsub_RUNTIME_LIMIT']
    outhandle.output_pbs(command_list)

    #call qsub script
    cmd = 'qsub -pe single_node '+cores+' '+qsub_filename
    args = shlex.split(cmd)
    output,error = subprocess.Popen(args,stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()



class PBS_File_class():
    def __init__(self,name,path=os.getcwd()):
        self.filename=name

    def output_handle_gen(self,FILEPATH=os.getcwd()):
        if '/' in self.filename:
            complete_path=self.filename
        else:
            complete_path=FILEPATH + '/' + self.filename
        self.handle=open(complete_path,'w')

    def output_pbs(self,command_line_list):
        self.output_handle_gen()

        if self.MACHINE=="scc":
            self.handle.write("source ~/.bashrc\n")
            self.handle.write("#!bin/bash\n")
            self.handle.write("#$ -l h_rt="+self.RUNTIME_LIMIT+'\n')
            self.handle.write("\n")
        else:
            self.handle.write("#!bin/bash\n")
            self.handle.write("#\n")
            self.handle.write("\n")

        self.handle.write("#Specify which shell to use\n")
        self.handle.write("#$ -S /bin/bash\n")
        self.handle.write("\n")

        self.handle.write("#Run on the current working folder\n")
        self.handle.write("#$ -cwd\n")
        self.handle.write("\n")

        self.handle.write("#Give this job a name\n")
        self.handle.write("#$ -N "+self.job_id+'\n')
        self.handle.write("\n")

        self.handle.write("#Join standard output and error to a single file\n")
        self.handle.write("#$ -j y\n")
        self.handle.write("\n")

        self.handle.write("# Name the file where to redict standard output and error\n")
        if self.filename.count("/")>=1:
            filename_info_list=self.filename.split("/")
            filename_info=filename_info_list[-1]
        else:
            filename_info=self.filename
        self.handle.write("#$ -o "+ self.qsub_dir + filename_info +".qlog\n")
        self.handle.write("\n")

        self.handle.write("# Project this job belongs to \n")
        self.handle.write("#$ -P " + self.PROJECT+ " \n")
        self.handle.write("\n")

        if  self.email!="" and self.send_email:
            self.handle.write("# Send an email when the job begins and when it ends running\n")
            self.handle.write("#$ -m be\n")
            self.handle.write("\n")

            self.handle.write("# Whom to send the email to\n")
            self.handle.write("#$ -M "+self.email+ "\n")
            self.handle.write("\n")

        self.handle.write("# memory usage\n")
        self.handle.write("#$ -l mem_free="+self.memory+ "\n")
        self.handle.write("\n")

        self.handle.write("# Now let's Keep track of some information just in case anything go wrong\n")
        self.handle.write("echo "+'"'+"========================================" + '"'+'\n')
        self.handle.write("echo "+'"'+"Starting on : $(date)"+'"'+ "\n")
        self.handle.write("echo "+'"'+"Running on node : $(hostname)"+'"'+"\n")
        self.handle.write("echo "+'"'+"Current directory : $(pwd)"+'"'+"\n")
        self.handle.write("echo "+'"'+"Current job ID : $JOB_ID"+'"'+"\n")
        self.handle.write("echo "+'"'+"Current job name : $JOB_NAME"+'"'+"\n")
        self.handle.write("echo "+'"'+"Task index number : $TASK_ID"+'"'+"\n")
        self.handle.write("echo "+'"'+"========================================" + '"'+'\n')
        self.handle.write("\n")

        for command_line in command_line_list:
            self.handle.write(command_line)
            self.handle.write('\n')

        self.handle.write("\n")
        self.handle.write("echo "+'"'+"========================================" + '"'+'\n')
        self.handle.write("echo "+'"'+"Finished on : $(date)"+'"'+ "\n")
        self.handle.write("echo "+'"'+"========================================" + '"'+'\n')
        self.handle.close()