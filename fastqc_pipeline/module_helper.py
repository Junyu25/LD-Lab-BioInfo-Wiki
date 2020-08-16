import os, sys

def checkParameter(param,key,dType,allowed=[],checkFile=False,optional=False):
    #generic function that checks if a parameter was in the parameter file, casts to the right data type and if the parameter is a file/ directory checks if it actually exists

    if optional and key not in param:
        param[key]=''
    else:
        #check if key is present
        if key not in param: print('Parameter '+key+' is missing in the parameter file.'); sys.exit(0)

        #cast to correct data type
        if dType == bool:
            param[key] = param[key] in ['True','TRUE','true','T','1']
        else:
            param[key]=dType(param[key])

        #check if the value for the key is allowed
        if len(allowed)>0 and param[key] not in allowed:
            print('Parameter '+key+' can only be one of the following: '+allowed); sys.exit(0)
    
        #if file or directory check if it exists
        if checkFile and not os.path.exists(param[key]):
            print('The file in '+key+' = ',param[key],' does not exist'); sys.exit(0)
    
    
def getPercentage(a,b,n):
    percent=[0.0]*n    
    for idx in range(n):
        percent[idx]=round(float(a[idx])/float(b[idx])*100,2)  
    return [str(round(pc,1))+'%' for pc in percent]
            
def initialize_module():
    #Import modules
    import sys, getopt, json

    working_dir='./'
    #Check arguments
    if len(sys.argv)<5:
        print('ERROR: Specify the index of the file the parameter should be run on.')
        sys.exit(0)
    optlist,cmdlist = getopt.getopt(sys.argv[1:],'i:n:d:')
    for opt in optlist:
        if opt[0]=='-i': file_index=opt[1]
        if opt[0]=='-n': num_processors=opt[1]
        if opt[0]=='-d': working_dir=opt[1]
    print('###########################')
    print(sys.argv)
    print(working_dir)

    #Read and initialize parameters
    with open(working_dir+'results/parameters.json') as f:
        param = json.load(f)
    param['file_index']=int(file_index)
    param['num_processors']=num_processors
    
    #use the input files that were specified in the pipeline call
    param['working_file']=param[param['input_files']][param['file_index']]
    if param['paired'] and param['input_files']+'2' in param:
       param['working_file2']=param[param['input_files']+'2'][param['file_index']]

    #name of the log file
    log_file = param['working_dir']+'results/log/'+param['stub'][param['file_index']]+'.log'

    #check if it is not a clean install and if module already completed
    param['resume_module'] = False
    if not param['clean_run']:
        param['file_handle'] = open(log_file)
        #check if the module already finished
        lines_end = [line for line in param['file_handle'].readlines() if 'ENDING %s |' %(param['current_flag']) in line.rstrip()]
        if len(lines_end)>0:  
            param['file_handle'].close()
            #open file for writing 
            param['file_handle'] = open(log_file,'a')
            param['file_handle'].write(param['current_flag']+' module already run on this file .. SKIPPING\n')  
            param['file_handle'].close()
            sys.exit(0)
        #check if the module was started, but not finished the enable resuming
        lines_start = [line for line in param['file_handle'].readlines() if 'STARTING %s |' %(param['current_flag']) in line.rstrip()]
        if len(lines_start)>0:
            param['resume_module'] = True
            
    #start process log
    param['file_handle'] = open(log_file,'a')
    param['file_handle'].write('STARTING '+param['current_flag']+'\n')
    return param
     
def outputPhenotype(param,pheno_file):
    #write sample info only for the files that successfully completed
    header=param['raw_file_header']
    index=0
    out = open(pheno_file, "w")
    f = open (param['raw_filenames'] , 'r')
    for line in f:
        #take header into account
        if header:
            out.write(line)
            header=False
        else:
            #write only samples that successfully finished
            if param['run_log'][-1][index]:
                out.write(line)
            index+=1    
    out.close()
    f.close()

def output_sample_info(param):
    #create a file with the phenotype data of the samples that actually made it through HTSeq
    param['pheno_file']=param['working_dir']+'deliverables/sample_info.txt'
    outputPhenotype(param,param['pheno_file'])


def rotate_word(word,deg=270):
    return '<div style="float: center;position: relative;-moz-transform: rotate(270deg);  /* FF3.5+ */-o-transform: rotate(270deg);  /* Opera 10.5 */ -webkit-transform: rotate('+str(deg)+'deg);  /* Saf3.1+, Chrome */ filter:  progid:DXImageTransform.Microsoft.BasicImage(rotation=3);  /* IE6,IE7 */ -ms-filter: progid:DXImageTransform.Microsoft.BasicImage(rotation=3); /* IE8 */">' +word+'</div>'         

def writeHTMLtable(param,table,out,fCol_width=200,cell_width=50,initial_breaks=8,deg=315):
    out.write('<table id="one-column-emphasis" class="fixed"><col width="'+str(fCol_width)+'px"/>\n')
    out.write(''.join(['<br>']*initial_breaks)+'\n')
    out.write(''.join(['<col width="'+str(cell_width)+'px"/>\n']*len(table[0])))
    #write header
    param['report'].write('<thead><tr><th></th>'+''.join(['<th>'+rotate_word(stub.replace('-','_'),deg)+'</th>\n' for stub in table[0]])+'</tr></thead>\n')
           
    #write the pass and fail for each module
    for idx in range(1,len(table)):
        out.write('<tr>')      
        if len(table[idx])==1:
            out.write('<th colspan="'+str(1+len(table[0]))+'">'+table[idx][0]+'</th>\n')
        else:
            for i in range(len(table[idx])):
                out.write('<td>'+str(table[idx][i])+'</td>\n')         
        param['report'].write('</tr>\n')  
       
    #close table
    out.write('</table><br>')

def wrapup_module(param, new_working_file=[]):    
    #end process log
    param['file_handle'].write('ENDING '+param['current_flag']+' | ')
    #if there was an actual output file specified
    if len(new_working_file)>0:
        param['file_handle'].write(';'.join([w for w in new_working_file]))
    param['file_handle'].write('\n\n')
    param['file_handle'].close()