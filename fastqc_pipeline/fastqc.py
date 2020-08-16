import os, subprocess, csv
import module_helper



def extractTables(param):
    #get the rownames
    csv_file = open(param['fastqc_dir']+param['fastqc_stub'][0]+'_fastqc/summary.txt')
    csv_reader = csv.reader(csv_file, delimiter='\t')
    #get the rownames    
    checks = [[row[1]] for row in csv_reader]
    csv_file.close()
  
    #links to the icons
    PASS='<img src="Icons/tick.png">'
    FAIL='<img src="Icons/error.png">'
    WARN='<img src="Icons/warning.png">'

    #get the values for each sample (icons for pass, faile or warning)
    for idx in range(len(param['fastqc_stub'])):
        csv_file = open(param['fastqc_dir']+param['fastqc_stub'][idx]+'_fastqc/summary.txt')
        overview_file='fastqc/'+param['fastqc_stub'][idx]+'_fastqc/fastqc_report.html#M'
        csv_reader = csv.reader(csv_file, delimiter='\t')
        i=0
        for row in csv_reader:
            cell='<a href="'+overview_file+str(i)+'">'
            if row[0]=='PASS': cell=cell+PASS
            if row[0]=='FAIL': cell=cell+FAIL
            if row[0]=='WARN': cell=cell+WARN
            cell=cell+'</a>'
            checks[i].append(cell)
            i+=1
        csv_file.close()
    return checks


def createOverviewTable(param):
        
    #create a table
    table=[]
    #put in headers
    table.append([stub for stub in param['fastqc_stub']])
    #link to summary files
    table.append(['Summary files']+['<a href="fastqc/'+stub+'_fastqc/fastqc_data.txt">raw</a>' for stub in param['fastqc_stub']])
    #link to overview files
    table.append(['Full report']+['<a href="fastqc/'+stub+'_fastqc/fastqc_report.html"><img src="Icons/fastqc_icon.png"></a>' for stub in param['fastqc_stub']])
    #extract check marks
    table=table+extractTables(param)    
    #write the table as html     
    module_helper.writeHTMLtable(param,table,out=param['report'])
     
def readRawfastqc(param):
    summary_files = [param['fastqc_dir']+param['fastqc_stub'][idx]+'_fastqc'+'/summary.txt' for idx in range(len(param['fastqc_stub']))]
    
    fastqc=dict()
    #add entries into fastqc dictionary
    fh=open(summary_files[0])
    fastqc = dict([(name.split('\t')[1].strip(), []) for name in fh.readlines()])
    fh.close()

    #fill fastqc dictionary with information from the summary files
    for sum_file in summary_files:
        fh=open(sum_file)
        for name in (fh.readlines()):
            fastqc[name.split('\t')[1].rstrip()].append(name.split('\t')[0].strip())
        fh.close()
    
    key_list = fastqc.keys()
    
    #fill fastqc dictionary with information from the data file
    data_files = [param['fastqc_dir']+param['fastqc_stub'][idx]+'_fastqc'+'/fastqc_data.txt' for idx in range(len(param['fastqc_stub']))]
    labels = ['Encoding', 'Total Sequences', 'Filtered Sequences', 'Sequence length',	'%GC']
    
    fastqc.update(dict([(l, []) for l in labels]))
    key_list.extend(labels)
    
    for d_file in data_files:
        fh=open(d_file)
        for name in (fh.readlines()):
            if name.split('\t')[0].strip() in fastqc.keys():
                fastqc[name.split('\t')[0].strip()].append(name.split('\t')[1].rstrip())
        fh.close()
    
    param['fast_qc_summary']=fastqc

    # write overview file 
    fh=open(param['fastqc_dir']+'overview.txt','w')
    fh.write(' \t'+'\t'.join(param['fastqc_stub'])+'\n')
    for nam in key_list: fh.write(nam+'\t'+'\t'.join([str(vv) for vv in param['fast_qc_summary'][nam]])+'\n') 
    fh.close()


def copyFiles(param):
    #if there is no fastqc directory in the report make one
    param['fastqc_dir']=param['working_dir']+'report/fastqc/'
    if not os.path.exists(param['fastqc_dir']):
        os.makedirs(param['fastqc_dir'])
    
    call='ls ' + param['working_dir']+'results/fastqc/'
    output,error = subprocess.Popen(call.split(),stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()
    param['fastqc_stub']=[line.replace('_fastqc','') for line in output.split('\n') if '.zip' not in line and line != '' and '.html' not in line] 

    #copy the unpacked directories   
    for fastqc_file in param['fastqc_stub']:
        call='cp -R '+ param['working_dir']+'results/fastqc/'+fastqc_file+'_fastqc/ '+param['fastqc_dir']
        output,error = subprocess.Popen(call.split(),stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()

def plotNumberOfReads(param):
    import numpy as np
    import matplotlib.pyplot as plt
    import math, pylab
    #extract number of reads, these are needed not only here but also for the bamqc
    summary_files = [param['fastqc_dir']+param['fastqc_stub'][idx]+'_fastqc/fastqc_data.txt' for idx in range(len(param['fastqc_stub']))]
    num_total_reads=[open(sum_file).readlines()[6].split('\t')[1].strip().rstrip() for sum_file in summary_files]
    num_total_reads=[int(num) for num in num_total_reads]
    
    #if we deal with unpaired data the total number of reads used downstream is equivalent to the fastqc reads, but for paired it is the sum of the two paried files
    if not param['paired']:
        param['num_total_reads']=num_total_reads
    else:
        param['num_total_reads']=[0]*param['num_samples']
        samples_1=[param['raw_files'][idx].split('/')[-1].replace('.fastq','').replace('.fq','').replace('.gz','') for idx in range(param['num_samples'])]
        samples_2=[param['raw_files2'][idx].split('/')[-1].replace('.fastq','').replace('.fq','').replace('.gz','') for idx in range(param['num_samples'])]
        
        for idx in range(param['num_samples']):
            index1=[i for i in range(len(param['fastqc_stub'])) if param['fastqc_stub'][i]==samples_1[idx]][0]
            index2=[i for i in range(len(param['fastqc_stub'])) if param['fastqc_stub'][i]==samples_2[idx]][0]           
            param['num_total_reads'][idx]=num_total_reads[index1]+num_total_reads[index2]
        
    #create plot 
    fig, ax = plt.subplots()
    fig.set_size_inches(3+len(param['fastqc_stub'])*0.4,8)
    index = np.arange(len(num_total_reads))
    bar_width = 0.8
    opacity = 0.4
    rects1 = plt.bar(index, num_total_reads, bar_width,
                     alpha=opacity,
                     color='b')
    plt.xlabel('Samples')
    plt.ylabel('Total number of reads')
    plt.title('Total number of reads across samples')
    ticks=param['fastqc_stub']
    plt.xticks(index + bar_width, ticks,rotation='vertical')
    plt.tight_layout()
    
    #put it into the report
    filename='report/fastqc/total_reads.png'
    pylab.savefig(param['working_dir']+filename)
    param['report'].write('<img src="fastqc/total_reads.png" alt="total number of reads"><br><br>\n')

def plotGCContent(param):
    import numpy as np
    import matplotlib.pyplot as plt
    import math, pylab

    #extract number of reads, these are needed not only here but also for the bamqc
    summary_files = [param['fastqc_dir']+param['fastqc_stub'][idx]+'_fastqc/fastqc_data.txt' for idx in range(len(param['fastqc_stub']))]
    gc_content=[open(sum_file).readlines()[9].split('\t')[1].strip().rstrip() for sum_file in summary_files]
    gc_content=[int(num) for num in gc_content]
    
    #create plot 
    fig, ax = plt.subplots()
    fig.set_size_inches(3+len(param['fastqc_stub'])*0.4,8)
    index = np.arange(len(gc_content))
    bar_width = 0.8
    opacity = 0.4
    rects1 = plt.bar(index, gc_content, bar_width,
                     alpha=opacity,
                     color='b')
    plt.xlabel('Samples')
    plt.ylabel('%GC content')
    plt.title('GC content across samples')
    ticks=param['fastqc_stub']
    plt.xticks(index + bar_width, ticks,rotation='vertical')
    plt.tight_layout()
    ax.set_ylim(0, 100)
    
    #put it into the report
    filename='report/fastqc/gc_content.png'
    pylab.savefig(param['working_dir']+filename)
    param['report'].write('<img src="fastqc/gc_content.png" alt="GC content"><br><br>\n')


def plotAgvLengthOfReads(param):
    import numpy as np
    import matplotlib.pyplot as plt
    import math, pylab

    #extract number of reads, these are needed not only here but also for the bamqc
    summary_files = [param['fastqc_dir']+param['fastqc_stub'][idx]+'_fastqc/fastqc_data.txt' for idx in range(len(param['fastqc_stub']))]
    seq_length=[open(sum_file).readlines()[8].split('\t')[1].strip().rstrip() for sum_file in summary_files]
    for idx in range(len(seq_length)):
        if "-" in seq_length[idx]: 
            num_max=seq_length[idx].split('-')[1]
            seq_length[idx]=int(num_max)
        else:
            num=seq_length[idx]
            seq_length[idx]=int(num)

    #create plot 
    fig, ax = plt.subplots()
    fig.set_size_inches(3+len(param['fastqc_stub'])*0.4,6)
    index = np.arange(len(seq_length))
    bar_width = 0.8
    opacity = 0.4
    rects1 = plt.bar(index, seq_length, bar_width,
                     alpha=opacity,
                     color='b')
    plt.xlabel('Samples')
    plt.ylabel('Read length')
    plt.title('Read length across samples')
    ticks=param['fastqc_stub']
    plt.xticks(index + bar_width, ticks,rotation='vertical')
    plt.tight_layout()
    
    #put it into the report
    filename='report/fastqc/read_length.png'
    pylab.savefig(param['working_dir']+filename)
    param['report'].write('<img src="fastqc/read_length.png" alt="read length"><br><br>\n')


def report(param):  
    #assemble the full fastqc report
    param['report'].write('<center><br><h2>FastQC results</h2>')    
    copyFiles(param)
    if len(param['fastqc_stub'])>0:
        readRawfastqc(param)
        createOverviewTable(param)
        param['report'].write('<a href="fastqc/overview.txt">Table as tab delimited file</a><br><br><br>')
        plotNumberOfReads(param)
        plotGCContent(param)
        plotAgvLengthOfReads(param)
    else:
        param['report'].write('There were no results to show.')    


def init(param):
   module_helper.checkParameter(param,key='fastqc_exec',dType=str,checkFile=True)

def run_fastqc(filename,param):
    param['module_dir']=param['working_dir']+'results/fastqc/'
    if not os.path.exists(param['module_dir']):
        os.makedirs(param['module_dir'])
    call = param['fastqc_exec']+' '+" ".join(param[filename])+' -o '+param['module_dir']+' --extract'
    log_file=param['working_dir']+'results/fqc_log.txt'
    #if not os.path.exists(log_file):
    #    open(log_file,'a').close()
    param['file_handle'] = open(log_file,'a')
    param['file_handle'].write('CALL: '+call+'\n')
    output,error = subprocess.Popen(call.split() ,stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()
    param['file_handle'].write(error)
    param['file_handle'].write(output)
    if ('Analysis complete for' not in output):
       sys.exit()


if __name__ == "__main__":
    import subprocess, sys
    param=module_helper.initialize_module()

    #run fastqc
    run_fastqc('working_file',param)
    
    if not param['paired']:
        module_helper.wrapup_module(param)
    #calling it on the second fastq file if it is paired    
    else:
        run_fastqc('working_file2',param)
        module_helper.wrapup_module(param) 


