import module_helper
import helper
import fastqc
import sys

    
def initialize_all(param):
    # this function calls the initialize functions of every module
    fastqc.init(param)


    
def run_all(param):
    #fastqc
    #helper.submit_job(param, 'fastqc.py',    input_files='fastq_files') 
    fastqc.run_fastqc('fastq_files',param)
    module_helper.wrapup_module(param)
    
    
def report_all(param):
    fastqc.report(param)
        
    
    