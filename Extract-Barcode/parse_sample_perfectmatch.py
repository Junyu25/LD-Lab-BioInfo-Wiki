#coding:utf-8
"""
@author: Yuxiang Tan
updated: 01/05/2020

summary: this script parses Illumina Paired-End data (FASTQ)

input: ./raw/*.fastq
Illumina HiSeq PE 100

output: ./parse/*_R1, ./parse/*_R2
format: see function readrecord

"""
"""
#   Copyright {2019} Yuxiang Tan
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

#This script will parse the fastq file with multiple indexes into seperated samples by the sample-index table 
#The output folder will be generated in the current working folder (it is pointed by the workfd in qsub)
#Usage rule: python /home/yxtan/qiime2_custom_scripts/parse_sample_perfectmatch -i input_fastq -t sample_index_table -f legnth_of_the_forward_primer -r length_of_the_reverse_primer

#This script starts from the joined fastq. 
#First, it will use the index table to genearte sample dictionary (the sample must follow this order: sample name;forward index;revesr index, tab-deliminated in three columns with no header.)
#Second, for each read, check the first Nbp for the forward index (N = the length of the index, and all the forward index must have the same length)
#Third, if the read past the second step, check the last Mbp for the reverse index (M = the length of the index, and all the reverse index must have the same length)
#Four, if both ends of the read have the indexes, assign this read to the fastq file of the corresponding sample. (trunc the index and primers at the mean time.)

#Library file requirement:
#biopython

updated: 04/05/2020


input: the fastq with multiple indexes

output: ./parse/*.fastq
*is the sample name corresponding to the index pair.
"""

import os
import sys
import glob
import string
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import time

#conver sample_index_table to two lists and a dictionary
def index_2_sample(sample_index_table):
  f_index = []
  r_index = []
  sample_index = {}
  File = open(sample_index_table)
  for line in File: 
    line = line.strip()
    Fields = line.split("\t")
    if len(Fields) != 3:
      print('The input sample table has incorret number of columns; Exit now.')
      exit(0)
    sample_name = Fields[0]
    f_i = Fields[1]
    r_i = Fields[2]
    if f_i not in f_index:
      f_index.append(f_i)
    if r_i not in r_index:
      r_index.append(r_i)
    dict_index = f_i+"_"+r_i
    if dict_index in sample_index:
      print('The input sample table has duplicated samples with the same index-pair, please check first; Exit now.')
      exit(0)
    sample_index[dict_index] = sample_name
  return (sample_index, f_index, r_index)

#parse FASTQ file using SeqIO#parse FASTQ file using SeqIO
def readrecord(input_fastq,sample_index, f_index, r_index, l_fp, l_rp):
  start_time = time.time()
  dict_record={}
  dict_record["r_index_not_mapped"]=[]
  dict_record["f_index_not_mapped"]=[]
  l_i_f = len(f_index[0])
  l_i_r = len(r_index[0])
  count_file=0
  count_used=0
  for record in SeqIO.parse(input_fastq, "fastq"):
    seq  = record.seq
    i_f = seq[0:l_i_f]
    i_r = seq[-l_i_r:].reverse_complement()
    #Second, for each read, check the first Nbp for the forward index (N = the length of the index, and all the forward index must have the same length)
    if i_f in f_index: #or i_r.reverse in r_index
      #Third, if the read past the second step, check the last Mbp for the reverse index (M = the length of the index, and all the reverse index must have the same length)
      if i_r in r_index: #or i_f.reverse in f_index
        record_in = record[(l_i_f+int(l_fp)):-(l_i_r+int(l_rp))]
        dict_index = str(i_f)+"_"+str(i_r)
        if dict_index in sample_index.keys():
          sample_name = sample_index[dict_index]
          count_used+=1
          #Four, if both ends of the read have the indexes, assign this read to the fastq file of the corresponding sample. (trunc the index and primers at the mean time.)
          if dict_record.get(sample_name)==None:
            dict_record[sample_name] = [record_in]
          else:
            dict_record[sample_name].append(record_in)
      else:
        dict_record["r_index_not_mapped"].append(record)
    elif i_f in r_index:
      if i_r in f_index:
        record_in = record[(l_i_f+int(l_fp)):-(l_i_r+int(l_rp))].reverse_complement()
        dict_index = str(i_r)+"_"+str(i_f)
        if dict_index in sample_index.keys():
          sample_name = sample_index[dict_index]
          count_used+=1
          #Four, if both ends of the read have the indexes, assign this read to the fastq file of the corresponding sample. (trunc the index and primers at the mean time.)
          if dict_record.get(sample_name)==None:
            dict_record[sample_name] = [record_in]
          else:
            dict_record[sample_name].append(record_in)
      else:
        dict_record["r_index_not_mapped"].append(record)
    else:
      dict_record["f_index_not_mapped"].append(record)
    count_file+=1

  outfile = open("reads_distr_of_samples.txt",'w')
  for sample_key in dict_record.keys():
    SeqIO.write(dict_record[sample_key],"parse/"+sample_key+".fastq", "fastq")
    outfile.write(sample_key+"\t"+str(len(dict_record[sample_key]))+"\n")
  outfile.close()

  #check files with emplty reads by using the keys in sample_index and dict_record
  outfile = open("empty_samples.txt",'w')
  for sample_v in sample_index.values():
    if sample_v not in dict_record.keys():
      outfile.write(sample_v+"\n")

  outfile.write("Total number of reads:"+str(count_file)+"\n"+"Total number of reads with appropriated index:"+str(count_used)+"\n"+"The percentage of useful reads:"+str(count_used*100/count_file)+"%\n")
  outfile.close()
  print("Parsing time: %s seconds " % (time.time() - start_time))
  return (count_file, count_used)


if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', '--input', dest='fq', type=str, required=True,
                      help="the joined fastq file with all multi-index reads together")
  parser.add_argument('-t', '--table', dest='tb', type=str, required=True, 
                      help="the tabular-table contains sample name;forward index;revesr index")
  parser.add_argument('-f', '--len_fp',dest='lf', type=str, required=True, 
                      help="the legnth_of_the_forward_primer")
  parser.add_argument('-r', '--len_rp', dest='lr', type=str, required=True, 
                      help="the length_of_the_reverse_primer")
  
  args = parser.parse_args()
  print('Usage example: python /home/yxtan/qiime2_custom_scripts/parse_sample_perfectmatch -i Undetermined_merged.fq -t sample_indexes.txt -f 20 -r 18')
  input_fastq = os.path.abspath(args.fq)
  sample_index_table = os.path.abspath(args.tb)
  l_fp = args.lf
  l_rp = args.lr

  #########check files and folders  
  if not os.path.isfile(input_fastq):
    print('Input fastq is not exist; Exit now.')
    exit(0)
  if not os.path.isfile(sample_index_table):
    print('The tabular-table is not exist; Exit now.')
    exit(0)
  if not os.path.exists("parse"):
    os.mkdir("parse")

  #step1:conver sample_index_table to two lists and a dictionary
  (sample_index, f_index, r_index) = index_2_sample(sample_index_table)

  #step2: parse the fastq
  (count_file, count_used) = readrecord(input_fastq,sample_index, f_index, r_index, l_fp, l_rp)