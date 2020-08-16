#coding:utf-8
from __future__ import print_function
import glob,os,re,argparse


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

#This script will scan the folder and generate 3 column paired/single end sample table automatically, the output will be in the current working folder.

#This script starts from raw fqs in a single folder
#Assumption, the files in the folder must have a unique recognizer for forward or reverse end in fqs,default is "R1", otherwise, you can specify.
#Usage example:python /home/yxtan/YT_scripts/fq_sample_table_generate.py -i abs_path_raw_fq_folder
"""

def pair_pro(rd_dir,sp1_str, sp2_str, spN_str, OutN):
    manifest_f = open(OutN,"w")
    manifest_f.write('#SampleID\tforward-absolute-filepath\treverse-absolute-filepath\n')
    file_names = os.listdir(rd_dir)
    cwd = os.getcwd()
    
    for files in file_names:
    	if files.find(sp1_str)>0:
    		files_part = files.strip().split(sp1_str,1)
    		p1=files_part[0]
    		p2=files_part[1]
    		files_head = files.split(".fq")[0].split(spN_str)[0]
    		out_str = "%s\t%s/%s%s%s\t%s/%s%s%s\n" % (files_head, rd_dir, p1,sp1_str, p2, rd_dir, p1,sp2_str, p2)
    		manifest_f.write(out_str)
    manifest_f.close()


def single_pro(rd_dir,sp1_str, spN_str, OutN):
    manifest_f = open(OutN,"w")
    manifest_f.write('#absolute-filepath\tSampleID\n')
    file_names = os.listdir(rd_dir)
    cwd = os.getcwd()
    for files in file_names:
    	if files.find(sp1_str)>0:
    		files_part = files.split(sp1_str,1)
    		p1=files_part[0]
    		p2=files_part[1]
    		files_head = files.split(".fq")[0].split(spN_str)[0]
    		out_str = "%s/%s%s%s\t%s\n" % (rd_dir, p1,sp1_str, p2, files_head)
    		manifest_f.write(out_str)
    manifest_f.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='rd', type=str, required=True,
                        help="the path of folder with raw fqs")
    parser.add_argument('-e', '--pair', dest='pair', type=str, required=False, default='True',
                        help="Is it pair-end seq data? Default is 'True'; Any other strings will be considered False")
    parser.add_argument('-F', '--sepF', dest='sp1', type=str, required=False, default='R1',
                        help="It is the seperator to recognize the forward info, default='R1'; for single end, should use 'fq' or similar.")
    parser.add_argument('-R', '--sepR', dest='sp2', type=str, required=False, default='R2',
                        help="It is the seperator to recognize the reverse info")
    parser.add_argument('-N', '--sepN', dest='spN', type=str, required=False, default='_',
                        help="It is the seperator to recognize the sample name")
    parser.add_argument('-o', '--output', dest='OutN', type=str, required=False, default='sample_table.txt',
                        help="The name of outputfile")
    args = parser.parse_args()
    print('Usage example with minimum parameters:   -i abs_path_raw_fq_folder')
    rd_dir = os.path.abspath(args.rd)
    sp1_str = args.sp1
    sp2_str = args.sp2
    spN_str = args.spN
    pair_end = args.pair
    OutN = args.OutN
    print(rd_dir)
    print(sp1_str)
    print(sp2_str)
    print(spN_str)
    print(OutN)
    #需要检查输入的参数是否正确，主要是路径是否存在    
    if not os.path.isdir(rd_dir):
        print('Input sample folder is not exist; Exit now.')
        exit(0)
    
    """
    :rd_dir: the dir where the raw data are '/'.
    :return: None
    """    
    if pair_end == 'True':
        pair_pro(rd_dir,sp1_str, sp2_str, spN_str, OutN)
    else:
        single_pro(rd_dir,sp1_str, spN_str, OutN)

