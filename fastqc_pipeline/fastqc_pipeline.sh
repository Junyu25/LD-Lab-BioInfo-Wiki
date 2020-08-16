#!/bin/bash
#This script will run the fastqc pipeline.
#This script starts from a input folder, and the path of a 
#Library file requirement: python2.7
if [ $# -ne 2 ]
then
    echo ""
    echo "Usage: bash fastqc_pipeline.sh input_fd script_path"
    echo "Example:bash /home/yxtan/fastqc_pipeline/fastqc_pipeline.sh /mnt/data2/LD_lab/yxtan/StrainAnal_tools/MockData/Xtrain_bench/Fprau3/final_fq/ /home/yxtan/fastqc_pipeline/"
    echo ""
    echo "input_fd - the folder with splitted fastqs."
    echo "script_path - the path of scripts."
    echo "Note: path of fastqc (/home/hrzhang/software/FastQC/fastqc) is harcoded in the parameter_file_example.txt in step2, modify it if necessary."
    exit 1
fi

#set para
input_fd=$1 
script_path=$2 

#Check folders
if [ ! -d $input_fd ] 
then
  echo ""
  echo "Warning: The directory $input_fd does not exist in fastqc_pipeline.sh, exit."
  echo ""
  exit 1
fi

#Check folders
if [ ! -d $script_path ] 
then
  echo ""
  echo "Warning: The directory $script_path does not exist in fastqc_pipeline.sh, exit."
  echo ""
  exit 1
fi

#Step 1:run fq_sample_table_generate.py to generate sample_table.txt in the working folder
python $script_path"/fq_sample_table_generate.py" -i $input_fd -F .fq -e False -N . -o $input_fd'/sample_table.txt'

#check out file from previous step
if [ ! -s $input_fd'/sample_table.txt' ] 
then
  echo ""
  echo "Warning: The file "$input_fd"/sample_table.txt does not exist in fastqc_pipeline.sh; There must be no fq files in the folder, exit."
  echo ""
  exit 1
fi

#Step 2:generate parameter_file.txt
cp $script_path"/parameter_file_example.txt" $input_fd"/parameter_file.txt"
#replace the paramters
sed -i -e "s|\/home\/yxtan\/final_fq_test\/|$input_fd\/|g" $input_fd"/parameter_file.txt"
sed -i -e "s|\/home\/yxtan\/fastqc_pipeline\/|$script_path\/|g" $input_fd"/parameter_file.txt"


#step 3: run fastqc and summary
python2.7 $script_path"/fastqc_pipeline.py" -p $input_fd"/parameter_file.txt"
#cp the Icons to the target folder
cp -rf $script_path"/Icons" $input_fd"/report/"