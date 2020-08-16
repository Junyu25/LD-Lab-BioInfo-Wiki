#coding:utf-8
#!/bin/bash

#   Copyright {2015} Yuxiang Tan
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

#This script will run bowtie2 on target folder with fastqs

#This script starts from bowtie2 on each fastq and ref.fa; then the sam output will be sorted into bam and index; finally, bedtools multicov will be used to extract coverage.
#Library file requirement:
#bowtie2\samtools\bedtools should be in PATH

#updated: 04/05/2020

#input1: the folder of fastqs
#input2: The bed file (must be customized into strains in specific format) and the fasta with bowtie2 index of the refernce strains


#output_path: current folder
#final output: count_matrix.txt
#intermediate outputs: intermediate folder with sam\bam\bai of each sample.



if [ $# -ne 3 ]
then
  echo ""
    echo "Usage: Bowtie2align_and_feature_extract.sh input_folder input_ref.fa input_ref.bed "
    echo "Example: bash Bowtie2align_and_feature_extract.sh /thinker/dstore/gene2/r3data/ytan/microbiome_proj/Shi_QIPING_Diabetes/grouped.fasta /thinker/dstore/gene2/r3data/ytan/microbiome_proj/Shi_QIPING_Diabetes/Map.txt /mnt/fastdata/ytan/microbio/Shi_QIPING_Diabetes/ /mnt/dstore/gene2/r3data/ytan/microbiome_proj/program/"
    echo ""
    echo "folder of input_fastqs - The folder with all targeted fastq files."
    echo "input_ref.fa - The fasta file of the refernce region of each strain; the name of strains should be the same as the ones in the bed file"
    echo "input_ref.bed - The bed file of the refernce region of each strain; the name of strains should be the same as the ones in the fasta file. The format of lines should be: strain_name\t0\tlenth\tstrain_name\n"
    exit 1
fi

#name the parameters
fd_fq=$1
ref_fa=$2
ref_bed=$3
fd_imd="intermediate_fd"
count_mx="count_matrix.txt"
#Check folders
if [ ! -d $fd_fq ] 
then
  echo ""
  echo "Warning: The directory $fd_fq does not exist in Bowtie2align_and_feature_extract.sh, exit."
  echo ""
  exit 1
fi

if [ ! -d $fd_imd ] 
then
  echo ""
  echo "Warning: The directory $fd_imd does not exist in Bowtie2align_and_feature_extract.sh, generate it."
  echo ""
  mkdir -p "$fd_imd"
fi


#check files
if [ ! -s $ref_fa ] 
then
  echo ""
  echo "Warning: The file $ref_fa does not exist in Bowtie2align_and_feature_extract.sh, exit."
  echo ""
  exit 1
fi

if [ ! -s $ref_bed ] 
then
  echo ""
  echo "Warning: The file $ref_bed does not exist in Bowtie2align_and_feature_extract.sh, exit."
  echo ""
  exit 1
fi

bowtie2_log="bowtie2.log"
echo "Step1: bowtie2 align and format conversion"
echo `date`
#generate bowtie2 index
bowtie2-build $ref_fa $ref_fa
#In a for loop, generate header, bowtie2 for fastqs and generate bam 
bam_line=""


sample_list=$(ls $fd_fq/*.fastq)
array=(${sample_list// / })
fq_num=${#array[@]}
echo "The total number of fastqs are: "$fq_num
group_size=500

group_count=1
if [ ! -d $fd_imd"/"$group_count ] 
then
  echo ""
  echo "Warning: The directory $fd_imd"/"$group_count does not exist in Bowtie2align_and_feature_extract.sh, generate it."
  echo ""
  mkdir -p $fd_imd"/"$group_count
fi

fastq_count=1
echo -e "#StrainName\tbegin\tEnd\tStrainName\c" > $group_count"_"$count_mx
echo "Group number "$group_count
for i in ${array[*]}
do
  if [ $fastq_count -gt $(( $group_count*$group_size)) ]
  then
    echo -e "\n\c" >> $group_count"_"$count_mx
    group_count=$(( $group_count+1))
    echo -e "#StrainName\tbegin\tEnd\tStrainName\c" > $group_count"_"$count_mx
    if [ ! -d $fd_imd"/"$group_count ] 
    then
      echo ""
      echo "Warning: The directory $fd_imd"/"$group_count does not exist in Bowtie2align_and_feature_extract.sh, generate it."
      echo ""
      mkdir -p $fd_imd"/"$group_count
    fi
    echo "Group number "$group_count
  fi
  fastq_name=${i##*/} 
  #bam_line=$bam_line$fd_imd"/"$sample_name".bam " 
  sample_name=${fastq_name%%.*} 
  echo -e "\t"$sample_name"\c" >> $group_count"_"$count_mx 
  echo $i >> $bowtie2_log
  #can clean up the low qual alignment by bowtie or samtools.
  bowtie2 -x $ref_fa -U $i -S $fd_imd"/"$group_count"/"$sample_name".sam" --very-sensitive 1>>$bowtie2_log 2>>$bowtie2_log
  #filter by MAPQ30, not 42, but it's wrong, because the MAPQ various among samples.
  #samtools view -hSq 30 $fd_imd"/"$group_count"/"$sample_name".sam" -o $fd_imd"/"$group_count"/"$sample_name"_filtered.sam"
  #filter by NM:i:0 and 1 should be better
      samtools view -HS $fd_imd"/"$group_count"/"$sample_name".sam" -o $fd_imd"/"$group_count"/"$sample_name"_filtered.sam"
      grep "NM:i:0" $fd_imd"/"$group_count"/"$sample_name".sam" >> $fd_imd"/"$group_count"/"$sample_name"_filtered.sam"
      grep "NM:i:1"$'\t' $fd_imd"/"$group_count"/"$sample_name".sam" >> $fd_imd"/"$group_count"/"$sample_name"_filtered.sam"
    samtools sort $fd_imd"/"$group_count"/"$sample_name"_filtered.sam" -o $fd_imd"/"$group_count"/"$sample_name".bam"
  samtools index $fd_imd"/"$group_count"/"$sample_name".bam"
  fastq_count=$(( $fastq_count+1)) 
done

if [ $(( $fq_num-($group_count-1)*$group_size)) -lt $group_size ]
  then
      echo -e "\n\c" >> $group_count"_"$count_mx
fi


#echo $bam_line
group_num=$group_count

echo "Step2: bedtools multicov to get raw count"
echo `date`
group_count=1
while [ $group_count -le $group_num ]
do
  echo "Group number "$group_count
  bedtools multicov -bams $(ls $fd_imd/$group_count/*.bam) -bed $ref_bed >> $group_count"_"$count_mx
  group_count=$(( $group_count+1))
done

echo "Step3: merge the matrixs into one if necessary"
group_count=1
cut $group_count"_"$count_mx -f 4- > $group_count"_"Tem.txt
temp_line=$group_count"_"Tem.txt
if [ $group_num -gt 1 ]
then
  group_count=2
  while [ $group_count -le $group_num ]
  do
    cut "${group_count}"_"${count_mx}" -f 5- > "${group_count}"_Tem.txt
    temp_line=$temp_line" "$group_count"_"Tem.txt   
    group_count=$(( $group_count+1))
  done
fi
echo $temp_line
paste -d "\t" $temp_line > $count_mx
rm $temp_line