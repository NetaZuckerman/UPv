#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate gatk4
bam_path=$1
vcf_path=$2
reference=$3
picard="java -jar /mnt/project2/PersonalFolders/hagar/picard/build/libs/picard.jar"

samtools faidx "$reference"
 
$picard CreateSequenceDictionary -R $reference -O "$(dirname "${reference}")"/`basename "$reference" .fasta`.dict
#add reading group manually. it is required by gatk and we normally didnt generate it so i add it menually
for file in "$bam_path"*mapped.sorted.bam; do
  #echo $bam
  RG_file="$bam_path"`basename "$file" .mapped.sorted.bam`.RG.bam
  #echo $RG_file
  $picard AddOrReplaceReadGroups I= "$file" O= "$RG_file" SORT_ORDER=coordinate RGID=foo RGLB=bar RGPL=illumina RGSM=Sample1 RGPU=unit1
  samtools index "$RG_file"
  gatk --java-options -Xmx7g HaplotypeCaller -R  "$reference" -I "$RG_file" -O "$vcf_path"/`basename "$RG_file" .RG.bam`.vcf
  gatk VariantsToTable -V "$vcf_path"/`basename "$RG_file" .RG.bam`.vcf -F CHROM -F POS -F TYPE -GF AD -O "$vcf_path"/`basename "$RG_file" .RG.bam`.table
done

conda deactivate