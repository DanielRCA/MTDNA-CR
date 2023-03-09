##Throughout the text, ${name} refers to sample's name 
##${PMDTOOLS_ROUTE} is the route where pmdtools.0.60.py (PMDtools script) is saved

#Sequencing quality control
fastqc -f fastq *fastq.gz

#Duplicates are removed with fastp
fastp --verbose \
        --dedup \
        --disable_length_filtering \
        --disable_quality_filtering \
        --json ${name}_dedup_temp.json \
        --thread 8 \
        --out1 ${name}_fastp_dedup_R1_temp.fq.gz \
        --out2 ${name}_fastp_dedup_R2_temp.fq.gz \
        --in1 *R1*.fastq.gz \
        --in2 *R2*.fastq.gz

#Adapters are trimmed with fastp
fastp --verbose \
  --detect_adapter_for_pe \
  --trim_poly_x \
  --length_required 30 \
  --low_complexity_filter --complexity_threshold 55 \
  --json ${name}_int_fastp.json \
  --html ${name}_int_fastp.html \
  --report_title ${name}_int  \
  --thread 8 \
  --out1 ${name}_fastp_filtered_R1.fq.gz \
  --out2 ${name}_fastp_filtered_R2.fq.gz \
  --in1 ${name}_fastp_dedup_R1_temp.fq.gz \
  --in2 ${name}_fastp_dedup_R2_temp.fq.gz 
       
#25 bp (primers) are trimmed with another fastp
fastp --verbose \
  --trim_front1 25 --trim_tail1 25 \
  --trim_front2 25 --trim_tail2 25 \
  --json ${name}_fastp.json \
  --html ${name}_fastp.html \
  --report_title ${name} \
  --out1 ${name}_nonprimer_R1_temp.fq.gz \
  --out2 ${name}_nonprimer_R2_temp.fq.gz \
  --in1 ${name}_fastp_filtered_R1.fq.gz \
  --in2 ${name}_fastp_filtered_R2.fq.gz

#Paired End FastQ files are merged
fastp --verbose \
      --cut_tail --cut_tail_mean_quality 20 \
      --trim_poly_x \
      --merge --overlap_len_require 10 \
      --correction \
      --html ${name}_temp.html \
      --json ${name}_temp.json \
      --merged_out ${d}_M.fq.gz \
      --out1 ${name}_R1.fq.gz \
      --out2 ${name}_R2.fq.gz \
      --in1 ${name}_nonprimer_R1_temp.fq.gz \
      --in2 ${name}_nonprimer_R2_temp.fq.gz

#Trimmed FastQ files quality control
fastqc -f fastq *.fq.gz

#Trimmed FastQ files are aligned with the refernece genome and only those with a quality higher than 30 are kept
bwa aln -l 1024 -n 0.01 -o 2 reference_mt.fa ${name}_R1.fq.gz > aln_R1.sai
bwa aln -l 1024 -n 0.01 -o 2 reference_mt.fa ${name}_R2.fq.gz > aln_R2.sai
bwa aln -l 1024 -n 0.01 -o 2 reference_mt.fa ${name}_M.fq.gz > aln_M.sai
bwa sampe reference_mt.fa aln_R1.sai aln_R2.sai ${name}_R1.fq.gz ${name}_R2.fq.gz | samtools view -Sbh -F4 -q30 | samtools sort -O bam -o tem_${name}_NM_map.bam
bwa samse reference_mt.fa aln_M.sai ${name}_M.fq.gz | samtools view -Sbh -F4 -q30 | samtools sort -O bam -o tem_${name}_M_map.bam

#Aligned bam files (merged and non-merged Paired End reads) are merged
samtools merge -o tem_${name}_map_unsort.bam tem_${name}_NM_map.bam tem_${name}_M_map.bam

#BAM file is sorted
samtools sort -O bam -o ${name}.bam tem_${name}_map_unsort.bam

#BAM file is indexed
samtools index ${name}.bam

#Quality of BAM file is checked 
qualimap bamqc -bam ${name}.bam -c -gd hg19 -outdir ${name}_qualimap_bamqc  

#Reads with Post-Mortem Damage are selected
samtools view -h ${name}.bam | python ${PMDTOOLS_ROUTE}/pmdtools.0.60.py--threshold 1 --header | samtools view -Sb - > pmd_${name}.bam

#Polimorphisms are obteined with freebayes
freebayes -f reference_mt.fa -i -X -F 0.1 -C 3 --min-coverage 30 ${name}.bam | vcfallelicprimitives -kg > ${name}.vcf


