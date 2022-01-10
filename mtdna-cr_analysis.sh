##Throughout the text, ${name} refers to sample's name 
##${PMDTOOLS_ROUTE} is the route where pmdtools.0.60.py (PMDtools script) is saved

#Sequencing quality control
fastqc -f fastq *fq.gz

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
  --in1 *R1*.fastq*.gz \
  --in2 *R2*.fastq*.gz
       
#25 bp (primers) are trimmed with another fastp
fastp --verbose \
  --trim_front1 25 --trim_tail1 25 \
  --trim_front2 25 --trim_tail2 25 \
  --json ${name}_fastp.json \
  --html ${name}_fastp.html \
  --report_title ${name} \
  --out1 ${name}_R1.fq.gz \
  --out2 ${name}_R2.fq.gz \
  --in1 ${name}_fastp_filtered_R1.fq.gz \
  --in2 ${name}_fastp_filtered_R2.fq.gz

#Both final trimmed fastQ files are aligned
bwa aln reference_mt.fa ${name}_R1.fq.gz > aln_R1.sai
bwa aln reference_mt.fa ${name}_R2.fq.gz > aln_R2.sai
bwa sampe reference_mt.fa aln_R1.sai aln_R2.sai ${name}_R1.fq.gz ${name}_R2.fq.gz | samtools view -Sbho ${name}_aln.bam

#Reds that did not aligned with the reference genome with a quality higher than 30 are discarded
samtools view -bh -F4 -q 30 ${name}_aln.bam > ${name}_map.bam

#BAM file is sorted
samtools sort -O bam -o ${name}.bam ${name}_map.bam

#BAM file is indexed
samtools index ${name}.bam

#Quality of BAM file is checked 
qualimap bamqc -bam ${name}.bam -c -gd hg19 -outdir ${name}_qualimap_bamqc  

#Reads with Post-Mortem Damage are selected
samtools view -h ${name}.bam | python ${PMDTOOLS_ROUTE}/pmdtools.0.60.py--threshold 0 --header | samtools view -Sb - > pmd_${name}.bam

#Polimorphisms are obteined with freebayes
freebayes -f reference_mt.fa -m 20 -q 30 -i -X -u -F 0.3 -C 10 ${name}.bam > ${name}.vcf 
