# MTDNA-CR
Raw code to get polymorphisms from samples analysed by using the PowerSeq® CRM Nested System kit (Promega Corporation)

## Pre-requesites
A group of tools must be pre-installed. For each tool, version used by our group is shown in brackets, but another version might work as well:
- FastQC (v0.11.9)
- fastp (v0.23.2)
- BWA (v0.7.17)
- SAMtools (v1.13)
- QualiMap (v2.2.2a)
-  PMDtools (v0.60)
-  Freebayes (v0.9.21.7)

## Pipeline Summary

- Sequencing quality control (FastQC)
- Adapter removal (fastp)
- Primers removal (fastp)
- Reads mapping to reference (bwa aln, bwa sampe, samtools)
- Post-mapping processing (samtools)
- BAM quality control (QualiMap)
- Selection of reads with Post-Mortem Damage (PMDtools)
- Obtention of polimorphisms (freebayes)

## Citations

This pipeline is not already published (it is under submission)

> Diana C. Vinueza-Espinosa, Daniel R. Cuesta-Aguirre, Assumpció Malgosa, Cristina Santos; MtDNA-CR Typing from Highly Degraded Skeletal Remains by Single-Multiplex Massively Parallel Sequencing, Unpublished

References of tools used:
- **FastQC** S. Andrews (2010), FASTQC. A quality control tool for high throughput sequence data, https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- **fastp** S. Chen, Y. Zhou, Y. Chen, J. Gu; fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, 1 September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560
- **BWA** H. Li, R. Durbin, Fast and accurate short read alignment with Burrows-Wheeler transform, Bioinformatics. 25 (2009) 1754–1760. https://doi.org/10.1093/bioinformatics/btp324
- **SAMtools** P. Danecek, JK. Bonfield, J. Liddle, J. Marshall, V. Ohan, MO. Pollard, A. Whitwham, T. Keane, SA. McCarthy, RM. Davies, H. Li (2021), Twelve years of SAMtools and BCFtools. GigaScience, Volume 10, Issue 2, February 2021, giab008, https://doi.org/10.1093/gigascience/giab008
- **QualiMap** K. Okonechnikov, A. Conesa, F. García-Alcalde (2016), Qualimap 2: advanced multi-sample quality control for high-throughput sequencing data.  Bioinformatics, Volume 32, Issue 2, 15 January 2016, Pages 292–294, https://doi.org/10.1093/bioinformatics/btv566
- **PMDtools** P. Skoglund, BH. Northoff, MV. Shunkov, A. Derevianko, S. Pääbo, J. Krause, M. Jakobsson (2014), Separating ancient DNA from modern contamination in a Siberian Neandertal. Proceedings of the National Academy of Sciences USA doi:10.1073/pnas.1318934111
- **Freebayes** E. Garrison, G. Marth (2012), Haplotype-based variant detection from short-read sequencing. arXiv preprint arXiv:1207.3907 [q-bio.GN]
