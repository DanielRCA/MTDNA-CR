# MTDNA-CR
Raw code to get polymorphisms from samples analysed by using the PowerSeq® CRM Nested System kit (Promega Corporation)

## Pre-requesites
A group of tools must be pre-installed. For each tool, version used by our group is shown in brackets, but another version might work as well:
- FastQC (v0.11.9)
- fastp (v0.23.2)
- BWA (v0.7.17)
- SAMtools (v1.16.1)
- QualiMap (v2.2.2a)
- PMDtools (v0.60)
- freebayes (v1.3.6)
- vcflib (v1.0.3)
- vcftools (v0.1.16)

## Pipeline Summary

- Sequencing quality control (FastQC)
- Duplicates removal (fastp)
- Adapter removal (fastp)
- Primers removal (fastp)
- Merge of overlapping Paired End reads (fastp)
- Reads mapping to reference (bwa aln, bwa sampe, samtools)
- Post-mapping processing (samtools)
- BAM quality control (QualiMap)
- Selection of reads with Post-Mortem Damage (PMDtools)
- Obtention of polimorphisms (freebayes + vcfallelicprimitives)

## Citations

If you use `DanielRCA/MTDNA-CR` for your analysis, please cite the `MTDNA-CR` paper as follows:

> Vinueza-Espinosa, DC, Cuesta-Aguirre, DR, Malgosa, A, Santos, C. Mitochondrial DNA control region typing from highly degraded skeletal remains by single-multiplex next-generation sequencing. Electrophoresis. 2023; 44: 1423–1434. https://doi.org/10.1002/elps.202200052

References of tools used:
- **FastQC** S. Andrews (2010), FASTQC. A quality control tool for high throughput sequence data, https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- **fastp** S. Chen, Y. Zhou, Y. Chen, J. Gu; fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, 1 September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560
- **BWA** H. Li, R. Durbin, Fast and accurate short read alignment with Burrows-Wheeler transform, Bioinformatics. 25 (2009) 1754–1760. https://doi.org/10.1093/bioinformatics/btp324
- **SAMtools** P. Danecek, JK. Bonfield, J. Liddle, J. Marshall, V. Ohan, MO. Pollard, A. Whitwham, T. Keane, SA. McCarthy, RM. Davies, H. Li (2021), Twelve years of SAMtools and BCFtools. GigaScience, Volume 10, Issue 2, February 2021, giab008, https://doi.org/10.1093/gigascience/giab008
- **QualiMap** K. Okonechnikov, A. Conesa, F. García-Alcalde (2016), Qualimap 2: advanced multi-sample quality control for high-throughput sequencing data.  Bioinformatics, Volume 32, Issue 2, 15 January 2016, Pages 292–294, https://doi.org/10.1093/bioinformatics/btv566
- **PMDtools** P. Skoglund, BH. Northoff, MV. Shunkov, A. Derevianko, S. Pääbo, J. Krause, M. Jakobsson (2014), Separating ancient DNA from modern contamination in a Siberian Neandertal. Proceedings of the National Academy of Sciences USA https://doi.org/10.1073/pnas.1318934111
- **freebayes** E. Garrison, G. Marth (2012), Haplotype-based variant detection from short-read sequencing. arXiv preprint arXiv:1207.3907 [q-bio.GN]
- **vcflib** E. Garrison, ZN. Kronenberg, ET. Dawson, BS Pedersen, P Prins (2022), A spectrum of free software tools for processing the VCF variant call format: vcflib, bio-vcf, cyvcf2, hts-nim and slivar. PLoS Comput Biol 18(5): e1009123, https://doi.org/10.1371/journal.pcbi.1009123
