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

UNDER CONSTRUCTION

> Diana C. Vinueza-Espinosa, Daniel R. Cuesta-Aguirre, Assumpció Malgosa, Cristina Santos; MtDNA-CR Typing from Highly Degraded Skeletal Remains by Single-Multiplex Massively Parallel Sequencing, Unpublished
