#!/bin/bash

set -x

### folders
REFDIR=${PWD}/Configuration
OUTDIR=${PATHTOOUTPUT}/${SAMPLE}
TMPDIR=${PATHTOOUTPUT}/${SAMPLE}/${SAMPLE}_TEMP
HLA=${OUTDIR}/${SAMPLE}_HLA_REPORT
HLA_ADDITIONAL=${OUTDIR}/${SAMPLE}_HLA_SUPPLEMENTARY

mkdir -p ${OUTDIR}
mkdir -p ${TMPDIR}
mkdir -p ${HLA}
mkdir -p ${HLA_ADDITIONAL}


##### HLAVBSeq

{
  samtools view ${BAM} chr6:28500120-33490577 | cut -f1 | sort | uniq > ${TMPDIR}/${SAMPLE}_readnames.list

  zcat < ${FASTQ1} | seqtk subseq - ${TMPDIR}/${SAMPLE}_readnames.list > \
    ${TMPDIR}/${SAMPLE}_Extracted_R1.fastq

  zcat < ${FASTQ2} | seqtk subseq - ${TMPDIR}/${SAMPLE}_readnames.list > \
    ${TMPDIR}/${SAMPLE}_Extracted_R2.fastq


  bwa index ${REFDIR}/hla_all_v2.fasta

  bwa  mem -t 20 -P -L 10000 -a  ${REFDIR}/hla_all_v2.fasta \
    ${TMPDIR}/${SAMPLE}_Extracted_R1.fastq ${TMPDIR}/${SAMPLE}_Extracted_R2.fastq > \
    ${TMPDIR}/${SAMPLE}_extracted.sam

  # final result
  java -jar ${REFDIR}/HLAVBSeq.jar ${REFDIR}/hla_all_v2.fasta \
    ${TMPDIR}/${SAMPLE}_extracted.sam ${HLA_ADDITIONAL}/${SAMPLE}'_results.txt' --alpha_zero 0.01 --is_paired

  # accounts for homozygous alleles
  python ${REFDIR}/call_hla_digits.py -v ${HLA_ADDITIONAL}/${SAMPLE}'_results.txt' -a \
    ${REFDIR}/Allelelist_v2.txt -r 90 -d 8 --ispaired > ${HLA}/${SAMPLE}'_report.d8.txt'

  # allele and read depth
  perl ${REFDIR}/parse_result.pl ${REFDIR}/Allelelist_v2.txt \
    ${HLA_ADDITIONAL}/${SAMPLE}'_results.txt' > ${HLA_ADDITIONAL}/${SAMPLE}'_Allele_Depth.txt'


} 2>&1 | tee ${HLA}/HLA.log
