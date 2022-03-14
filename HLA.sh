#!/bin/bash

set -x

while test $# -gt 0; do
  case "$1" in
    -i)
      shift
      if test $# -gt 0; then
        SAMPLE=$1
      else
        echo 'Missing -i ID' && exit 1
      fi
      shift
      ;;
    -o)
      shift
      if test $# -gt 0; then
        PATHTOOUTPUT=$1
      else
        echo 'Missing -o Output directory' && exit 1
      fi
      shift
      ;;
    -b)
      shift
      if test $# -gt 0; then
        BAM=$1
      else
        echo 'Missing -b BAM file' && exit 1
      fi
      shift
      ;;
    --f1*)
      shift
      if test $# -gt 0; then
        FASTQ1=$1
      else
        echo 'Missing --f1 FASTQ1 file' && exit 1
      fi
      shift
      ;;
    --f2*)
      shift
      if test $# -gt 0; then
        FASTQ2=$1
      else
        echo 'Missing --f2 FASTQ2 file' && exit 1
      fi
      shift
      ;;
    *)
      echo "An invalid option has been entered" && exit 0
      ;;
  esac
done


echo $SAMPLE
echo $BAM
echo $PATHTOOUTPUT
echo $FASTQ1
echo $FASTQ2


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


  echo "" > ${HLA_ADDITIONAL}/${SAMPLE}'_Top_2_Results.txt'

  for gene in "A" "B" "C" "DQA1" "DQB1" "DRB1" "DPB1"; do
    perl ${REFDIR}/parse_result.pl ${REFDIR}/Allelelist_v2.txt \
      ${HLA_ADDITIONAL}/${SAMPLE}_results.txt | grep '^'${gene}'\*' | sort -k2 -n -r > ${TMPDIR}'/HLAVBSeq_HLA_'${gene}'.txt'
    head -n 2 ${TMPDIR}'/HLAVBSeq_HLA_'${gene}'.txt' > ${TMPDIR}'/HLAVBSeq_HLA_Final_'${gene}'_Top2.txt'
    cat ${TMPDIR}'/HLAVBSeq_HLA_Final_'${gene}'_Top2.txt' >> ${HLA_ADDITIONAL}/${SAMPLE}'_Top_2_Results.txt'
  done

  rm -r ${TMPDIR}

  # HLA tables
  python calling_hla_lookup_tables.py ${HLA} ${SAMPLE} ${HLA_ADDITIONAL}

} 2>&1 | tee ${HLA}/HLA.log
