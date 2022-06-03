# HLA

This repository outlines an end-to-end pipeline using HLAVBSeq and an in-house Python script to assess disease risks based on HLA genotyping.

## Built With
- HLAVBSeq
- Python

## Getting Started
### Installation
```
git clone https://github.com/jlelabs/HLA
```
### Configuration Files
Download the following from HLAVBSeq http://nagasakilab.csml.org/hla/ and transfer them into your Configuration folder within the HLA repository. 
- HLAVBSeq.jar
- parse_result.pl
- call_hla_digits.py
- hla_all_v2.fasta
- Allelelist_v2.txt

### Usage

Download and configure NA12878 sample for ancestry estimation
```
./NA12878.sh
```

Run main HLA script
```
./HLA.sh -i ${SAMPLEID} -v ${VCF} -b ${BAM} --f1 ${FASTQ1} --f2 ${FASTQ2} -o ${OUTPUT_DIRECTORY}
```
