#!/bin/bash

set -e

echo "正在创建数据目录..."
mkdir -p data/genome data/raw

echo "下载 E. coli K12 MG1655 参考基因组..."
cd data/genome
wget -O Ecoli_K12_MG1655.fa.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
gunzip -f Ecoli_K12_MG1655.fa.gz

echo "直接下载 SRR5168216 paired-end FastQ 文件..."
cd ../raw
wget -c https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR516/006/SRR5168216/SRR5168216_1.fastq.gz
wget -c https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR516/006/SRR5168216/SRR5168216_2.fastq.gz

echo "下载完成，跳过 SRA 解压，进入 FastQC 步骤。"
