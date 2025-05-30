#!/bin/bash

set -e

echo "创建结果目录 results/fastqc"
mkdir -p results/fastqc

echo "正在运行 FastQC..."
fastqc data/raw/SRR5168216_1.fastq.gz data/raw/SRR5168216_2.fastq.gz -o results/fastqc -t 4

echo "FastQC 分析完成，结果保存在 results/fastqc/"
