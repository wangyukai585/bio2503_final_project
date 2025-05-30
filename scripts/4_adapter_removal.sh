#!/bin/bash

set -e

echo "创建清洗后的数据目录 data/clean"
mkdir -p data/clean

echo "正在执行 AdapterRemoval 去除接头序列和低质量 reads..."

AdapterRemoval \
  --file1 data/raw/SRR5168216_1.fastq.gz \
  --file2 data/raw/SRR5168216_2.fastq.gz \
  --output1 data/clean/SRR5168216_trimmed_1.fastq \
  --output2 data/clean/SRR5168216_trimmed_2.fastq \
  --trimns --trimqualities --minquality 20 --minlength 50 \
  --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
  --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

echo "AdapterRemoval 处理完成，输出文件保存在 data/clean/"
