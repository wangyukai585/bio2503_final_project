#!/bin/bash

set -e

# 创建输出目录
echo "创建比对结果目录 results/bam"
mkdir -p results/bam

# 构建参考基因组索引（只需一次）
if [ ! -f data/genome/Ecoli_K12_MG1655.fa.bwt ]; then
  echo "构建 BWA 索引..."
  bwa index data/genome/Ecoli_K12_MG1655.fa
else
  echo "索引已存在，跳过"
fi

# 执行 BWA 比对
echo "正在进行 BWA 比对..."
bwa mem -t 4 data/genome/Ecoli_K12_MG1655.fa \
  data/clean/SRR5168216_trimmed_1.fastq \
  data/clean/SRR5168216_trimmed_2.fastq |
  samtools view -Sb - > results/bam/SRR5168216.bam

# BAM 文件排序
echo "排序 BAM 文件..."
samtools sort -@ 4 -o results/bam/SRR5168216.sorted.bam results/bam/SRR5168216.bam

# 创建 BAM 索引
echo "索引 BAM 文件..."
samtools index results/bam/SRR5168216.sorted.bam

echo "比对完成，已生成排序并索引的 BAM 文件"
