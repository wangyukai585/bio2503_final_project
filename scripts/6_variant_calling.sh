#!/bin/bash

set -e

# 创建输出目录
echo "创建 VCF 目录 results/vcf"
mkdir -p results/vcf

# Step 1: 构建参考基因组 fa 的索引（如不存在）
if [ ! -f data/genome/Ecoli_K12_MG1655.fa.fai ]; then
  echo "构建参考基因组索引..."
  samtools faidx data/genome/Ecoli_K12_MG1655.fa
fi

# Step 2: 生成 BCF 文件（调用变异）
echo "生成 BCF 文件..."
bcftools mpileup -Ou -f data/genome/Ecoli_K12_MG1655.fa results/bam/SRR5168216.sorted.bam |
  bcftools call -mv -Ob -o results/vcf/SRR5168216_raw.bcf

# Step 3: 转换为 VCF 格式
echo "转换 BCF 为 VCF 格式..."
bcftools view results/vcf/SRR5168216_raw.bcf -Oz -o results/vcf/SRR5168216_raw.vcf.gz

# Step 4: 创建 VCF 索引
echo "索引 VCF 文件..."
bcftools index results/vcf/SRR5168216_raw.vcf.gz

echo "变异检测完成，结果保存为 VCF 文件：results/vcf/SRR5168216_raw.vcf.gz"
