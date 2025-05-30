#!/bin/bash

set -e

echo "创建 VCF 过滤后输出目录 results/vcf"
mkdir -p results/vcf

echo "过滤低质量和低深度变异..."

bcftools filter -e 'INFO/DP<10 || QUAL<20' \
  results/vcf/SRR5168216_raw.vcf.gz -Oz -o results/vcf/SRR5168216_filtered.vcf.gz

echo "索引过滤后的 VCF..."
bcftools index results/vcf/SRR5168216_filtered.vcf.gz

echo "过滤完成，输出文件：results/vcf/SRR5168216_filtered.vcf.gz"
