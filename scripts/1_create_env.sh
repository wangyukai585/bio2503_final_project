#!/bin/bash

set -e  # 遇到错误立即退出脚本
ENV_NAME="bioinfo_course_project"

echo "🚧 正在准备 Conda 环境: $ENV_NAME"

# 添加 Bioconda 和 Conda-forge 频道（只需添加一次）
if ! grep -q "bioconda" ~/.condarc 2>/dev/null; then
  echo "首次添加 Conda 频道"
  conda config --add channels defaults
  conda config --add channels bioconda
  conda config --add channels conda-forge
  conda config --set channel_priority strict
else
  echo "Conda 频道已配置，无需重复添加"
fi

# 检查环境是否已存在
if conda info --envs | grep -q "^$ENV_NAME"; then
  echo "环境 $ENV_NAME 已存在，先删除..."
  conda remove -y --name $ENV_NAME --all
fi

# 创建环境
echo "开始创建环境 $ENV_NAME..."
conda create -y -n $ENV_NAME fastqc bwa samtools seqtk adapterremoval picard sra-tools bcftools

echo "Conda 环境创建成功，请运行以下命令激活："
echo ""
echo "    conda activate $ENV_NAME"
