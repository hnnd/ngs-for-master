#!/bin/bash
# 第1次课实验环境设置脚本
# 作者：王运生
# 日期：2025年
# 用法：bash setup.sh

set -e  # 遇到错误立即退出
set -u  # 使用未定义变量时报错

echo "开始设置第1次课实验环境..."
echo "================================"

# 检查Python环境
echo "检查Python环境..."
if ! command -v python3 &> /dev/null; then
    echo "错误：未找到Python3，请先安装Python"
    exit 1
fi

python_version=$(python3 --version | cut -d' ' -f2)
echo "Python版本: $python_version"

# 安装必要的Python包
echo "安装Python依赖包..."
pip3 install --user pandas matplotlib seaborn biopython numpy

# 创建目录结构
echo "创建工作目录结构..."
mkdir -p ~/ngs-analysis/lesson-01/{data,scripts,results,logs}

# 复制脚本文件
echo "复制分析脚本..."
cp *.py ~/ngs-analysis/lesson-01/scripts/

# 设置脚本执行权限
chmod +x ~/ngs-analysis/lesson-01/scripts/*.py

# 生成模拟数据（如果真实数据不可用）
echo "准备实验数据..."
cd ~/ngs-analysis/lesson-01/data

# 创建模拟FASTQ数据的Python脚本
cat > generate_sample_data.py << 'EOF'
#!/usr/bin/env python3
import random
import string

def generate_fastq_record(seq_id, length, platform='illumina'):
    """生成FASTQ记录"""
    # 生成随机DNA序列
    bases = 'ATCG'
    sequence = ''.join(random.choices(bases, k=length))
    
    # 生成质量分数
    if platform == 'illumina':
        # Illumina高质量
        quality_scores = [random.randint(25, 40) for _ in range(length)]
    elif platform == 'pacbio':
        # PacBio中等质量
        quality_scores = [random.randint(10, 30) for _ in range(length)]
    else:  # nanopore
        # Nanopore较低质量
        quality_scores = [random.randint(5, 25) for _ in range(length)]
    
    # 转换为ASCII字符
    quality_string = ''.join([chr(score + 33) for score in quality_scores])
    
    return f"@{seq_id}\n{sequence}\n+\n{quality_string}\n"

# 生成Illumina样本数据
print("生成Illumina样本数据...")
with open('illumina_sample.fastq', 'w') as f:
    for i in range(10000):
        length = random.randint(140, 160)  # Illumina典型读长
        record = generate_fastq_record(f"illumina_read_{i}", length, 'illumina')
        f.write(record)

# 生成PacBio样本数据
print("生成PacBio样本数据...")
with open('pacbio_sample.fastq', 'w') as f:
    for i in range(1000):
        length = random.randint(5000, 20000)  # PacBio长读长
        record = generate_fastq_record(f"pacbio_read_{i}", length, 'pacbio')
        f.write(record)

# 生成Nanopore样本数据
print("生成Nanopore样本数据...")
with open('nanopore_sample.fastq', 'w') as f:
    for i in range(500):
        length = random.randint(1000, 50000)  # Nanopore超长读长
        record = generate_fastq_record(f"nanopore_read_{i}", length, 'nanopore')
        f.write(record)

print("样本数据生成完成")
EOF

python3 generate_sample_data.py
rm generate_sample_data.py

echo "================================"
echo "环境设置完成！"
echo "工作目录: ~/ngs-analysis/lesson-01"
echo "可以开始实验了"