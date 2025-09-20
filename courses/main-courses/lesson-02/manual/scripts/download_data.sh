#!/bin/bash

# 测序数据下载脚本
# 课程：高通量测序数据分析 - 测序数据质量控制与预处理
# 作者：王运生
# 日期：2025-01-20
# 用法：bash download_data.sh

set -e  # 遇到错误立即退出

echo "=========================================="
echo "下载测序数据用于质量控制练习"
echo "=========================================="

# 切换到数据目录
cd ~/ngs-analysis/lesson-02/data/

# 检查网络连接
echo "检查网络连接..."
if ping -c 1 google.com &> /dev/null; then
    echo "网络连接正常"
else
    echo "警告：网络连接可能有问题，将使用本地生成的测试数据"
    USE_LOCAL_DATA=true
fi

# 检查是否已安装SRA Toolkit
if command -v fastq-dump &> /dev/null; then
    echo "SRA Toolkit已安装"
    SRA_AVAILABLE=true
else
    echo "SRA Toolkit未安装，将使用其他方式获取数据"
    SRA_AVAILABLE=false
fi

# 函数：生成模拟测序数据
generate_mock_data() {
    echo "生成模拟测序数据..."
    
    python3 << 'EOF'
import gzip
import random
import sys

def generate_quality_score(position, total_length):
    """根据位置生成质量分数，模拟真实测序数据的质量分布"""
    if position < total_length * 0.1:
        # 前10%位置，质量稍低（测序起始不稳定）
        return random.randint(25, 35)
    elif position < total_length * 0.8:
        # 中间70%位置，质量较高
        return random.randint(30, 40)
    else:
        # 后20%位置，质量逐渐下降
        decline_factor = (position - total_length * 0.8) / (total_length * 0.2)
        max_quality = int(35 - decline_factor * 20)
        return random.randint(max(15, max_quality - 5), max_quality)

def add_sequencing_errors(sequence, error_rate=0.01):
    """在序列中添加测序错误"""
    bases = ['A', 'T', 'G', 'C']
    sequence_list = list(sequence)
    
    for i in range(len(sequence_list)):
        if random.random() < error_rate:
            # 随机选择一个不同的碱基
            current_base = sequence_list[i]
            other_bases = [b for b in bases if b != current_base]
            sequence_list[i] = random.choice(other_bases)
    
    return ''.join(sequence_list)

def add_adapter_contamination(sequence, adapter_seq="AGATCGGAAGAGC", contamination_rate=0.05):
    """模拟接头污染"""
    if random.random() < contamination_rate:
        # 在序列末端添加接头序列
        cut_point = random.randint(len(sequence) - 30, len(sequence) - 10)
        return sequence[:cut_point] + adapter_seq[:len(sequence) - cut_point]
    return sequence

def generate_fastq_read(read_id, length=150, gc_content=0.42):
    """生成一条FASTQ格式的序列"""
    # 根据GC含量生成序列
    sequence = []
    for _ in range(length):
        if random.random() < gc_content / 2:
            sequence.append('G')
        elif random.random() < gc_content / 2:
            sequence.append('C')
        elif random.random() < 0.5:
            sequence.append('A')
        else:
            sequence.append('T')
    
    sequence_str = ''.join(sequence)
    
    # 添加测序错误
    sequence_str = add_sequencing_errors(sequence_str)
    
    # 添加接头污染（少量序列）
    sequence_str = add_adapter_contamination(sequence_str)
    
    # 生成质量分数
    qualities = []
    for i in range(len(sequence_str)):
        q = generate_quality_score(i, length)
        qualities.append(chr(q + 33))  # Phred+33编码
    
    quality_string = ''.join(qualities)
    
    return f"@{read_id}\n{sequence_str}\n+\n{quality_string}\n"

def generate_duplicate_reads(base_read, num_duplicates=5):
    """生成重复序列"""
    duplicates = []
    lines = base_read.strip().split('\n')
    base_id = lines[0][1:]  # 去掉@符号
    
    for i in range(num_duplicates):
        new_id = f"@{base_id}_dup_{i+1}"
        duplicate_read = f"{new_id}\n{lines[1]}\n+\n{lines[3]}\n"
        duplicates.append(duplicate_read)
    
    return duplicates

# 生成R1文件
print("生成sample_R1.fastq.gz (包含各种质量问题)...")
with gzip.open('sample_R1.fastq.gz', 'wt') as f:
    duplicate_reads = []
    
    for i in range(50000):  # 生成50000条序列
        read = generate_fastq_read(f"SRR_mock_{i+1:06d}/1")
        f.write(read)
        
        # 每1000条序列生成一些重复序列
        if i % 1000 == 0 and i > 0:
            duplicates = generate_duplicate_reads(read, 3)
            for dup in duplicates:
                f.write(dup)

# 生成R2文件
print("生成sample_R2.fastq.gz (配对的第二端)...")
with gzip.open('sample_R2.fastq.gz', 'wt') as f:
    for i in range(50000):  # 生成50000条序列
        read = generate_fastq_read(f"SRR_mock_{i+1:06d}/2")
        f.write(read)
        
        # 对应的重复序列
        if i % 1000 == 0 and i > 0:
            duplicates = generate_duplicate_reads(read, 3)
            for dup in duplicates:
                f.write(dup)

print("模拟数据生成完成")
print("数据特征：")
print("- 序列长度：150bp")
print("- 序列数量：约53000条/文件")
print("- 包含质量下降、测序错误、接头污染、重复序列等问题")
EOF
}

# 函数：从SRA下载真实数据
download_sra_data() {
    echo "从SRA数据库下载真实测序数据..."
    
    # 使用一个较小的公共数据集进行演示
    SRR_ID="SRR1234567"  # 这是一个示例ID，实际使用时需要替换为真实的SRR ID
    
    echo "下载SRA数据：$SRR_ID"
    
    # 检查是否已经下载
    if [ -f "${SRR_ID}_1.fastq.gz" ] && [ -f "${SRR_ID}_2.fastq.gz" ]; then
        echo "数据已存在，跳过下载"
        return
    fi
    
    # 下载并转换为FASTQ格式
    if [[ "$SRA_AVAILABLE" == true ]]; then
        echo "使用fastq-dump下载数据..."
        fastq-dump --split-files --gzip --readids --read-filter pass \
                   --dumpbase --skip-technical --clip $SRR_ID
        
        # 重命名文件
        if [ -f "${SRR_ID}_1.fastq.gz" ]; then
            mv ${SRR_ID}_1.fastq.gz sample_R1.fastq.gz
            mv ${SRR_ID}_2.fastq.gz sample_R2.fastq.gz
            echo "SRA数据下载完成"
        else
            echo "SRA下载失败，使用模拟数据"
            generate_mock_data
        fi
    else
        echo "SRA Toolkit不可用，使用模拟数据"
        generate_mock_data
    fi
}

# 函数：从URL下载示例数据
download_example_data() {
    echo "从示例数据源下载..."
    
    # 示例：从公共FTP服务器下载小型测试数据
    # 这里使用模拟的URL，实际使用时需要替换为真实的数据源
    
    EXAMPLE_URL_R1="ftp://ftp.example.com/pub/data/sample_R1.fastq.gz"
    EXAMPLE_URL_R2="ftp://ftp.example.com/pub/data/sample_R2.fastq.gz"
    
    if command -v wget &> /dev/null; then
        echo "使用wget下载数据..."
        wget -O sample_R1.fastq.gz $EXAMPLE_URL_R1 || echo "wget下载失败"
        wget -O sample_R2.fastq.gz $EXAMPLE_URL_R2 || echo "wget下载失败"
    elif command -v curl &> /dev/null; then
        echo "使用curl下载数据..."
        curl -o sample_R1.fastq.gz $EXAMPLE_URL_R1 || echo "curl下载失败"
        curl -o sample_R2.fastq.gz $EXAMPLE_URL_R2 || echo "curl下载失败"
    else
        echo "wget和curl都不可用，使用模拟数据"
        generate_mock_data
    fi
    
    # 检查下载是否成功
    if [ ! -f "sample_R1.fastq.gz" ] || [ ! -f "sample_R2.fastq.gz" ]; then
        echo "下载失败，使用模拟数据"
        generate_mock_data
    fi
}

# 主要下载逻辑
if [ -f "sample_R1.fastq.gz" ] && [ -f "sample_R2.fastq.gz" ]; then
    echo "数据文件已存在，跳过下载"
    echo "如需重新下载，请删除现有文件："
    echo "rm sample_R1.fastq.gz sample_R2.fastq.gz"
else
    echo "选择数据源："
    echo "1. 生成模拟数据（推荐，快速且包含典型问题）"
    echo "2. 从SRA下载真实数据（需要网络连接和SRA Toolkit）"
    echo "3. 从示例服务器下载（需要网络连接）"
    
    # 自动选择最佳选项
    if [[ "${USE_LOCAL_DATA:-false}" == true ]]; then
        echo "网络连接问题，使用模拟数据"
        generate_mock_data
    elif [[ "$SRA_AVAILABLE" == true ]] && [[ "${USE_LOCAL_DATA:-false}" != true ]]; then
        echo "尝试从SRA下载数据..."
        # download_sra_data  # 注释掉，因为需要真实的SRR ID
        echo "SRA下载需要真实的SRR ID，使用模拟数据"
        generate_mock_data
    else
        echo "使用模拟数据（默认选择）"
        generate_mock_data
    fi
fi

# 验证数据文件
echo "验证下载的数据..."
if [ -f "sample_R1.fastq.gz" ] && [ -f "sample_R2.fastq.gz" ]; then
    echo "数据文件验证："
    echo "sample_R1.fastq.gz: $(ls -lh sample_R1.fastq.gz | awk '{print $5}')"
    echo "sample_R2.fastq.gz: $(ls -lh sample_R2.fastq.gz | awk '{print $5}')"
    
    # 检查文件内容
    echo "检查文件格式..."
    zcat sample_R1.fastq.gz | head -4
    echo "..."
    
    # 统计序列数量
    R1_COUNT=$(zcat sample_R1.fastq.gz | wc -l)
    R2_COUNT=$(zcat sample_R2.fastq.gz | wc -l)
    R1_READS=$((R1_COUNT / 4))
    R2_READS=$((R2_COUNT / 4))
    
    echo "序列统计："
    echo "R1文件：$R1_READS 条序列"
    echo "R2文件：$R2_READS 条序列"
    
    if [ $R1_READS -eq $R2_READS ]; then
        echo "✓ 双端数据配对正确"
    else
        echo "⚠ 警告：双端数据数量不匹配"
    fi
    
    # 生成MD5校验和
    echo "生成校验和..."
    md5sum sample_R1.fastq.gz sample_R2.fastq.gz > checksums.md5
    echo "校验和已保存到 checksums.md5"
    
else
    echo "错误：数据文件下载失败"
    exit 1
fi

# 创建数据说明文件
cat > data_info.txt << EOF
测序数据信息
============

文件列表：
- sample_R1.fastq.gz: 双端测序数据第一端
- sample_R2.fastq.gz: 双端测序数据第二端
- TruSeq3-PE.fa: Illumina TruSeq接头序列文件
- checksums.md5: 文件校验和

数据特征：
- 测序平台：模拟Illumina HiSeq数据
- 读长：150bp
- 序列数量：约50,000条配对reads
- 数据质量：包含典型的质量问题用于练习

质量问题包括：
1. 读长末端质量下降
2. 少量接头污染
3. 重复序列存在
4. 随机测序错误
5. GC含量轻微偏差

使用说明：
这些数据专门为质量控制练习设计，包含了真实测序数据中
常见的各种质量问题，适合用于学习FastQC分析和数据清洗。

生成时间：$(date)
EOF

echo "=========================================="
echo "数据下载完成！"
echo "数据位置：$(pwd)"
echo "数据说明：data_info.txt"
echo "可以开始进行质量控制分析了"
echo "=========================================="