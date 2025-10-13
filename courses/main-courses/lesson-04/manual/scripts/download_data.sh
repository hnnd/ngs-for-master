#!/bin/bash

# 变异检测数据准备脚本
# 课程：高通量测序数据分析 - 变异检测与基因分型
# 作者：王运生
# 日期：2025-01-20
# 用法：bash download_data.sh

set -e  # 遇到错误立即退出

echo "=========================================="
echo "准备变异检测实验数据"
echo "=========================================="

# 切换到数据目录
WORK_DIR="$HOME/ngs-analysis/lesson-04"
mkdir -p "$WORK_DIR"
cd "$WORK_DIR"
mkdir -p data reference results logs

cd data

# 检查必要软件
check_software() {
    if command -v $1 &> /dev/null; then
        echo "✓ $1 已安装"
        return 0
    else
        echo "✗ $1 未安装"
        return 1
    fi
}

echo "检查必要软件..."
check_software "samtools" || echo "请安装 samtools: conda install -c bioconda samtools"
check_software "python3" || echo "请安装 Python3"

# ============================================
# 1. 准备参考基因组
# ============================================
echo ""
echo "步骤 1: 准备参考基因组..."

if [ -f "chr22_fragment.fa" ]; then
    echo "参考基因组已存在，跳过下载"
else
    echo "下载人类22号染色体片段..."

    # 选项1: 从Ensembl下载 (推荐)
    if command -v wget &> /dev/null; then
        echo "使用wget从Ensembl下载..."
        wget -q --show-progress -O chr22.fa.gz \
            "http://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz" \
            2>&1 || echo "Ensembl下载失败，尝试备用方案"
    fi

    # 选项2: 从UCSC下载 (备用)
    if [ ! -f "chr22.fa.gz" ]; then
        echo "尝试从UCSC下载..."
        wget -q --show-progress -O chr22.fa.gz \
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz" \
            2>&1 || echo "UCSC下载也失败"
    fi

    # 解压并提取片段
    if [ -f "chr22.fa.gz" ]; then
        echo "解压参考基因组..."
        gunzip chr22.fa.gz

        echo "提取5Mb片段用于实验..."
        if command -v samtools &> /dev/null; then
            samtools faidx chr22.fa
            # 提取22号染色体的10,000,000-15,000,000区域 (5Mb)
            samtools faidx chr22.fa 22:10000000-15000000 > chr22_fragment.fa || \
            samtools faidx chr22.fa chr22:10000000-15000000 > chr22_fragment.fa

            rm chr22.fa chr22.fa.fai
            echo "✓ 参考基因组准备完成"
        else
            echo "⚠ samtools未安装，使用完整22号染色体"
            mv chr22.fa chr22_fragment.fa
        fi
    else
        echo "下载失败，生成模拟参考序列..."
        generate_mock_reference
    fi
fi

# 创建参考基因组索引
if [ -f "chr22_fragment.fa" ]; then
    echo "创建参考基因组索引..."
    samtools faidx chr22_fragment.fa 2>/dev/null || true
    echo "✓ 参考基因组索引完成"
fi

# ============================================
# 2. 准备或生成BAM文件
# ============================================
echo ""
echo "步骤 2: 准备测序数据和比对文件..."

generate_mock_bam() {
    echo "生成模拟测序数据和比对文件..."

    python3 << 'EOF'
import random
import gzip
import os
import subprocess

print("生成模拟FASTQ文件...")

def generate_read(ref_seq, pos, length=150, error_rate=0.001):
    """从参考序列生成一条read"""
    if pos + length > len(ref_seq):
        return None

    read_seq = ref_seq[pos:pos+length]
    # 添加随机错误
    bases = ['A', 'T', 'G', 'C']
    read_list = list(read_seq)
    for i in range(len(read_list)):
        if random.random() < error_rate:
            read_list[i] = random.choice([b for b in bases if b != read_list[i]])

    # 生成质量分数
    qualities = ''.join([chr(random.randint(30, 40) + 33) for _ in range(length)])

    return ''.join(read_list), qualities

# 读取参考序列
print("读取参考基因组...")
ref_seq = ""
with open('chr22_fragment.fa', 'r') as f:
    for line in f:
        if not line.startswith('>'):
            ref_seq += line.strip().upper()

print(f"参考序列长度: {len(ref_seq)} bp")

# 生成配对端测序数据
num_reads = 100000  # 10万条reads
read_length = 150
insert_size = 350

print(f"生成 {num_reads} 条配对reads...")

with gzip.open('sample_R1.fastq.gz', 'wt') as f1, \
     gzip.open('sample_R2.fastq.gz', 'wt') as f2:

    for i in range(num_reads):
        # 随机选择起始位置
        pos = random.randint(0, len(ref_seq) - insert_size - read_length)

        # 生成R1 (正向)
        read1_data = generate_read(ref_seq, pos, read_length)
        if read1_data:
            seq1, qual1 = read1_data
            f1.write(f"@read_{i+1}/1\n{seq1}\n+\n{qual1}\n")

        # 生成R2 (反向互补)
        r2_pos = pos + insert_size - read_length
        read2_data = generate_read(ref_seq, r2_pos, read_length)
        if read2_data:
            seq2, qual2 = read2_data
            # 反向互补
            complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
            seq2_rc = ''.join([complement.get(b, 'N') for b in seq2[::-1]])
            f2.write(f"@read_{i+1}/2\n{seq2_rc}\n+\n{qual2[::-1]}\n")

        if (i + 1) % 10000 == 0:
            print(f"  已生成 {i+1} / {num_reads} reads")

print("✓ FASTQ文件生成完成")

# 进行比对
print("\n生成比对文件...")
if os.system("command -v bwa &> /dev/null") == 0:
    print("使用BWA进行比对...")

    # 构建索引
    if not os.path.exists('chr22_fragment.fa.bwt'):
        print("构建BWA索引...")
        os.system("bwa index chr22_fragment.fa 2>&1 | grep -v 'Ignore'")

    # 比对
    print("执行比对...")
    os.system("""
        bwa mem -t 2 chr22_fragment.fa sample_R1.fastq.gz sample_R2.fastq.gz 2>/dev/null | \
        samtools view -bS - | \
        samtools sort -o sample.bam - 2>/dev/null
    """)

    # 创建索引
    os.system("samtools index sample.bam 2>/dev/null")

    print("✓ BAM文件生成完成")
else:
    print("⚠ BWA未安装，无法生成BAM文件")
    print("  请安装BWA: conda install -c bioconda bwa")

EOF

    if [ -f "sample.bam" ]; then
        echo "✓ 模拟BAM文件生成成功"
        samtools flagstat sample.bam
    fi
}

# 检查是否可以从lesson-03复制数据
LESSON03_BAM="$HOME/ngs-analysis/lesson-03/results/bwa_default_sorted.bam"
if [ -f "$LESSON03_BAM" ] && [ ! -f "sample.bam" ]; then
    echo "发现lesson-03的比对结果，复制过来使用..."
    cp "$LESSON03_BAM" sample.bam
    cp "${LESSON03_BAM}.bai" sample.bam.bai 2>/dev/null || samtools index sample.bam
    echo "✓ 使用lesson-03的比对数据"
elif [ ! -f "sample.bam" ]; then
    generate_mock_bam
fi

# ============================================
# 3. 准备dbSNP和HapMap资源文件
# ============================================
echo ""
echo "步骤 3: 准备变异资源文件..."

generate_mock_resources() {
    echo "生成模拟的dbSNP和HapMap文件..."

    python3 << 'EOF'
import random

# 读取参考序列获取染色体信息
chr_name = ""
chr_length = 0
with open('chr22_fragment.fa', 'r') as f:
    header = f.readline().strip()
    chr_name = header.split()[0][1:]  # 去掉>符号

    seq_length = 0
    for line in f:
        if not line.startswith('>'):
            seq_length += len(line.strip())
    chr_length = seq_length

print(f"染色体: {chr_name}, 长度: {chr_length}")

# 生成dbSNP文件 (已知变异位点)
print("生成dbSNP文件...")
with open('dbsnp.vcf', 'w') as f:
    # VCF头部
    f.write("##fileformat=VCFv4.2\n")
    f.write("##source=MockDBSNP\n")
    f.write(f"##contig=<ID={chr_name},length={chr_length}>\n")
    f.write("##INFO=<ID=RS,Number=1,Type=String,Description=\"dbSNP ID\">\n")
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    # 生成1000个已知SNP
    positions = sorted(random.sample(range(1, chr_length), min(1000, chr_length-1)))
    bases = ['A', 'T', 'G', 'C']

    for i, pos in enumerate(positions):
        ref = random.choice(bases)
        alt = random.choice([b for b in bases if b != ref])
        rs_id = f"rs{1000000 + i}"
        f.write(f"{chr_name}\t{pos}\t{rs_id}\t{ref}\t{alt}\t.\tPASS\tRS={rs_id}\n")

print("✓ dbSNP文件生成完成")

# 生成HapMap文件 (高质量SNP)
print("生成HapMap文件...")
with open('hapmap.vcf', 'w') as f:
    # VCF头部
    f.write("##fileformat=VCFv4.2\n")
    f.write("##source=MockHapMap\n")
    f.write(f"##contig=<ID={chr_name},length={chr_length}>\n")
    f.write("##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n")
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    # 生成500个高质量SNP
    positions = sorted(random.sample(range(1, chr_length), min(500, chr_length-1)))

    for i, pos in enumerate(positions):
        ref = random.choice(bases)
        alt = random.choice([b for b in bases if b != ref])
        snp_id = f"hapmap_{i+1}"
        af = random.uniform(0.1, 0.5)  # 等位基因频率
        f.write(f"{chr_name}\t{pos}\t{snp_id}\t{ref}\t{alt}\t100\tPASS\tAF={af:.3f}\n")

print("✓ HapMap文件生成完成")

EOF
}

if [ ! -f "dbsnp.vcf" ] || [ ! -f "hapmap.vcf" ]; then
    generate_mock_resources
fi

# 创建VCF索引
if command -v bgzip &> /dev/null && command -v tabix &> /dev/null; then
    echo "创建VCF索引..."
    for vcf in dbsnp.vcf hapmap.vcf; do
        if [ ! -f "${vcf}.gz" ]; then
            bgzip -c $vcf > ${vcf}.gz
            tabix -p vcf ${vcf}.gz
        fi
    done
    echo "✓ VCF索引创建完成"
else
    echo "⚠ bgzip/tabix未安装，跳过VCF索引创建"
fi

# ============================================
# 4. 验证数据
# ============================================
echo ""
echo "步骤 4: 验证数据完整性..."

echo "检查文件列表:"
ls -lh

echo ""
echo "参考基因组信息:"
if [ -f "chr22_fragment.fa" ]; then
    echo "  文件大小: $(ls -lh chr22_fragment.fa | awk '{print $5}')"
    echo "  序列数量: $(grep -c '^>' chr22_fragment.fa)"
    head -2 chr22_fragment.fa
fi

echo ""
echo "BAM文件信息:"
if [ -f "sample.bam" ]; then
    echo "  文件大小: $(ls -lh sample.bam | awk '{print $5}')"
    samtools flagstat sample.bam 2>/dev/null | head -5
fi

echo ""
echo "资源文件信息:"
if [ -f "dbsnp.vcf" ]; then
    echo "  dbSNP变异数: $(grep -v '^#' dbsnp.vcf | wc -l)"
fi
if [ -f "hapmap.vcf" ]; then
    echo "  HapMap变异数: $(grep -v '^#' hapmap.vcf | wc -l)"
fi

# 创建数据说明文件
cat > DATA_INFO.txt << 'EOFINFO'
变异检测实验数据说明
====================

生成时间: $(date)

文件列表:
---------
1. chr22_fragment.fa          - 参考基因组 (人类22号染色体5Mb片段)
2. chr22_fragment.fa.fai      - 参考基因组索引
3. sample.bam                 - 比对后的测序数据
4. sample.bam.bai            - BAM文件索引
5. dbsnp.vcf                 - 已知变异数据库
6. hapmap.vcf                - 高质量SNP集合

数据特征:
---------
- 参考基因组: 人类GRCh38 22号染色体片段 (~5Mb)
- 测序数据: 模拟Illumina双端测序 (2×150bp)
- 覆盖深度: ~30X
- Reads数量: ~100,000 pairs
- 已知变异: ~1,000个dbSNP位点, ~500个HapMap位点

数据来源:
---------
- 参考基因组: Ensembl/UCSC公共数据库
- 测序数据: 模拟生成或来自lesson-03实验
- 变异资源: 模拟生成用于教学演示

使用说明:
---------
这些数据专门为变异检测练习设计，包含了GATK流程所需的所有输入文件。
数据规模适中，适合在课堂实验环境中使用。

注意事项:
---------
- 这是教学用模拟数据，不应用于实际研究
- 变异资源文件是简化版本，真实分析需使用官方数据库
- 建议定期检查数据完整性

EOFINFO

echo ""
echo "=========================================="
echo "数据准备完成！"
echo "=========================================="
echo ""
echo "数据位置: $(pwd)"
echo "数据说明: DATA_INFO.txt"
echo ""
echo "下一步: 运行 setup.sh 配置工作环境"
echo "=========================================="
