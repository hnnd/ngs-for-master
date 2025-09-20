#!/bin/bash
# 第3次课环境设置脚本
# 作者：王运生
# 日期：2025-01-20
# 用法：bash setup.sh

set -e  # 遇到错误立即退出
set -u  # 使用未定义变量时报错

echo "=========================================="
echo "高通量测序序列比对实验环境设置"
echo "=========================================="

# 检查当前目录
if [[ ! -d "scripts" ]]; then
    echo "错误：请在实验根目录运行此脚本"
    exit 1
fi

# 创建必要的目录
echo "创建实验目录结构..."
mkdir -p data reference results logs scripts

# 检查必需软件
echo "检查软件环境..."

check_software() {
    local software=$1
    local command=$2
    
    if command -v $command &> /dev/null; then
        echo "✓ $software 已安装"
        $command --version 2>&1 | head -1 || $command 2>&1 | head -1
    else
        echo "✗ $software 未安装"
        echo "请使用以下命令安装："
        case $software in
            "BWA")
                echo "  conda install -c bioconda bwa"
                ;;
            "Bowtie2")
                echo "  conda install -c bioconda bowtie2"
                ;;
            "samtools")
                echo "  conda install -c bioconda samtools"
                ;;
            "Python")
                echo "  conda install python>=3.8"
                ;;
        esac
        return 1
    fi
}

# 检查所有必需软件
software_ok=true
check_software "BWA" "bwa" || software_ok=false
check_software "Bowtie2" "bowtie2" || software_ok=false
check_software "samtools" "samtools" || software_ok=false
check_software "Python" "python" || software_ok=false

if [ "$software_ok" = false ]; then
    echo ""
    echo "请先安装缺失的软件，然后重新运行此脚本"
    exit 1
fi

# 检查Python包
echo ""
echo "检查Python包..."
python -c "
import sys
required_packages = ['pandas', 'matplotlib', 'seaborn', 'numpy']
missing_packages = []

for package in required_packages:
    try:
        __import__(package)
        print(f'✓ {package} 已安装')
    except ImportError:
        print(f'✗ {package} 未安装')
        missing_packages.append(package)

if missing_packages:
    print(f'请安装缺失的包: pip install {\" \".join(missing_packages)}')
    sys.exit(1)
"

if [ $? -ne 0 ]; then
    echo "Python包检查失败，请安装缺失的包"
    exit 1
fi

# 设置环境变量
echo ""
echo "设置环境变量..."
export OMP_NUM_THREADS=4
export TMPDIR=$(pwd)/tmp
mkdir -p $TMPDIR

# 创建配置文件
echo "创建配置文件..."
cat > config.txt << EOF
# 实验配置文件
REFERENCE_FILE=reference/chr22_fragment.fa
SAMPLE_R1=data/sample_R1.fastq
SAMPLE_R2=data/sample_R2.fastq
THREADS=4
MEMORY=8G
OUTPUT_DIR=results
LOG_DIR=logs
EOF

# 创建示例参考基因组（如果不存在）
if [ ! -f "reference/chr22_fragment.fa" ]; then
    echo ""
    echo "创建示例参考基因组..."
    mkdir -p reference
    
    # 生成一个5Mb的模拟基因组片段
    python << 'EOF'
import random
import os

# 设置随机种子以确保可重现性
random.seed(42)

# 生成5Mb的随机DNA序列
sequence_length = 5000000  # 5Mb
bases = ['A', 'T', 'G', 'C']

print("生成5Mb模拟基因组序列...")
sequence = ''.join(random.choices(bases, k=sequence_length))

# 写入FASTA文件
os.makedirs('reference', exist_ok=True)
with open('reference/chr22_fragment.fa', 'w') as f:
    f.write('>chr22_fragment simulated 5Mb sequence\n')
    
    # 每行80个字符
    for i in range(0, len(sequence), 80):
        f.write(sequence[i:i+80] + '\n')

print(f"参考基因组已生成: reference/chr22_fragment.fa ({len(sequence)} bp)")
EOF
fi

# 创建示例测序数据（如果不存在）
if [ ! -f "data/sample_R1.fastq" ]; then
    echo ""
    echo "创建示例测序数据..."
    mkdir -p data
    
    python << 'EOF'
import random
import os

# 设置随机种子
random.seed(123)

# 读取参考基因组
print("读取参考基因组...")
with open('reference/chr22_fragment.fa', 'r') as f:
    lines = f.readlines()
    reference = ''.join(line.strip() for line in lines[1:])  # 跳过header

print(f"参考基因组长度: {len(reference)} bp")

# 生成测序reads参数
num_reads = 500000  # 50万对reads
read_length = 100
insert_size_mean = 300
insert_size_std = 50
error_rate = 0.01

bases = ['A', 'T', 'G', 'C']
complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

def reverse_complement(seq):
    return ''.join(complement[base] for base in reversed(seq))

def introduce_errors(seq, error_rate):
    """在序列中引入随机错误"""
    seq_list = list(seq)
    for i in range(len(seq_list)):
        if random.random() < error_rate:
            seq_list[i] = random.choice(bases)
    return ''.join(seq_list)

def generate_quality_scores(length):
    """生成质量分数（Phred+33格式）"""
    # 生成质量分数，大部分为高质量
    scores = []
    for i in range(length):
        if i < 10 or i > length - 10:  # 两端质量稍低
            q = random.randint(20, 35)
        else:  # 中间质量较高
            q = random.randint(30, 40)
        scores.append(chr(q + 33))
    return ''.join(scores)

print("生成paired-end测序数据...")

os.makedirs('data', exist_ok=True)

with open('data/sample_R1.fastq', 'w') as f1, open('data/sample_R2.fastq', 'w') as f2:
    for read_id in range(num_reads):
        # 随机选择起始位置
        max_start = len(reference) - insert_size_mean - 2 * read_length
        if max_start <= 0:
            continue
            
        start_pos = random.randint(0, max_start)
        
        # 生成插入片段大小
        insert_size = max(read_length * 2, 
                         int(random.gauss(insert_size_mean, insert_size_std)))
        
        # 提取R1和R2序列
        r1_seq = reference[start_pos:start_pos + read_length]
        r2_start = start_pos + insert_size - read_length
        r2_seq = reference[r2_start:r2_start + read_length]
        r2_seq = reverse_complement(r2_seq)
        
        # 引入测序错误
        r1_seq = introduce_errors(r1_seq, error_rate)
        r2_seq = introduce_errors(r2_seq, error_rate)
        
        # 生成质量分数
        r1_qual = generate_quality_scores(len(r1_seq))
        r2_qual = generate_quality_scores(len(r2_seq))
        
        # 写入FASTQ格式
        read_name = f"read_{read_id:06d}"
        
        # R1
        f1.write(f"@{read_name}/1\n")
        f1.write(f"{r1_seq}\n")
        f1.write(f"+\n")
        f1.write(f"{r1_qual}\n")
        
        # R2
        f2.write(f"@{read_name}/2\n")
        f2.write(f"{r2_seq}\n")
        f2.write(f"+\n")
        f2.write(f"{r2_qual}\n")
        
        if (read_id + 1) % 50000 == 0:
            print(f"已生成 {read_id + 1} 对reads...")

print(f"测序数据生成完成:")
print(f"  R1: data/sample_R1.fastq")
print(f"  R2: data/sample_R2.fastq")
print(f"  总reads数: {num_reads} pairs")
EOF
fi

# 检查生成的文件
echo ""
echo "检查生成的文件..."
if [ -f "reference/chr22_fragment.fa" ]; then
    ref_size=$(wc -c < reference/chr22_fragment.fa)
    echo "✓ 参考基因组: reference/chr22_fragment.fa (${ref_size} bytes)"
fi

if [ -f "data/sample_R1.fastq" ] && [ -f "data/sample_R2.fastq" ]; then
    r1_lines=$(wc -l < data/sample_R1.fastq)
    r2_lines=$(wc -l < data/sample_R2.fastq)
    reads_count=$((r1_lines / 4))
    echo "✓ 测序数据: data/sample_R*.fastq (${reads_count} reads each)"
fi

# 创建日志文件
echo ""
echo "初始化日志文件..."
cat > logs/experiment.log << EOF
实验开始时间: $(date)
参考基因组: reference/chr22_fragment.fa
测序数据: data/sample_R1.fastq, data/sample_R2.fastq
线程数: $OMP_NUM_THREADS
临时目录: $TMPDIR
EOF

echo ""
echo "=========================================="
echo "环境设置完成！"
echo "=========================================="
echo "配置信息："
echo "- 工作目录: $(pwd)"
echo "- 线程数: $OMP_NUM_THREADS"
echo "- 临时目录: $TMPDIR"
echo "- 配置文件: config.txt"
echo "- 日志文件: logs/experiment.log"
echo ""
echo "现在可以开始进行序列比对实验了！"
echo "请按照实验手册的步骤继续操作。"