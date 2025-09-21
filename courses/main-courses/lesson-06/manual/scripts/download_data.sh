#!/bin/bash
# ChIP-seq数据下载脚本
# 作者：王运生
# 日期：2025-01-21
# 用法：bash download_data.sh

set -e  # 遇到错误立即退出
set -u  # 使用未定义变量时报错

# 加载配置
source ~/ngs-analysis/lesson-06/scripts/config.sh

echo "=== ChIP-seq数据下载 ==="
echo "目标目录：$DATA_DIR"

# 创建数据目录
mkdir -p $DATA_DIR
cd $DATA_DIR

# 数据服务器URL（示例）
SERVER_URL="http://example.server.com/chipseq_data"

# 定义下载文件列表
declare -A files=(
    ["H3K4me3_ChIP.fastq.gz"]="2.1GB"
    ["Input_control.fastq.gz"]="2.0GB"
    ["hg38.fa"]="3.0GB"
    ["hg38.fa.fai"]="3.2KB"
    ["gencode.v38.gtf"]="51MB"
)

# MD5校验文件
declare -A md5sums=(
    ["H3K4me3_ChIP.fastq.gz"]="a1b2c3d4e5f6789012345678901234567"
    ["Input_control.fastq.gz"]="b2c3d4e5f6789012345678901234567a1"
    ["hg38.fa"]="c3d4e5f6789012345678901234567a1b2"
    ["hg38.fa.fai"]="d4e5f6789012345678901234567a1b2c3"
    ["gencode.v38.gtf"]="e5f6789012345678901234567a1b2c3d4"
)

# 检查网络连接
echo "检查网络连接..."
if ! ping -c 1 google.com &> /dev/null; then
    echo "⚠️  网络连接异常，请检查网络设置"
    echo "如果在内网环境，请联系管理员获取数据"
    exit 1
fi

# 下载函数
download_file() {
    local filename=$1
    local expected_size=$2
    local url="$SERVER_URL/$filename"
    
    echo ""
    echo "下载文件：$filename (预期大小：$expected_size)"
    
    # 检查文件是否已存在
    if [ -f "$filename" ]; then
        echo "文件已存在，检查完整性..."
        if verify_file "$filename"; then
            echo "✓ 文件完整，跳过下载"
            return 0
        else
            echo "文件损坏，重新下载..."
            rm -f "$filename"
        fi
    fi
    
    # 使用wget下载，支持断点续传
    if command -v wget &> /dev/null; then
        wget -c -t 3 -T 30 "$url" -O "$filename"
    elif command -v curl &> /dev/null; then
        curl -C - -L -o "$filename" "$url"
    else
        echo "错误：未找到wget或curl下载工具"
        exit 1
    fi
    
    # 验证下载的文件
    if verify_file "$filename"; then
        echo "✓ $filename 下载完成并验证成功"
    else
        echo "✗ $filename 下载失败或文件损坏"
        return 1
    fi
}

# 文件验证函数
verify_file() {
    local filename=$1
    
    # 检查文件是否存在
    if [ ! -f "$filename" ]; then
        return 1
    fi
    
    # 检查文件大小（简单检查）
    local file_size=$(stat -c%s "$filename" 2>/dev/null || stat -f%z "$filename" 2>/dev/null)
    if [ "$file_size" -lt 1000 ]; then
        echo "文件大小异常：$file_size bytes"
        return 1
    fi
    
    # MD5校验（如果有校验值）
    if [ -n "${md5sums[$filename]:-}" ]; then
        local expected_md5="${md5sums[$filename]}"
        local actual_md5
        
        if command -v md5sum &> /dev/null; then
            actual_md5=$(md5sum "$filename" | cut -d' ' -f1)
        elif command -v md5 &> /dev/null; then
            actual_md5=$(md5 -q "$filename")
        else
            echo "警告：无法进行MD5校验"
            return 0
        fi
        
        if [ "$actual_md5" != "$expected_md5" ]; then
            echo "MD5校验失败："
            echo "  期望：$expected_md5"
            echo "  实际：$actual_md5"
            return 1
        fi
    fi
    
    return 0
}

# 主下载流程
echo "开始下载ChIP-seq分析所需数据..."
echo "注意：总下载量约7GB，请确保网络稳定"

# 检查磁盘空间
available_space=$(df -BG "$DATA_DIR" | awk 'NR==2 {print $4}' | sed 's/G//')
if [ "$available_space" -lt 10 ]; then
    echo "错误：磁盘空间不足，需要至少10GB空间"
    exit 1
fi

# 下载所有文件
failed_files=()
for filename in "${!files[@]}"; do
    if ! download_file "$filename" "${files[$filename]}"; then
        failed_files+=("$filename")
    fi
done

# 检查下载结果
if [ ${#failed_files[@]} -eq 0 ]; then
    echo ""
    echo "=== 所有文件下载完成 ==="
    echo "下载位置：$DATA_DIR"
    echo ""
    echo "文件列表："
    ls -lh $DATA_DIR
    
    # 创建数据清单
    echo ""
    echo "创建数据清单..."
    cat > $DATA_DIR/data_manifest.txt << EOF
ChIP-seq数据清单
下载时间：$(date)
下载位置：$DATA_DIR

文件列表：
$(ls -lh $DATA_DIR | grep -v "^total")

文件说明：
- H3K4me3_ChIP.fastq.gz: H3K4me3 ChIP-seq测序数据
- Input_control.fastq.gz: Input对照测序数据
- hg38.fa: 人类参考基因组序列（GRCh38/hg38）
- hg38.fa.fai: 参考基因组索引文件
- gencode.v38.gtf: GENCODE基因注释文件

数据来源：ENCODE项目公开数据
用途：ChIP-seq数据分析教学实验
EOF
    
    echo "✓ 数据清单已创建：$DATA_DIR/data_manifest.txt"
    
else
    echo ""
    echo "=== 部分文件下载失败 ==="
    echo "失败文件："
    for file in "${failed_files[@]}"; do
        echo "  - $file"
    done
    echo ""
    echo "请检查网络连接后重新运行此脚本"
    exit 1
fi

# 准备参考基因组索引
echo ""
echo "准备参考基因组索引..."
if [ ! -f "hg38.fa.bwt" ]; then
    echo "构建BWA索引（这可能需要几分钟）..."
    bwa index hg38.fa
    echo "✓ BWA索引构建完成"
else
    echo "✓ BWA索引已存在"
fi

# 验证数据完整性
echo ""
echo "最终验证..."
all_good=true

# 检查FASTQ文件
for fastq in H3K4me3_ChIP.fastq.gz Input_control.fastq.gz; do
    if ! zcat "$fastq" | head -4 | grep -q "^@"; then
        echo "✗ $fastq 格式异常"
        all_good=false
    else
        echo "✓ $fastq 格式正常"
    fi
done

# 检查FASTA文件
if ! head -1 hg38.fa | grep -q "^>"; then
    echo "✗ hg38.fa 格式异常"
    all_good=false
else
    echo "✓ hg38.fa 格式正常"
fi

# 检查GTF文件
if ! head -1 gencode.v38.gtf | grep -q -E "^(#|chr)"; then
    echo "✗ gencode.v38.gtf 格式异常"
    all_good=false
else
    echo "✓ gencode.v38.gtf 格式正常"
fi

if $all_good; then
    echo ""
    echo "=== 数据下载和验证完成 ==="
    echo "所有文件已准备就绪，可以开始ChIP-seq分析"
    echo ""
    echo "下一步："
    echo "1. 运行质量控制：bash quality_control.sh"
    echo "2. 或直接开始完整分析：bash chipseq_pipeline.sh"
else
    echo ""
    echo "=== 数据验证失败 ==="
    echo "请检查失败的文件并重新下载"
    exit 1
fi