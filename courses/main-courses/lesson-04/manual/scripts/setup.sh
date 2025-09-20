#!/bin/bash
# 第4次课环境设置脚本
# 作者：王运生
# 日期：2025-01-20
# 用法：bash setup.sh

set -e  # 遇到错误立即退出
set -u  # 使用未定义变量时报错

echo "=== 第4次课：变异检测与基因分型 - 环境设置 ==="

# 检查并创建工作目录
WORK_DIR="$HOME/ngs-analysis/lesson-04"
echo "创建工作目录: $WORK_DIR"
mkdir -p "$WORK_DIR"/{data,scripts,results,logs,plots}

cd "$WORK_DIR"

# 检查必要软件
echo "检查软件环境..."

check_software() {
    local software=$1
    local version_cmd=$2
    
    if command -v "$software" &> /dev/null; then
        echo "✓ $software 已安装"
        eval "$version_cmd" 2>/dev/null || echo "  版本信息获取失败"
    else
        echo "✗ $software 未安装"
        return 1
    fi
}

# 检查各个软件
check_software "gatk" "gatk --version | head -1"
check_software "samtools" "samtools --version | head -1"
check_software "bcftools" "bcftools --version | head -1"
check_software "vep" "vep --help | head -1"
check_software "python3" "python3 --version"

# 检查Java环境（GATK需要）
if command -v java &> /dev/null; then
    echo "✓ Java 已安装"
    java -version 2>&1 | head -1
else
    echo "✗ Java 未安装，GATK需要Java环境"
fi

# 设置环境变量
echo "设置环境变量..."
export GATK_JAVA_OPTIONS="-Xmx8g -XX:+UseParallelGC"
echo "GATK_JAVA_OPTIONS=$GATK_JAVA_OPTIONS"

# 检查数据文件
echo "检查数据文件..."
DATA_SOURCES=(
    "/data/ngs/lesson04/sample.bam"
    "/data/reference/hg38/reference.fa"
    "/data/reference/dbsnp/dbsnp.vcf"
    "/data/reference/hapmap/hapmap.vcf"
)

for file in "${DATA_SOURCES[@]}"; do
    if [[ -f "$file" ]]; then
        echo "✓ 找到数据文件: $file"
    else
        echo "✗ 缺少数据文件: $file"
        echo "  请联系管理员准备数据文件"
    fi
done

# 创建符号链接（如果数据文件存在）
echo "创建数据文件链接..."
if [[ -f "/data/ngs/lesson04/sample.bam" ]]; then
    ln -sf /data/ngs/lesson04/sample.bam* data/
    echo "✓ 链接测序数据文件"
fi

if [[ -f "/data/reference/hg38/reference.fa" ]]; then
    ln -sf /data/reference/hg38/reference.fa* data/
    echo "✓ 链接参考基因组文件"
fi

if [[ -f "/data/reference/dbsnp/dbsnp.vcf" ]]; then
    ln -sf /data/reference/dbsnp/dbsnp.vcf* data/
    echo "✓ 链接dbSNP数据库"
fi

if [[ -f "/data/reference/hapmap/hapmap.vcf" ]]; then
    ln -sf /data/reference/hapmap/hapmap.vcf* data/
    echo "✓ 链接HapMap数据库"
fi

# 创建日志文件
touch logs/analysis.log
echo "✓ 创建日志文件"

# 检查磁盘空间
echo "检查磁盘空间..."
AVAILABLE_SPACE=$(df -h . | awk 'NR==2 {print $4}')
echo "可用空间: $AVAILABLE_SPACE"

# 创建配置文件
cat > config.txt << EOF
# 第4次课配置文件
WORK_DIR=$WORK_DIR
REFERENCE=data/reference.fa
DBSNP=data/dbsnp.vcf
HAPMAP=data/hapmap.vcf
INPUT_BAM=data/sample.bam
THREADS=4
MEMORY=8g
EOF

echo "✓ 创建配置文件: config.txt"

# 显示目录结构
echo "工作目录结构:"
tree . 2>/dev/null || ls -la

echo ""
echo "=== 环境设置完成 ==="
echo "工作目录: $WORK_DIR"
echo "请确保所有软件都已正确安装"
echo "如有问题，请联系: wangys@hunau.edu.cn"