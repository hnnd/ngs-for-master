#!/bin/bash
# ChIP-seq分析环境设置脚本
# 作者：王运生
# 日期：2025-01-21
# 用法：bash setup.sh

set -e  # 遇到错误立即退出
set -u  # 使用未定义变量时报错

echo "=== ChIP-seq分析环境设置 ==="

# 检查并创建目录结构
echo "创建目录结构..."
mkdir -p ~/ngs-analysis/lesson-06/{data,scripts,results,logs,figures}
mkdir -p ~/ngs-analysis/lesson-06/results/{fastqc,alignment,peaks,annotation}

# 检查必要软件
echo "检查软件环境..."
check_software() {
    local software=$1
    local command=$2
    
    if command -v $command &> /dev/null; then
        echo "✓ $software 已安装"
        $command --version 2>/dev/null || $command -v 2>/dev/null || echo "  版本信息不可用"
    else
        echo "✗ $software 未安装"
        return 1
    fi
}

# 检查软件列表
software_list=(
    "MACS2:macs2"
    "SAMtools:samtools"
    "BWA:bwa"
    "FastQC:fastqc"
    "deepTools:bamCoverage"
    "bedtools:bedtools"
    "R:R"
)

missing_software=()
for item in "${software_list[@]}"; do
    IFS=':' read -r name command <<< "$item"
    if ! check_software "$name" "$command"; then
        missing_software+=("$name")
    fi
done

# 如果有缺失软件，提供安装建议
if [ ${#missing_software[@]} -gt 0 ]; then
    echo ""
    echo "缺失软件安装建议："
    for software in "${missing_software[@]}"; do
        case $software in
            "MACS2")
                echo "  pip install MACS2"
                ;;
            "SAMtools"|"BWA"|"FastQC"|"bedtools")
                echo "  conda install $software"
                ;;
            "deepTools")
                echo "  pip install deepTools"
                ;;
            "R")
                echo "  conda install r-base"
                ;;
        esac
    done
    echo ""
    echo "请安装缺失软件后重新运行此脚本。"
    exit 1
fi

# 检查R包
echo ""
echo "检查R包..."
R --slave --no-restore --no-save << 'EOF'
required_packages <- c("ChIPseeker", "TxDb.Hsapiens.UCSC.hg38.knownGene", 
                       "clusterProfiler", "org.Hs.eg.db", "ggplot2")

missing_packages <- c()
for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        missing_packages <- c(missing_packages, pkg)
        cat("✗", pkg, "未安装\n")
    } else {
        cat("✓", pkg, "已安装\n")
    }
}

if (length(missing_packages) > 0) {
    cat("\n缺失R包安装命令：\n")
    cat("if (!requireNamespace('BiocManager', quietly = TRUE))\n")
    cat("    install.packages('BiocManager')\n")
    cat("BiocManager::install(c(")
    cat(paste0("'", missing_packages, "'", collapse = ", "))
    cat("))\n")
    quit(status = 1)
} else {
    cat("\n所有R包已安装完成。\n")
}
EOF

if [ $? -ne 0 ]; then
    echo "请在R中安装缺失的包后重新运行此脚本。"
    exit 1
fi

# 设置环境变量
echo ""
echo "设置环境变量..."
export LESSON_DIR=~/ngs-analysis/lesson-06
echo "export LESSON_DIR=~/ngs-analysis/lesson-06" >> ~/.bashrc

# 创建有用的别名
echo ""
echo "创建便捷别名..."
cat >> ~/.bashrc << 'EOF'

# ChIP-seq分析别名
alias chipseq-cd='cd ~/ngs-analysis/lesson-06'
alias chipseq-results='ls -la ~/ngs-analysis/lesson-06/results/'
alias chipseq-peaks='head ~/ngs-analysis/lesson-06/results/peaks/*.narrowPeak'
EOF

# 检查磁盘空间
echo ""
echo "检查磁盘空间..."
available_space=$(df -h ~ | awk 'NR==2 {print $4}' | sed 's/G//')
if [ "${available_space%.*}" -lt 20 ]; then
    echo "⚠️  警告：可用磁盘空间不足20GB，建议清理磁盘空间"
else
    echo "✓ 磁盘空间充足"
fi

# 检查内存
echo ""
echo "检查系统内存..."
total_mem=$(free -g | awk 'NR==2{print $2}')
if [ "$total_mem" -lt 8 ]; then
    echo "⚠️  警告：系统内存少于8GB，可能影响分析性能"
else
    echo "✓ 系统内存充足"
fi

# 创建示例配置文件
echo ""
echo "创建配置文件..."
cat > ~/ngs-analysis/lesson-06/scripts/config.sh << 'EOF'
#!/bin/bash
# ChIP-seq分析配置文件

# 基本路径
WORK_DIR=~/ngs-analysis/lesson-06
DATA_DIR=$WORK_DIR/data
RESULTS_DIR=$WORK_DIR/results
SCRIPTS_DIR=$WORK_DIR/scripts

# 参考基因组
GENOME_FA=$DATA_DIR/hg38.fa
GENOME_GTF=$DATA_DIR/gencode.v38.gtf

# 分析参数
THREADS=4
MACS2_QVALUE=0.05
MACS2_GSIZE="hs"

# 样本信息
CHIP_SAMPLE="H3K4me3_ChIP"
INPUT_SAMPLE="Input_control"

# 输出前缀
OUTPUT_PREFIX="H3K4me3"
EOF

echo ""
echo "=== 环境设置完成 ==="
echo "工作目录：~/ngs-analysis/lesson-06"
echo "配置文件：~/ngs-analysis/lesson-06/scripts/config.sh"
echo ""
echo "下一步："
echo "1. 运行 'source ~/.bashrc' 加载新的环境变量"
echo "2. 运行 'chipseq-cd' 进入工作目录"
echo "3. 运行数据下载脚本下载实验数据"
echo ""
echo "如有问题，请联系：wangys@hunau.edu.cn"