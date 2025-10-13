#!/bin/bash

################################################################################
# 实验数据准备脚本
# 课程：高通量测序数据分析 - 第3次课
# 功能：自动下载参考基因组和测序数据
################################################################################

set -e  # 遇到错误立即退出

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# 打印信息函数
print_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# 检查必需的工具
check_dependencies() {
    print_info "检查必需的软件..."

    local missing_tools=()

    for tool in wget samtools; do
        if ! command -v $tool &> /dev/null; then
            missing_tools+=($tool)
        fi
    done

    if [ ${#missing_tools[@]} -ne 0 ]; then
        print_error "缺少以下工具: ${missing_tools[*]}"
        print_info "请先安装: conda install -c bioconda ${missing_tools[*]}"
        exit 1
    fi

    print_info "✓ 所有必需工具已安装"
}

# 下载参考基因组
download_reference() {
    print_info "========================================="
    print_info "步骤1: 下载参考基因组"
    print_info "========================================="

    if [ -f "reference/chr22_fragment.fa" ]; then
        print_warning "参考基因组已存在，跳过下载"
        return
    fi

    mkdir -p reference

    # 尝试从Ensembl下载
    print_info "从Ensembl下载人类22号染色体..."
    if wget -q --show-progress -O reference/chr22.fa.gz \
        "http://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz"; then

        print_info "下载成功，正在解压..."
        gunzip reference/chr22.fa.gz

        print_info "提取前5Mb序列..."
        samtools faidx reference/chr22.fa
        samtools faidx reference/chr22.fa 22:1-5000000 > reference/chr22_fragment.fa

        # 清理临时文件
        rm reference/chr22.fa reference/chr22.fa.fai

        print_info "✓ 参考基因组准备完成"
    else
        print_error "从Ensembl下载失败"

        # 尝试备用源
        print_info "尝试从UCSC下载..."
        if wget -q --show-progress -O reference/chr22.fa.gz \
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz"; then

            gunzip reference/chr22.fa.gz
            samtools faidx reference/chr22.fa
            samtools faidx reference/chr22.fa chr22:1-5000000 > reference/chr22_fragment.fa
            rm reference/chr22.fa reference/chr22.fa.fai

            print_info "✓ 参考基因组准备完成（使用UCSC源）"
        else
            print_error "所有数据源均无法访问，请检查网络连接"
            exit 1
        fi
    fi
}

# 下载或生成测序数据
prepare_sequencing_data() {
    print_info "========================================="
    print_info "步骤2: 准备测序数据"
    print_info "========================================="

    if [ -f "data/sample_R1.fastq" ] && [ -f "data/sample_R2.fastq" ]; then
        print_warning "测序数据已存在，跳过准备"
        return
    fi

    mkdir -p data

    # 检查是否有SRA toolkit
    if command -v fastq-dump &> /dev/null; then
        print_info "检测到SRA toolkit，尝试下载真实数据..."

        if fastq-dump --split-files -X 1000000 SRR622461 -O data/ 2>/dev/null; then
            mv data/SRR622461_1.fastq data/sample_R1.fastq
            mv data/SRR622461_2.fastq data/sample_R2.fastq
            print_info "✓ 真实测序数据下载完成"
            return
        else
            print_warning "SRA数据下载失败，将生成模拟数据"
        fi
    fi

    # 检查是否有wgsim
    if command -v wgsim &> /dev/null; then
        print_info "使用wgsim生成模拟测序数据..."

        if [ ! -f "reference/chr22_fragment.fa" ]; then
            print_error "参考基因组不存在，请先运行参考基因组下载步骤"
            exit 1
        fi

        wgsim -N 1000000 -1 100 -2 100 -r 0.001 -R 0.15 -X 0.3 \
            reference/chr22_fragment.fa \
            data/sample_R1.fastq \
            data/sample_R2.fastq

        print_info "✓ 模拟测序数据生成完成"
    else
        print_error "无法获取测序数据"
        print_info "请安装以下工具之一："
        print_info "  - SRA toolkit: conda install -c bioconda sra-tools"
        print_info "  - wgsim: conda install -c bioconda wgsim"
        exit 1
    fi
}

# 验证数据
validate_data() {
    print_info "========================================="
    print_info "步骤3: 验证数据完整性"
    print_info "========================================="

    # 检查参考基因组
    if [ -f "reference/chr22_fragment.fa" ]; then
        ref_size=$(wc -c < reference/chr22_fragment.fa)
        print_info "✓ 参考基因组: chr22_fragment.fa ($ref_size bytes)"

        # 显示FASTA头部
        header=$(head -1 reference/chr22_fragment.fa)
        print_info "  序列名称: $header"
    else
        print_error "✗ 参考基因组文件缺失"
        exit 1
    fi

    # 检查测序数据
    if [ -f "data/sample_R1.fastq" ] && [ -f "data/sample_R2.fastq" ]; then
        r1_reads=$(($(wc -l < data/sample_R1.fastq) / 4))
        r2_reads=$(($(wc -l < data/sample_R2.fastq) / 4))

        print_info "✓ 测序数据准备完成"
        print_info "  R1 reads数量: $r1_reads"
        print_info "  R2 reads数量: $r2_reads"

        if [ $r1_reads -ne $r2_reads ]; then
            print_warning "R1和R2的reads数量不匹配"
        fi
    else
        print_error "✗ 测序数据文件缺失"
        exit 1
    fi

    print_info "========================================="
    print_info "✓ 所有数据准备完成，可以开始实验"
    print_info "========================================="
}

# 显示帮助信息
show_help() {
    cat << EOF
用法: $0 [选项]

选项:
    -h, --help          显示此帮助信息
    -r, --reference     仅下载参考基因组
    -s, --sequencing    仅准备测序数据
    -v, --validate      仅验证数据

示例:
    $0                  # 下载所有数据
    $0 -r              # 仅下载参考基因组
    $0 -s              # 仅准备测序数据

EOF
}

# 主函数
main() {
    echo ""
    print_info "高通量测序数据分析 - 实验数据准备脚本"
    echo ""

    # 解析命令行参数
    case "${1:-all}" in
        -h|--help)
            show_help
            exit 0
            ;;
        -r|--reference)
            check_dependencies
            download_reference
            ;;
        -s|--sequencing)
            check_dependencies
            prepare_sequencing_data
            ;;
        -v|--validate)
            validate_data
            ;;
        all|*)
            check_dependencies
            download_reference
            prepare_sequencing_data
            validate_data
            ;;
    esac

    echo ""
}

# 运行主函数
main "$@"
