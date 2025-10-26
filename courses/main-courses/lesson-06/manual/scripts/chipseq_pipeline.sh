#!/bin/bash

###############################################################################
#          ChIP-seq完整分析流程脚本                                         #
# 集成所有分析步骤：质量控制、比对、Peak calling、注释和可视化           #
###############################################################################

set -e

# ============================================================================
# 配置参数
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
REFERENCE_DIR="$WORK_DIR/reference"
DATA_DIR="$WORK_DIR/data"
RESULTS_DIR="$WORK_DIR/../results"
QC_DIR="$WORK_DIR/../qc"
LOGS_DIR="$WORK_DIR/../logs"

# 创建必要的目录
mkdir -p "$RESULTS_DIR"/{fastqc,alignment,peaks,annotation,plots} "$QC_DIR" "$LOGS_DIR"

# 参数设置
THREADS=4
REFERENCE="$REFERENCE_DIR/hg38_chr22_demo.fa"
GTF="$REFERENCE_DIR/hg38_genes_demo.gtf"
CHIP_FASTQ_R1="$DATA_DIR/H3K4me3_ChIP_R1.fastq.gz"
CHIP_FASTQ_R2="$DATA_DIR/H3K4me3_ChIP_R2.fastq.gz"
INPUT_FASTQ_R1="$DATA_DIR/Input_control_R1.fastq.gz"
INPUT_FASTQ_R2="$DATA_DIR/Input_control_R2.fastq.gz"
MACS2_QVALUE=0.05

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# ============================================================================
# 日志函数
# ============================================================================

log_info() {
    echo -e "${GREEN}[INFO]${NC} $1" | tee -a "$LOGS_DIR/pipeline.log"
}

log_section() {
    echo -e "${BLUE}========================================${NC}" | tee -a "$LOGS_DIR/pipeline.log"
    echo -e "${BLUE}$1${NC}" | tee -a "$LOGS_DIR/pipeline.log"
    echo -e "${BLUE}========================================${NC}" | tee -a "$LOGS_DIR/pipeline.log"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1" | tee -a "$LOGS_DIR/pipeline.log"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1" | tee -a "$LOGS_DIR/pipeline.log"
}

# ============================================================================
# 步骤1：质量控制
# ============================================================================

step_quality_control() {
    log_section "步骤1：数据质量控制"

    if [ ! -f "$CHIP_FASTQ_R1" ] || [ ! -f "$INPUT_FASTQ_R1" ]; then
        log_error "未找到FASTQ文件"
        return 1
    fi

    log_info "运行FastQC..."
    if command -v fastqc &> /dev/null; then
        fastqc "$CHIP_FASTQ_R1" "$CHIP_FASTQ_R2" "$INPUT_FASTQ_R1" "$INPUT_FASTQ_R2" \
            -o "$QC_DIR" -t "$THREADS" 2>&1 | tee -a "$LOGS_DIR/fastqc.log"
        log_success "FastQC分析完成"
    else
        log_error "FastQC未安装"
    fi
}

# ============================================================================
# 步骤2：序列比对
# ============================================================================

step_alignment() {
    log_section "步骤2：序列比对（BWA）"

    if [ ! -f "$REFERENCE" ]; then
        log_error "参考基因组文件不存在"
        return 1
    fi

    # 建立索引
    log_info "检查参考基因组索引..."
    if [ ! -f "${REFERENCE}.bwt" ]; then
        log_info "构建BWA索引..."
        if ! bwa index "$REFERENCE" 2>&1 | tee -a "$LOGS_DIR/bwa_index.log"; then
            log_error "索引构建失败"
            return 1
        fi
    fi

    # ChIP样本比对
    log_info "比对ChIP样本..."
    CHIP_SAM="$RESULTS_DIR/alignment/H3K4me3_ChIP.sam"
    bwa mem -t "$THREADS" "$REFERENCE" "$CHIP_FASTQ_R1" "$CHIP_FASTQ_R2" \
        > "$CHIP_SAM" 2>&1 | tee -a "$LOGS_DIR/bwa_chip.log"

    # Input样本比对
    log_info "比对Input对照..."
    INPUT_SAM="$RESULTS_DIR/alignment/Input_control.sam"
    bwa mem -t "$THREADS" "$REFERENCE" "$INPUT_FASTQ_R1" "$INPUT_FASTQ_R2" \
        > "$INPUT_SAM" 2>&1 | tee -a "$LOGS_DIR/bwa_input.log"

    # SAM转BAM并排序
    log_info "转换并排序BAM文件..."
    for sample in H3K4me3_ChIP Input_control; do
        samfile="$RESULTS_DIR/alignment/${sample}.sam"
        bamfile="$RESULTS_DIR/alignment/${sample}.bam"
        
        samtools view -@ "$THREADS" -bS "$samfile" | \
            samtools sort -@ "$THREADS" -o "$bamfile"
        samtools index -@ "$THREADS" "$bamfile"
        rm "$samfile"
        
        log_success "$sample处理完成"
    done
}

# ============================================================================
# 步骤3：Peak calling
# ============================================================================

step_peak_calling() {
    log_section "步骤3：Peak calling（MACS2）"

    CHIP_BAM="$RESULTS_DIR/alignment/H3K4me3_ChIP.bam"
    INPUT_BAM="$RESULTS_DIR/alignment/Input_control.bam"

    if [ ! -f "$CHIP_BAM" ] || [ ! -f "$INPUT_BAM" ]; then
        log_error "BAM文件不存在"
        return 1
    fi

    log_info "运行MACS2..."
    macs2 callpeak \
        -t "$CHIP_BAM" \
        -c "$INPUT_BAM" \
        -n H3K4me3 \
        -f BAM \
        -g hs \
        -q "$MACS2_QVALUE" \
        --outdir "$RESULTS_DIR/peaks" \
        2>&1 | tee -a "$LOGS_DIR/macs2.log"

    if [ -f "$RESULTS_DIR/peaks/H3K4me3_peaks.narrowPeak" ]; then
        peak_count=$(wc -l < "$RESULTS_DIR/peaks/H3K4me3_peaks.narrowPeak")
        log_success "Peak calling完成，检测到 $peak_count 个Peaks"
    else
        log_error "Peak calling失败"
        return 1
    fi
}

# ============================================================================
# 步骤4：可视化
# ============================================================================

step_visualization() {
    log_section "步骤4：生成可视化图表"

    CHIP_BAM="$RESULTS_DIR/alignment/H3K4me3_ChIP.bam"
    INPUT_BAM="$RESULTS_DIR/alignment/Input_control.bam"

    if [ ! -f "$CHIP_BAM" ]; then
        log_error "BAM文件不存在"
        return 1
    fi

    log_info "生成覆盖度文件..."
    if command -v bamCoverage &> /dev/null; then
        bamCoverage -b "$CHIP_BAM" -o "$RESULTS_DIR/plots/H3K4me3.bw" \
            --binSize 50 --numberOfProcessors "$THREADS" \
            2>&1 | tee -a "$LOGS_DIR/deeptools.log"
        
        bamCoverage -b "$INPUT_BAM" -o "$RESULTS_DIR/plots/Input.bw" \
            --binSize 50 --numberOfProcessors "$THREADS" \
            2>&1 | tee -a "$LOGS_DIR/deeptools.log"
        
        log_success "覆盖度文件已生成"
    else
        log_error "deepTools未安装"
    fi
}

# ============================================================================
# 主程序
# ============================================================================

usage() {
    cat << EOF
ChIP-seq分析流程脚本

使用方法: $0 [选项]

选项:
    --qc          仅质量控制
    --align       仅序列比对
    --peaks       仅Peak calling
    --visualize   仅可视化
    --all         完整流程（默认）
    --threads N   线程数（默认4）
    --help        显示帮助

EOF
}

main() {
    local action="all"

    while [[ $# -gt 0 ]]; do
        case $1 in
            --qc) action="qc"; shift ;;
            --align) action="align"; shift ;;
            --peaks) action="peaks"; shift ;;
            --visualize) action="visualize"; shift ;;
            --all) action="all"; shift ;;
            --threads) THREADS="$2"; shift 2 ;;
            --help) usage; exit 0 ;;
            *) echo "未知选项: $1"; usage; exit 1 ;;
        esac
    done

    log_section "ChIP-seq分析流程"
    log_info "开始时间: $(date)"

    case $action in
        qc) step_quality_control ;;
        align) step_alignment ;;
        peaks) step_peak_calling ;;
        visualize) step_visualization ;;
        all)
            step_quality_control
            step_alignment
            step_peak_calling
            step_visualization
            ;;
    esac

    log_info "结束时间: $(date)"
    log_success "分析完成！日志: $LOGS_DIR/pipeline.log"
}

main "$@"
