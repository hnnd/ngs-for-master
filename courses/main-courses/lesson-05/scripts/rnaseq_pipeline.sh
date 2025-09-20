#!/bin/bash

# RNA-seq数据分析流水线脚本
# 课程：高通量测序数据分析 - 第5次课
# 作者：王运生
# 邮箱：wangys@hunau.edu.cn

set -e  # 遇到错误立即退出

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 日志函数
log_info() {
    echo -e "${GREEN}[INFO]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

# 配置参数
WORK_DIR="$HOME/ngs-analysis/lesson-05"
THREADS=4
MEMORY="8G"

# 软件路径（如果不在PATH中，请修改为完整路径）
FASTQC="fastqc"
MULTIQC="multiqc"
TRIMMOMATIC="trimmomatic"
HISAT2="hisat2"
HISAT2_BUILD="hisat2-build"
SAMTOOLS="samtools"
FEATURECOUNTS="featureCounts"

# 文件路径
REFERENCE_DIR="$WORK_DIR/reference"
DATA_DIR="$WORK_DIR/data"
RESULTS_DIR="$WORK_DIR/results"
SCRIPTS_DIR="$WORK_DIR/scripts"
LOGS_DIR="$WORK_DIR/logs"

# 样本列表
SAMPLES=("sample1" "sample2" "sample3" "ctrl1" "ctrl2" "ctrl3")

# 帮助信息
show_help() {
    cat << EOF
RNA-seq数据分析流水线

用法: $0 [选项] [步骤]

选项:
    -h, --help          显示帮助信息
    -t, --threads NUM   设置线程数 (默认: 4)
    -w, --workdir DIR   设置工作目录 (默认: $HOME/ngs-analysis/lesson-05)
    -s, --samples LIST  设置样本列表，用逗号分隔 (默认: sample1,sample2,sample3,ctrl1,ctrl2,ctrl3)

步骤:
    setup              创建目录结构
    qc                 质量控制
    trim               数据清洗
    index              建立基因组索引
    align              序列比对
    quantify           基因定量
    all                运行所有步骤

示例:
    $0 setup                    # 只创建目录结构
    $0 -t 8 align              # 使用8线程进行比对
    $0 -s sample1,sample2 qc   # 只对指定样本进行质量控制
    $0 all                     # 运行完整流水线

EOF
}

# 检查软件依赖
check_dependencies() {
    log_info "检查软件依赖..."
    
    local missing_tools=()
    
    for tool in $FASTQC $MULTIQC $TRIMMOMATIC $HISAT2 $SAMTOOLS $FEATURECOUNTS; do
        if ! command -v $tool &> /dev/null; then
            missing_tools+=($tool)
        fi
    done
    
    if [ ${#missing_tools[@]} -ne 0 ]; then
        log_error "以下软件未找到: ${missing_tools[*]}"
        log_error "请安装缺失的软件后重试"
        exit 1
    fi
    
    log_info "所有依赖软件检查通过"
}

# 创建目录结构
setup_directories() {
    log_info "创建目录结构..."
    
    mkdir -p "$WORK_DIR"/{data,reference,scripts,results,logs,qc}
    mkdir -p "$RESULTS_DIR"/{alignment,quantification,deseq2,plots}
    mkdir -p "$DATA_DIR/clean"
    
    log_info "目录结构创建完成"
}

# 质量控制
run_quality_control() {
    log_info "开始质量控制分析..."
    
    # 检查输入文件
    local input_files=()
    for sample in "${SAMPLES[@]}"; do
        for read in R1 R2; do
            local file="$DATA_DIR/${sample}_${read}.fastq.gz"
            if [[ -f "$file" ]]; then
                input_files+=("$file")
            else
                log_warn "文件不存在: $file"
            fi
        done
    done
    
    if [[ ${#input_files[@]} -eq 0 ]]; then
        log_error "未找到输入的FASTQ文件"
        return 1
    fi
    
    # 运行FastQC
    log_info "运行FastQC..."
    $FASTQC "${input_files[@]}" -o "$WORK_DIR/qc" -t $THREADS 2>&1 | tee "$LOGS_DIR/fastqc.log"
    
    # 运行MultiQC
    log_info "运行MultiQC..."
    $MULTIQC "$WORK_DIR/qc" -o "$WORK_DIR/qc/multiqc_report" 2>&1 | tee "$LOGS_DIR/multiqc.log"
    
    log_info "质量控制分析完成"
}

# 数据清洗
run_trimming() {
    log_info "开始数据清洗..."
    
    for sample in "${SAMPLES[@]}"; do
        local r1_input="$DATA_DIR/${sample}_R1.fastq.gz"
        local r2_input="$DATA_DIR/${sample}_R2.fastq.gz"
        
        if [[ ! -f "$r1_input" ]] || [[ ! -f "$r2_input" ]]; then
            log_warn "跳过样本 $sample: 输入文件不存在"
            continue
        fi
        
        log_info "处理样本: $sample"
        
        $TRIMMOMATIC PE -threads $THREADS \
            "$r1_input" "$r2_input" \
            "$DATA_DIR/clean/${sample}_R1_clean.fastq.gz" "$DATA_DIR/clean/${sample}_R1_unpaired.fastq.gz" \
            "$DATA_DIR/clean/${sample}_R2_clean.fastq.gz" "$DATA_DIR/clean/${sample}_R2_unpaired.fastq.gz" \
            ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
            2>&1 | tee "$LOGS_DIR/trimmomatic_${sample}.log"
    done
    
    log_info "数据清洗完成"
}

# 建立基因组索引
build_genome_index() {
    log_info "建立基因组索引..."
    
    local genome_file="$REFERENCE_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    local index_prefix="$REFERENCE_DIR/grch38_index"
    
    if [[ ! -f "$genome_file" ]]; then
        log_error "参考基因组文件不存在: $genome_file"
        log_info "请先下载参考基因组文件"
        return 1
    fi
    
    # 检查索引是否已存在
    if [[ -f "${index_prefix}.1.ht2" ]]; then
        log_info "基因组索引已存在，跳过构建"
        return 0
    fi
    
    log_info "构建HISAT2索引（这可能需要较长时间）..."
    $HISAT2_BUILD -p $THREADS "$genome_file" "$index_prefix" 2>&1 | tee "$LOGS_DIR/hisat2_build.log"
    
    log_info "基因组索引构建完成"
}

# 序列比对
run_alignment() {
    log_info "开始序列比对..."
    
    local index_prefix="$REFERENCE_DIR/grch38_index"
    
    # 检查索引文件
    if [[ ! -f "${index_prefix}.1.ht2" ]]; then
        log_error "HISAT2索引不存在，请先运行 build_genome_index"
        return 1
    fi
    
    for sample in "${SAMPLES[@]}"; do
        local r1_input="$DATA_DIR/clean/${sample}_R1_clean.fastq.gz"
        local r2_input="$DATA_DIR/clean/${sample}_R2_clean.fastq.gz"
        local sam_output="$RESULTS_DIR/alignment/${sample}.sam"
        local bam_output="$RESULTS_DIR/alignment/${sample}_sorted.bam"
        
        if [[ ! -f "$r1_input" ]] || [[ ! -f "$r2_input" ]]; then
            log_warn "跳过样本 $sample: 清洗后的文件不存在"
            continue
        fi
        
        # 检查是否已经完成比对
        if [[ -f "$bam_output" ]]; then
            log_info "样本 $sample 已完成比对，跳过"
            continue
        fi
        
        log_info "比对样本: $sample"
        
        # HISAT2比对
        $HISAT2 -x "$index_prefix" \
            -1 "$r1_input" \
            -2 "$r2_input" \
            -S "$sam_output" \
            --threads $THREADS \
            --rna-strandness RF \
            --summary-file "$LOGS_DIR/${sample}_alignment_summary.txt" \
            2>&1 | tee "$LOGS_DIR/hisat2_${sample}.log"
        
        # SAM转BAM并排序
        log_info "转换和排序BAM文件: $sample"
        $SAMTOOLS view -bS "$sam_output" | \
        $SAMTOOLS sort -o "$bam_output" -@ $THREADS
        
        # 建立索引
        $SAMTOOLS index "$bam_output"
        
        # 删除SAM文件节省空间
        rm "$sam_output"
        
        # 生成比对统计
        $SAMTOOLS flagstat "$bam_output" > "$LOGS_DIR/${sample}_flagstat.txt"
    done
    
    log_info "序列比对完成"
}

# 基因定量
run_quantification() {
    log_info "开始基因定量..."
    
    local gtf_file="$REFERENCE_DIR/Homo_sapiens.GRCh38.104.gtf"
    local output_file="$RESULTS_DIR/quantification/gene_counts.txt"
    
    if [[ ! -f "$gtf_file" ]]; then
        log_error "GTF注释文件不存在: $gtf_file"
        return 1
    fi
    
    # 收集所有BAM文件
    local bam_files=()
    for sample in "${SAMPLES[@]}"; do
        local bam_file="$RESULTS_DIR/alignment/${sample}_sorted.bam"
        if [[ -f "$bam_file" ]]; then
            bam_files+=("$bam_file")
        else
            log_warn "BAM文件不存在: $bam_file"
        fi
    done
    
    if [[ ${#bam_files[@]} -eq 0 ]]; then
        log_error "未找到BAM文件"
        return 1
    fi
    
    log_info "使用featureCounts进行基因计数..."
    $FEATURECOUNTS -a "$gtf_file" \
        -o "$output_file" \
        -T $THREADS \
        -p \
        -B \
        -C \
        -s 2 \
        "${bam_files[@]}" \
        2>&1 | tee "$LOGS_DIR/featurecounts.log"
    
    log_info "基因定量完成"
}

# 生成分析报告
generate_report() {
    log_info "生成分析报告..."
    
    local report_file="$WORK_DIR/analysis_report.txt"
    
    cat > "$report_file" << EOF
RNA-seq数据分析报告
==================

分析时间: $(date)
工作目录: $WORK_DIR
线程数: $THREADS

样本信息:
EOF
    
    for sample in "${SAMPLES[@]}"; do
        echo "- $sample" >> "$report_file"
    done
    
    echo "" >> "$report_file"
    echo "文件统计:" >> "$report_file"
    
    # 统计各类文件数量
    local fastq_count=$(find "$DATA_DIR" -name "*.fastq.gz" | wc -l)
    local bam_count=$(find "$RESULTS_DIR/alignment" -name "*.bam" | wc -l)
    local qc_count=$(find "$WORK_DIR/qc" -name "*_fastqc.html" | wc -l)
    
    echo "- FASTQ文件: $fastq_count" >> "$report_file"
    echo "- BAM文件: $bam_count" >> "$report_file"
    echo "- QC报告: $qc_count" >> "$report_file"
    
    # 如果存在基因计数文件，统计基因数量
    local count_file="$RESULTS_DIR/quantification/gene_counts.txt"
    if [[ -f "$count_file" ]]; then
        local gene_count=$(tail -n +2 "$count_file" | wc -l)
        echo "- 检测基因数: $gene_count" >> "$report_file"
    fi
    
    echo "" >> "$report_file"
    echo "日志文件位置: $LOGS_DIR" >> "$report_file"
    echo "结果文件位置: $RESULTS_DIR" >> "$report_file"
    
    log_info "分析报告已生成: $report_file"
}

# 主函数
main() {
    local step="all"
    
    # 解析命令行参数
    while [[ $# -gt 0 ]]; do
        case $1 in
            -h|--help)
                show_help
                exit 0
                ;;
            -t|--threads)
                THREADS="$2"
                shift 2
                ;;
            -w|--workdir)
                WORK_DIR="$2"
                shift 2
                ;;
            -s|--samples)
                IFS=',' read -ra SAMPLES <<< "$2"
                shift 2
                ;;
            setup|qc|trim|index|align|quantify|all)
                step="$1"
                shift
                ;;
            *)
                log_error "未知参数: $1"
                show_help
                exit 1
                ;;
        esac
    done
    
    # 更新路径变量
    REFERENCE_DIR="$WORK_DIR/reference"
    DATA_DIR="$WORK_DIR/data"
    RESULTS_DIR="$WORK_DIR/results"
    SCRIPTS_DIR="$WORK_DIR/scripts"
    LOGS_DIR="$WORK_DIR/logs"
    
    log_info "开始RNA-seq数据分析流水线"
    log_info "工作目录: $WORK_DIR"
    log_info "线程数: $THREADS"
    log_info "样本: ${SAMPLES[*]}"
    log_info "执行步骤: $step"
    
    # 检查依赖
    check_dependencies
    
    # 根据步骤执行相应功能
    case $step in
        setup)
            setup_directories
            ;;
        qc)
            setup_directories
            run_quality_control
            ;;
        trim)
            setup_directories
            run_trimming
            ;;
        index)
            setup_directories
            build_genome_index
            ;;
        align)
            setup_directories
            build_genome_index
            run_alignment
            ;;
        quantify)
            setup_directories
            run_quantification
            ;;
        all)
            setup_directories
            run_quality_control
            run_trimming
            build_genome_index
            run_alignment
            run_quantification
            generate_report
            ;;
        *)
            log_error "未知步骤: $step"
            show_help
            exit 1
            ;;
    esac
    
    log_info "RNA-seq数据分析流水线完成"
}

# 运行主函数
main "$@"