#!/bin/bash

###############################################################################
#           RNA-seq完整分析流程脚本                                          #
# 集成所有分析步骤：质量控制、比对、定量、差异分析和可视化               #
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
SCRIPTS_DIR="$SCRIPT_DIR"

# 创建必要的目录
mkdir -p "$RESULTS_DIR"/{alignment,quantification,deseq2,plots} "$QC_DIR" "$LOGS_DIR"

# 参数设置
THREADS=4
REFERENCE_FASTA="$REFERENCE_DIR/demo_reference.fa"
REFERENCE_GTF="$REFERENCE_DIR/demo_annotation.gtf"
HISAT2_INDEX="$REFERENCE_DIR/demo_index"
PADJ_CUTOFF=0.05
FC_CUTOFF=1.0

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# ============================================================================
# 日志和输出函数
# ============================================================================

log_info() {
    echo -e "${GREEN}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} ${BLUE}[INFO]${NC} $1" | tee -a "$LOGS_DIR/pipeline.log"
}

log_section() {
    echo "" | tee -a "$LOGS_DIR/pipeline.log"
    echo -e "${BLUE}========================================${NC}" | tee -a "$LOGS_DIR/pipeline.log"
    echo -e "${BLUE}$1${NC}" | tee -a "$LOGS_DIR/pipeline.log"
    echo -e "${BLUE}========================================${NC}" | tee -a "$LOGS_DIR/pipeline.log"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1" | tee -a "$LOGS_DIR/pipeline.log"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1" | tee -a "$LOGS_DIR/pipeline.log"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1" | tee -a "$LOGS_DIR/pipeline.log"
}

# 打印帮助信息
usage() {
    cat << EOF
RNA-seq完整分析流程脚本

使用方法: $0 [选项]

选项:
    --qc             仅进行质量控制
    --align          仅进行序列比对
    --quantify       仅进行基因定量
    --deseq2         仅进行DESeq2分析
    --visualize      仅生成可视化图表
    --all            完整分析流程（默认）
    --threads N      指定线程数（默认4）
    --help           显示此帮助信息

示例:
    $0                      # 完整分析流程
    $0 --qc                 # 仅进行质量控制
    $0 --align --threads 8  # 使用8个线程进行比对

EOF
}

# ============================================================================
# 步骤1：质量控制
# ============================================================================

step_quality_control() {
    log_section "步骤1：数据质量控制（FastQC + MultiQC）"

    if [ ! -d "$DATA_DIR" ] || [ -z "$(ls -A $DATA_DIR/*.fastq.gz 2>/dev/null)" ]; then
        log_error "未找到FASTQ文件，请先运行数据准备脚本"
        return 1
    fi

    log_info "使用FastQC进行质量评估..."

    # 运行FastQC
    if command -v fastqc &> /dev/null; then
        fastqc "$DATA_DIR"/*.fastq.gz -o "$QC_DIR" -t "$THREADS" 2>&1 | tee -a "$LOGS_DIR/fastqc.log"
        log_success "FastQC分析完成"
    else
        log_warn "FastQC未安装，跳过质量控制"
        return 0
    fi

    # 运行MultiQC
    log_info "使用MultiQC生成汇总报告..."
    if command -v multiqc &> /dev/null; then
        multiqc "$QC_DIR" -o "$QC_DIR/multiqc_report" -f 2>&1 | tee -a "$LOGS_DIR/multiqc.log"
        log_success "MultiQC报告已生成: $QC_DIR/multiqc_report/multiqc_report.html"
    else
        log_warn "MultiQC未安装，跳过汇总报告生成"
    fi

    # 生成QC摘要
    cat > "$LOGS_DIR/qc_summary.txt" << EOF
FastQC质量控制摘要
==================
生成时间: $(date)
FASTQ文件数: $(ls -1 "$DATA_DIR"/*.fastq.gz | wc -l)
QC输出目录: $QC_DIR

检查项目:
  ✓ 序列质量评分
  ✓ GC含量分布
  ✓ 序列重复
  ✓ 接头污染检测
  ✓ 多个QC指标

建议:
  - 查看MultiQC报告: $QC_DIR/multiqc_report/multiqc_report.html
  - 如果质量问题严重，可使用Trimmomatic进行数据清洗

EOF

    log_info "QC摘要已保存到: $LOGS_DIR/qc_summary.txt"
}

# ============================================================================
# 步骤2：参考基因组索引构建
# ============================================================================

step_build_index() {
    log_section "步骤2：构建HISAT2索引"

    if [ ! -f "$REFERENCE_FASTA" ]; then
        log_error "参考基因组文件不存在: $REFERENCE_FASTA"
        return 1
    fi

    if [ -f "${HISAT2_INDEX}.1.ht2" ]; then
        log_info "索引文件已存在，跳过构建"
        return 0
    fi

    log_info "构建HISAT2索引..."
    if command -v hisat2-build &> /dev/null; then
        hisat2-build -p "$THREADS" "$REFERENCE_FASTA" "$HISAT2_INDEX" \
            2>&1 | tee -a "$LOGS_DIR/hisat2_build.log"

        if [ -f "${HISAT2_INDEX}.1.ht2" ]; then
            log_success "HISAT2索引构建成功"
        else
            log_error "HISAT2索引构建失败"
            return 1
        fi
    else
        log_error "HISAT2未安装"
        return 1
    fi
}

# ============================================================================
# 步骤3：序列比对
# ============================================================================

step_alignment() {
    log_section "步骤3：HISAT2序列比对"

    if [ ! -f "${HISAT2_INDEX}.1.ht2" ]; then
        log_info "索引文件不存在，先构建索引..."
        step_build_index || return 1
    fi

    if [ ! -f "$REFERENCE_GTF" ]; then
        log_warn "GTF文件不存在: $REFERENCE_GTF"
    fi

    # 获取所有样本
    samples=$(ls "$DATA_DIR"/*_R1.fastq.gz | xargs -I {} basename {} _R1.fastq.gz | sort -u)

    local sample_count=0
    local align_success=0

    for sample in $samples; do
        sample_count=$((sample_count + 1))

        log_info "[$sample_count/$(echo "$samples" | wc -w)] 处理样本: $sample"

        read1="$DATA_DIR/${sample}_R1.fastq.gz"
        read2="$DATA_DIR/${sample}_R2.fastq.gz"
        samfile="$RESULTS_DIR/alignment/${sample}.sam"
        bamfile="$RESULTS_DIR/alignment/${sample}_sorted.bam"

        if [ ! -f "$read1" ] || [ ! -f "$read2" ]; then
            log_warn "样本文件不完整: $sample"
            continue
        fi

        if [ -f "$bamfile" ]; then
            log_info "✓ 样本$sample已比对，跳过"
            align_success=$((align_success + 1))
            continue
        fi

        # 运行HISAT2比对
        log_info "  运行HISAT2比对..."
        if ! hisat2 -x "$HISAT2_INDEX" \
            -1 "$read1" \
            -2 "$read2" \
            -S "$samfile" \
            --threads "$THREADS" \
            --rna-strandness RF \
            --summary-file "$LOGS_DIR/${sample}_alignment_summary.txt" \
            2>&1 | tee -a "$LOGS_DIR/alignment.log"; then
            log_warn "  样本$sample比对失败"
            continue
        fi

        # SAM转BAM并排序
        log_info "  转换并排序BAM文件..."
        if ! samtools view -@ "$THREADS" -bS "$samfile" | \
            samtools sort -@ "$THREADS" -o "$bamfile" 2>&1 | tee -a "$LOGS_DIR/samtools.log"; then
            log_warn "  样本$sample BAM转换失败"
            rm "$samfile"
            continue
        fi

        # 建立BAM索引
        log_info "  建立BAM索引..."
        if samtools index -@ "$THREADS" "$bamfile" 2>&1; then
            log_success "  ✓ 样本$sample比对完成"
            rm "$samfile"
            align_success=$((align_success + 1))
        else
            log_warn "  样本$sample索引建立失败"
        fi
    done

    log_success "序列比对完成: $align_success/$sample_count 样本成功"

    # 生成比对统计摘要
    cat > "$LOGS_DIR/alignment_summary.txt" << EOF
HISAT2序列比对摘要
==================
比对完成时间: $(date)
参考基因组: $(basename "$REFERENCE_FASTA")
样本总数: $sample_count
成功比对: $align_success

比对统计详情:
EOF

    for summary_file in "$LOGS_DIR"/*_alignment_summary.txt; do
        if [ -f "$summary_file" ]; then
            echo "" >> "$LOGS_DIR/alignment_summary.txt"
            cat "$summary_file" >> "$LOGS_DIR/alignment_summary.txt"
        fi
    done

    log_info "比对摘要已保存到: $LOGS_DIR/alignment_summary.txt"
}

# ============================================================================
# 步骤4：基因表达定量
# ============================================================================

step_quantification() {
    log_section "步骤4：使用featureCounts进行基因定量"

    if [ ! -f "$REFERENCE_GTF" ]; then
        log_error "GTF注释文件不存在: $REFERENCE_GTF"
        return 1
    fi

    bam_files=$(ls "$RESULTS_DIR/alignment"/*_sorted.bam 2>/dev/null)
    if [ -z "$bam_files" ]; then
        log_error "未找到BAM文件，请先进行序列比对"
        return 1
    fi

    log_info "使用featureCounts统计基因表达量..."

    if command -v featureCounts &> /dev/null; then
        featureCounts \
            -a "$REFERENCE_GTF" \
            -o "$RESULTS_DIR/quantification/gene_counts.txt" \
            -T "$THREADS" \
            -p \
            -B \
            -C \
            -s 2 \
            $bam_files \
            2>&1 | tee -a "$LOGS_DIR/featurecounts.log"

        if [ -f "$RESULTS_DIR/quantification/gene_counts.txt" ]; then
            log_success "基因定量完成"

            # 查看计数结果
            log_info "计数矩阵摘要:"
            echo "---"
            head -5 "$RESULTS_DIR/quantification/gene_counts.txt"
            echo "..."
            echo "---"
            echo "总基因数: $(tail -n +2 "$RESULTS_DIR/quantification/gene_counts.txt" | wc -l)"

            # 查看计数统计
            if [ -f "$RESULTS_DIR/quantification/gene_counts.txt.summary" ]; then
                log_info "计数统计:"
                cat "$RESULTS_DIR/quantification/gene_counts.txt.summary"
            fi
        else
            log_error "基因定量失败"
            return 1
        fi
    else
        log_error "featureCounts未安装"
        return 1
    fi
}

# ============================================================================
# 步骤5：DESeq2差异表达分析
# ============================================================================

step_deseq2_analysis() {
    log_section "步骤5：DESeq2差异表达分析"

    if [ ! -f "$RESULTS_DIR/quantification/gene_counts.txt" ]; then
        log_error "基因计数矩阵不存在，请先进行基因定量"
        return 1
    fi

    log_info "创建DESeq2分析脚本..."

    cat > "$SCRIPTS_DIR/deseq2_main.R" << 'R_EOF'
#!/usr/bin/env Rscript

# 加载必要的包
suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
})

# 配置参数
setwd(Sys.getenv("WORK_DIR"))
counts_file <- Sys.getenv("COUNTS_FILE")
gtf_file <- Sys.getenv("GTF_FILE")
results_dir <- Sys.getenv("RESULTS_DIR")
padj_cutoff <- as.numeric(Sys.getenv("PADJ_CUTOFF", "0.05"))
fc_cutoff <- as.numeric(Sys.getenv("FC_CUTOFF", "1.0"))

cat("=== DESeq2分析 ===\n")
cat("计数矩阵:", counts_file, "\n")
cat("结果输出:", results_dir, "\n")
cat("显著性阈值 (padj <", padj_cutoff, ")\n")
cat("折叠变化阈值 (|log2FC| >", fc_cutoff, ")\n\n")

# 读取计数矩阵
cat("读取计数矩阵...\n")
countData <- read.table(counts_file, header=TRUE, row.names=1, skip=1)

# 提取计数列（去除注释列）
# featureCounts的前6列是注释信息
countData <- countData[,6:ncol(countData)]

# 修改列名（去除路径和扩展名）
colnames(countData) <- gsub(".*results.alignment.|_sorted.bam", "", colnames(countData))

# 样本分组信息
sampleInfo <- data.frame(
    row.names = colnames(countData),
    group = factor(c(
        rep("treatment", length(grep("^sample", colnames(countData)))),
        rep("control", length(grep("^ctrl", colnames(countData))))
    )),
    batch = factor(rep(1:3, 2))
)

cat("\n样本信息:\n")
print(sampleInfo)

# 创建DESeq2数据集
cat("\n创建DESeq2数据集...\n")
dds <- DESeqDataSetFromMatrix(
    countData = round(countData),
    colData = sampleInfo,
    design = ~ group
)

# 过滤低表达基因
cat("过滤低表达基因 (count sum >= 10)...\n")
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
cat("保留基因数:", nrow(dds), "\n\n")

# 运行DESeq2分析
cat("运行DESeq2分析...\n")
dds <- DESeq(dds)

# 获取结果
res <- results(dds, contrast=c("group", "treatment", "control"))

# 查看结果摘要
cat("\n=== DESeq2分析结果摘要 ===\n")
summary(res)

# 保存完整结果
cat("\n保存结果文件...\n")
write.csv(as.data.frame(res),
    file.path(results_dir, "deseq2", "deseq2_results.csv"))

# 识别显著差异基因
sig_genes <- subset(res, padj < padj_cutoff & abs(log2FoldChange) > fc_cutoff)
write.csv(as.data.frame(sig_genes),
    file.path(results_dir, "deseq2", "significant_genes.csv"))

cat("显著差异基因数:", nrow(sig_genes), "\n")
cat("  上调基因:", sum(sig_genes$log2FoldChange > 0), "\n")
cat("  下调基因:", sum(sig_genes$log2FoldChange < 0), "\n")

# 保存标准化计数
vsd <- vst(dds, blind=FALSE)
write.csv(assay(vsd),
    file.path(results_dir, "deseq2", "normalized_counts.csv"))

# 保存分析信息以供后续使用
saveRDS(dds, file.path(results_dir, "deseq2", "dds_object.rds"))
saveRDS(res, file.path(results_dir, "deseq2", "results_object.rds"))
saveRDS(vsd, file.path(results_dir, "deseq2", "vsd_object.rds"))
saveRDS(sampleInfo, file.path(results_dir, "deseq2", "sample_info.rds"))

cat("\n✓ DESeq2分析完成\n")
R_EOF

    # 运行DESeq2脚本
    log_info "运行DESeq2分析脚本..."

    if ! WORK_DIR="$WORK_DIR" \
         COUNTS_FILE="$RESULTS_DIR/quantification/gene_counts.txt" \
         GTF_FILE="$REFERENCE_GTF" \
         RESULTS_DIR="$RESULTS_DIR" \
         PADJ_CUTOFF="$PADJ_CUTOFF" \
         FC_CUTOFF="$FC_CUTOFF" \
         Rscript "$SCRIPTS_DIR/deseq2_main.R" 2>&1 | tee -a "$LOGS_DIR/deseq2.log"; then
        log_error "DESeq2分析失败"
        return 1
    fi

    if [ -f "$RESULTS_DIR/deseq2/significant_genes.csv" ]; then
        log_success "DESeq2分析完成"
    else
        log_error "DESeq2分析未生成结果文件"
        return 1
    fi
}

# ============================================================================
# 步骤6：结果可视化
# ============================================================================

step_visualization() {
    log_section "步骤6：生成可视化图表"

    if [ ! -f "$RESULTS_DIR/deseq2/dds_object.rds" ]; then
        log_error "DESeq2对象不存在，请先运行DESeq2分析"
        return 1
    fi

    log_info "创建可视化脚本..."

    cat > "$SCRIPTS_DIR/visualization.R" << 'R_EOF'
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
})

results_dir <- Sys.getenv("RESULTS_DIR")
plots_dir <- file.path(results_dir, "plots")

cat("=== 生成可视化图表 ===\n")
cat("图表输出目录:", plots_dir, "\n\n")

# 加载保存的对象
cat("加载分析对象...\n")
dds <- readRDS(file.path(results_dir, "deseq2", "dds_object.rds"))
res <- readRDS(file.path(results_dir, "deseq2", "results_object.rds"))
vsd <- readRDS(file.path(results_dir, "deseq2", "vsd_object.rds"))
sampleInfo <- readRDS(file.path(results_dir, "deseq2", "sample_info.rds"))
sig_genes <- read.csv(file.path(results_dir, "deseq2", "significant_genes.csv"), row.names=1)

# 1. MA图
cat("生成MA图...\n")
pdf(file.path(plots_dir, "MA_plot.pdf"), width=10, height=8)
plotMA(res, main="MA Plot: Treatment vs Control", ylim=c(-5,5))
dev.off()

# 2. PCA图
cat("生成PCA图...\n")
pdf(file.path(plots_dir, "PCA_plot.pdf"), width=10, height=8)
pcaData <- plotPCA(vsd, intgroup="group", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p <- ggplot(pcaData, aes(x = PC1, y = PC2, color = group)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle("PCA Plot") +
    theme_minimal()
print(p)
dev.off()

# 3. 样本距离热图
cat("生成样本距离热图...\n")
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$group, rownames(sampleInfo), sep="-")
colnames(sampleDistMatrix) <- NULL

pdf(file.path(plots_dir, "sample_distance_heatmap.pdf"), width=8, height=8)
pheatmap(sampleDistMatrix,
    clustering_distance_rows=sampleDists,
    clustering_distance_cols=sampleDists,
    main="Sample-to-Sample Distances",
    color=colorRampPalette(c("white", "blue"))(50))
dev.off()

# 4. 热图 - 显著差异基因
if(nrow(sig_genes) > 0) {
    cat("生成差异基因热图...\n")

    # 选择前20个差异基因
    top_genes <- rownames(sig_genes)[order(sig_genes$padj)][1:min(20, nrow(sig_genes))]

    pdf(file.path(plots_dir, "heatmap_top_genes.pdf"), width=10, height=12)
    pheatmap(assay(vsd)[top_genes,],
        cluster_rows=TRUE,
        show_rownames=TRUE,
        cluster_cols=TRUE,
        annotation_col=sampleInfo[,-3],
        scale="row",
        main=paste("Top", min(20, nrow(sig_genes)), "Differentially Expressed Genes"))
    dev.off()
}

cat("\n✓ 可视化图表生成完成\n")
R_EOF

    # 运行可视化脚本
    log_info "运行可视化脚本..."

    if ! RESULTS_DIR="$RESULTS_DIR" Rscript "$SCRIPTS_DIR/visualization.R" 2>&1 | tee -a "$LOGS_DIR/visualization.log"; then
        log_warn "可视化脚本运行出现错误"
    else
        log_success "可视化图表生成完成"
    fi

    # 列出生成的图表
    if [ -d "$RESULTS_DIR/plots" ]; then
        log_info "生成的图表文件:"
        ls -lh "$RESULTS_DIR/plots/"/*.pdf 2>/dev/null | awk '{print "  " $9, "(" $5 ")"}'
    fi
}

# ============================================================================
# 主程序
# ============================================================================

main() {
    local action="all"

    # 解析命令行参数
    while [[ $# -gt 0 ]]; do
        case $1 in
            --qc)
                action="qc"
                shift
                ;;
            --align)
                action="align"
                shift
                ;;
            --quantify)
                action="quantify"
                shift
                ;;
            --deseq2)
                action="deseq2"
                shift
                ;;
            --visualize)
                action="visualize"
                shift
                ;;
            --all)
                action="all"
                shift
                ;;
            --threads)
                THREADS="$2"
                shift 2
                ;;
            --help)
                usage
                exit 0
                ;;
            *)
                echo "未知选项: $1"
                usage
                exit 1
                ;;
        esac
    done

    log_section "RNA-seq完整分析流程"
    log_info "开始时间: $(date)"
    log_info "工作目录: $WORK_DIR"
    log_info "使用线程数: $THREADS"

    # 执行相应步骤
    case $action in
        qc)
            step_quality_control
            ;;
        align)
            step_build_index && step_alignment
            ;;
        quantify)
            step_quantification
            ;;
        deseq2)
            step_deseq2_analysis
            ;;
        visualize)
            step_visualization
            ;;
        all)
            step_quality_control
            step_build_index
            step_alignment
            step_quantification
            step_deseq2_analysis
            step_visualization
            ;;
    esac

    log_info "结束时间: $(date)"
    log_section "分析完成"
    log_success "详细日志已保存到: $LOGS_DIR/pipeline.log"
}

# 执行主程序
main "$@"
