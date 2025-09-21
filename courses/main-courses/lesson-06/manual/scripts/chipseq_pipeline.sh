#!/bin/bash
# ChIP-seq完整分析流程脚本
# 作者：王运生
# 日期：2025-01-21
# 用法：bash chipseq_pipeline.sh

set -e  # 遇到错误立即退出
set -u  # 使用未定义变量时报错

# 加载配置
source ~/ngs-analysis/lesson-06/scripts/config.sh

echo "=== ChIP-seq完整分析流程 ==="
echo "工作目录：$WORK_DIR"
echo "开始时间：$(date)"

# 创建日志文件
LOG_FILE="$WORK_DIR/logs/pipeline_$(date +%Y%m%d_%H%M%S).log"
mkdir -p $(dirname $LOG_FILE)

# 日志函数
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a $LOG_FILE
}

# 错误处理函数
error_exit() {
    log "错误：$1"
    exit 1
}

# 检查必要文件
check_files() {
    log "检查输入文件..."
    
    local required_files=(
        "$DATA_DIR/H3K4me3_ChIP.fastq.gz"
        "$DATA_DIR/Input_control.fastq.gz"
        "$DATA_DIR/hg38.fa"
        "$DATA_DIR/gencode.v38.gtf"
    )
    
    for file in "${required_files[@]}"; do
        if [ ! -f "$file" ]; then
            error_exit "缺少必要文件：$file"
        fi
        log "✓ 找到文件：$(basename $file)"
    done
}

# 步骤1：质量控制
quality_control() {
    log "=== 步骤1：质量控制 ==="
    
    local fastqc_dir="$RESULTS_DIR/fastqc"
    mkdir -p $fastqc_dir
    
    # FastQC分析
    log "运行FastQC..."
    fastqc $DATA_DIR/H3K4me3_ChIP.fastq.gz -o $fastqc_dir --threads $THREADS
    fastqc $DATA_DIR/Input_control.fastq.gz -o $fastqc_dir --threads $THREADS
    
    log "✓ FastQC分析完成"
}

# 步骤2：序列比对
alignment() {
    log "=== 步骤2：序列比对 ==="
    
    local align_dir="$RESULTS_DIR/alignment"
    mkdir -p $align_dir
    
    # 检查BWA索引
    if [ ! -f "$DATA_DIR/hg38.fa.bwt" ]; then
        log "构建BWA索引..."
        bwa index $DATA_DIR/hg38.fa
    fi
    
    # 比对ChIP样本
    log "比对ChIP样本..."
    bwa mem -t $THREADS $DATA_DIR/hg38.fa $DATA_DIR/H3K4me3_ChIP.fastq.gz | \
        samtools sort -@ $THREADS -o $align_dir/H3K4me3_ChIP.bam
    samtools index $align_dir/H3K4me3_ChIP.bam
    
    # 比对Input对照
    log "比对Input对照..."
    bwa mem -t $THREADS $DATA_DIR/hg38.fa $DATA_DIR/Input_control.fastq.gz | \
        samtools sort -@ $THREADS -o $align_dir/Input_control.bam
    samtools index $align_dir/Input_control.bam
    
    # 生成比对统计
    log "生成比对统计..."
    samtools flagstat $align_dir/H3K4me3_ChIP.bam > $align_dir/H3K4me3_ChIP.flagstat
    samtools flagstat $align_dir/Input_control.bam > $align_dir/Input_control.flagstat
    
    log "✓ 序列比对完成"
}

# 步骤3：Peak calling
peak_calling() {
    log "=== 步骤3：Peak calling ==="
    
    local peaks_dir="$RESULTS_DIR/peaks"
    mkdir -p $peaks_dir
    
    # MACS2 Peak calling
    log "运行MACS2..."
    macs2 callpeak \
        -t $RESULTS_DIR/alignment/H3K4me3_ChIP.bam \
        -c $RESULTS_DIR/alignment/Input_control.bam \
        -f BAM \
        -g $MACS2_GSIZE \
        -n $OUTPUT_PREFIX \
        -q $MACS2_QVALUE \
        --outdir $peaks_dir \
        --call-summits \
        2>&1 | tee $peaks_dir/macs2.log
    
    # 统计Peak数量
    local peak_count=$(wc -l < $peaks_dir/${OUTPUT_PREFIX}_peaks.narrowPeak)
    log "检测到 $peak_count 个peaks"
    
    # 计算FRiP score
    log "计算FRiP score..."
    local total_reads=$(samtools view -c $RESULTS_DIR/alignment/H3K4me3_ChIP.bam)
    local peak_reads=$(bedtools intersect -a $RESULTS_DIR/alignment/H3K4me3_ChIP.bam \
        -b $peaks_dir/${OUTPUT_PREFIX}_peaks.narrowPeak -u | samtools view -c)
    local frip=$(echo "scale=4; $peak_reads / $total_reads" | bc -l)
    
    log "FRiP Score: $frip"
    echo "FRiP Score: $frip" > $peaks_dir/frip_score.txt
    
    log "✓ Peak calling完成"
}

# 步骤4：数据可视化
visualization() {
    log "=== 步骤4：数据可视化 ==="
    
    local viz_dir="$WORK_DIR/figures"
    mkdir -p $viz_dir
    
    # 生成BigWig文件
    log "生成BigWig文件..."
    bamCoverage -b $RESULTS_DIR/alignment/H3K4me3_ChIP.bam \
        -o $RESULTS_DIR/H3K4me3_ChIP.bw \
        --normalizeUsing RPKM \
        --binSize 10 \
        --numberOfProcessors $THREADS
    
    bamCoverage -b $RESULTS_DIR/alignment/Input_control.bam \
        -o $RESULTS_DIR/Input_control.bw \
        --normalizeUsing RPKM \
        --binSize 10 \
        --numberOfProcessors $THREADS
    
    # 生成热图
    log "生成热图..."
    computeMatrix reference-point \
        -S $RESULTS_DIR/H3K4me3_ChIP.bw \
        -R $RESULTS_DIR/peaks/${OUTPUT_PREFIX}_summits.bed \
        --referencePoint center \
        -b 2000 -a 2000 \
        --skipZeros \
        -o $RESULTS_DIR/${OUTPUT_PREFIX}_matrix.gz \
        --numberOfProcessors $THREADS
    
    plotHeatmap -m $RESULTS_DIR/${OUTPUT_PREFIX}_matrix.gz \
        -out $viz_dir/${OUTPUT_PREFIX}_heatmap.png \
        --colorMap Blues \
        --whatToShow 'heatmap and colorbar' \
        --plotTitle "H3K4me3 Signal Heatmap"
    
    # 生成Profile图
    log "生成Profile图..."
    plotProfile -m $RESULTS_DIR/${OUTPUT_PREFIX}_matrix.gz \
        -out $viz_dir/${OUTPUT_PREFIX}_profile.png \
        --numPlotsPerRow 1 \
        --plotTitle "H3K4me3 Signal Profile"
    
    log "✓ 数据可视化完成"
}

# 步骤5：Peak注释
peak_annotation() {
    log "=== 步骤5：Peak注释 ==="
    
    local anno_dir="$RESULTS_DIR/annotation"
    mkdir -p $anno_dir
    
    # 运行R脚本进行注释
    log "运行Peak注释..."
    R --slave --no-restore --no-save << EOF
# 加载必要的包
suppressMessages({
    library(ChIPseeker)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(ggplot2)
})

# 设置工作目录
setwd("$WORK_DIR")

# 读取Peak文件
peaks <- readPeakFile("$RESULTS_DIR/peaks/${OUTPUT_PREFIX}_peaks.narrowPeak")
cat("读取到", length(peaks), "个peaks\n")

# 获取基因组注释
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Peak注释
peakAnno <- annotatePeak(peaks, tssRegion=c(-2000, 2000), TxDb=txdb)
cat("Peak注释完成\n")

# 保存注释结果
write.csv(as.data.frame(peakAnno), "$anno_dir/peak_annotation.csv", row.names=FALSE)

# 可视化Peak基因组分布
pdf("$viz_dir/peak_genomic_distribution.pdf", width=12, height=8)
plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
plotDistToTSS(peakAnno)
dev.off()

# 获取Peak相关基因
genes <- as.data.frame(peakAnno)\$geneId
genes <- genes[!is.na(genes)]
genes <- unique(genes)
cat("获得", length(genes), "个相关基因\n")

# GO富集分析
tryCatch({
    ego <- enrichGO(gene = genes,
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    readable = TRUE)
    
    if (nrow(as.data.frame(ego)) > 0) {
        write.csv(as.data.frame(ego), "$anno_dir/GO_enrichment.csv", row.names=FALSE)
        
        pdf("$viz_dir/GO_enrichment.pdf", width=12, height=8)
        print(dotplot(ego, showCategory=20))
        print(barplot(ego, showCategory=20))
        dev.off()
        cat("GO富集分析完成\n")
    } else {
        cat("GO富集分析无显著结果\n")
    }
}, error = function(e) {
    cat("GO富集分析出错:", e\$message, "\n")
})

# KEGG通路分析
tryCatch({
    kk <- enrichKEGG(gene = genes,
                     organism = 'hsa',
                     pvalueCutoff = 0.05)
    
    if (nrow(as.data.frame(kk)) > 0) {
        write.csv(as.data.frame(kk), "$anno_dir/KEGG_enrichment.csv", row.names=FALSE)
        
        pdf("$viz_dir/KEGG_enrichment.pdf", width=12, height=8)
        print(dotplot(kk, showCategory=20))
        dev.off()
        cat("KEGG通路分析完成\n")
    } else {
        cat("KEGG通路分析无显著结果\n")
    }
}, error = function(e) {
    cat("KEGG通路分析出错:", e\$message, "\n")
})

cat("Peak注释和功能分析完成\n")
EOF
    
    log "✓ Peak注释完成"
}

# 步骤6：生成分析报告
generate_report() {
    log "=== 步骤6：生成分析报告 ==="
    
    local report_file="$WORK_DIR/analysis_report.html"
    
    # 收集统计信息
    local chip_total_reads=$(samtools view -c $RESULTS_DIR/alignment/H3K4me3_ChIP.bam)
    local chip_mapped_reads=$(samtools view -c -F 4 $RESULTS_DIR/alignment/H3K4me3_ChIP.bam)
    local mapping_rate=$(echo "scale=2; $chip_mapped_reads * 100 / $chip_total_reads" | bc -l)
    local peak_count=$(wc -l < $RESULTS_DIR/peaks/${OUTPUT_PREFIX}_peaks.narrowPeak)
    local frip_score=$(cat $RESULTS_DIR/peaks/frip_score.txt | cut -d' ' -f3)
    
    # 生成HTML报告
    cat > $report_file << EOF
<!DOCTYPE html>
<html>
<head>
    <title>ChIP-seq分析报告</title>
    <meta charset="UTF-8">
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        h1 { color: #1e3a8a; }
        h2 { color: #374151; border-bottom: 2px solid #e5e7eb; }
        .metric { background: #f3f4f6; padding: 10px; margin: 10px 0; border-radius: 5px; }
        .good { color: #059669; font-weight: bold; }
        .warning { color: #f59e0b; font-weight: bold; }
        .error { color: #dc2626; font-weight: bold; }
        table { border-collapse: collapse; width: 100%; }
        th, td { border: 1px solid #d1d5db; padding: 8px; text-align: left; }
        th { background-color: #f9fafb; }
    </style>
</head>
<body>
    <h1>ChIP-seq分析报告</h1>
    
    <h2>基本信息</h2>
    <div class="metric">
        <strong>样本：</strong> H3K4me3 ChIP-seq<br>
        <strong>分析时间：</strong> $(date)<br>
        <strong>分析流程：</strong> 质量控制 → 序列比对 → Peak calling → 注释分析
    </div>
    
    <h2>数据质量统计</h2>
    <table>
        <tr><th>指标</th><th>数值</th><th>评价</th></tr>
        <tr>
            <td>总reads数</td>
            <td>$chip_total_reads</td>
            <td>$([ $chip_total_reads -gt 10000000 ] && echo '<span class="good">良好</span>' || echo '<span class="warning">偏低</span>')</td>
        </tr>
        <tr>
            <td>比对reads数</td>
            <td>$chip_mapped_reads</td>
            <td>-</td>
        </tr>
        <tr>
            <td>比对率</td>
            <td>${mapping_rate}%</td>
            <td>$([ ${mapping_rate%.*} -gt 80 ] && echo '<span class="good">良好</span>' || echo '<span class="warning">偏低</span>')</td>
        </tr>
    </table>
    
    <h2>Peak Calling结果</h2>
    <table>
        <tr><th>指标</th><th>数值</th><th>评价</th></tr>
        <tr>
            <td>Peak数量</td>
            <td>$peak_count</td>
            <td>$([ $peak_count -gt 10000 ] && [ $peak_count -lt 50000 ] && echo '<span class="good">正常</span>' || echo '<span class="warning">异常</span>')</td>
        </tr>
        <tr>
            <td>FRiP Score</td>
            <td>$frip_score</td>
            <td>$([ $(echo "$frip_score > 0.05" | bc -l) -eq 1 ] && echo '<span class="good">良好</span>' || echo '<span class="warning">偏低</span>')</td>
        </tr>
    </table>
    
    <h2>输出文件</h2>
    <ul>
        <li><strong>Peak文件：</strong> results/peaks/${OUTPUT_PREFIX}_peaks.narrowPeak</li>
        <li><strong>可视化文件：</strong> results/${OUTPUT_PREFIX}_ChIP.bw</li>
        <li><strong>注释结果：</strong> results/annotation/peak_annotation.csv</li>
        <li><strong>图表文件：</strong> figures/目录下的PNG和PDF文件</li>
    </ul>
    
    <h2>分析建议</h2>
    <div class="metric">
        $([ $peak_count -gt 10000 ] && [ $peak_count -lt 50000 ] && echo "✓ Peak数量正常，适合进一步分析" || echo "⚠ Peak数量异常，建议检查实验条件和参数设置")<br>
        $([ $(echo "$frip_score > 0.05" | bc -l) -eq 1 ] && echo "✓ FRiP score良好，ChIP实验成功" || echo "⚠ FRiP score偏低，可能存在实验问题")<br>
        $([ ${mapping_rate%.*} -gt 80 ] && echo "✓ 比对率良好，数据质量可靠" || echo "⚠ 比对率偏低，建议检查数据质量")
    </div>
    
    <h2>下一步分析</h2>
    <ul>
        <li>在IGV中可视化重要基因区域的ChIP-seq信号</li>
        <li>进行差异结合分析（如有多个条件）</li>
        <li>结合RNA-seq数据分析转录调控关系</li>
        <li>进行motif分析识别转录因子结合基序</li>
    </ul>
    
    <hr>
    <p><small>报告生成时间：$(date) | 联系方式：wangys@hunau.edu.cn</small></p>
</body>
</html>
EOF
    
    log "✓ 分析报告已生成：$report_file"
}

# 主流程
main() {
    log "开始ChIP-seq完整分析流程"
    
    # 检查输入文件
    check_files
    
    # 执行分析步骤
    quality_control
    alignment
    peak_calling
    visualization
    peak_annotation
    generate_report
    
    log "=== ChIP-seq分析流程完成 ==="
    log "结束时间：$(date)"
    log "分析报告：$WORK_DIR/analysis_report.html"
    log "日志文件：$LOG_FILE"
    
    echo ""
    echo "分析完成！主要结果："
    echo "- Peak文件：$RESULTS_DIR/peaks/${OUTPUT_PREFIX}_peaks.narrowPeak"
    echo "- 可视化文件：$RESULTS_DIR/${OUTPUT_PREFIX}_ChIP.bw"
    echo "- 分析报告：$WORK_DIR/analysis_report.html"
    echo ""
    echo "建议下一步："
    echo "1. 在浏览器中查看分析报告"
    echo "2. 使用IGV查看重要基因区域"
    echo "3. 分析功能富集结果"
}

# 运行主流程
main "$@"