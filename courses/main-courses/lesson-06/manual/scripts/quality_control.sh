#!/bin/bash
# ChIP-seq质量控制脚本
# 作者：王运生
# 日期：2025-01-21
# 用法：bash quality_control.sh

set -e  # 遇到错误立即退出
set -u  # 使用未定义变量时报错

# 加载配置
source ~/ngs-analysis/lesson-06/scripts/config.sh

echo "=== ChIP-seq质量控制分析 ==="
echo "工作目录：$WORK_DIR"

# 创建输出目录
QC_DIR="$RESULTS_DIR/quality_control"
mkdir -p $QC_DIR/{fastqc,alignment_stats,peak_qc}

# 日志函数
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# 1. FastQC质量评估
fastqc_analysis() {
    log "=== 步骤1：FastQC质量评估 ==="
    
    local fastqc_dir="$QC_DIR/fastqc"
    
    # 对原始数据进行FastQC分析
    log "分析ChIP样本..."
    fastqc $DATA_DIR/H3K4me3_ChIP.fastq.gz -o $fastqc_dir --threads $THREADS
    
    log "分析Input对照..."
    fastqc $DATA_DIR/Input_control.fastq.gz -o $fastqc_dir --threads $THREADS
    
    # 生成MultiQC报告（如果安装了MultiQC）
    if command -v multiqc &> /dev/null; then
        log "生成MultiQC报告..."
        multiqc $fastqc_dir -o $fastqc_dir --title "ChIP-seq Quality Control"
    else
        log "MultiQC未安装，跳过综合报告生成"
    fi
    
    log "✓ FastQC分析完成"
}

# 2. 比对质量统计
alignment_stats() {
    log "=== 步骤2：比对质量统计 ==="
    
    local align_dir="$RESULTS_DIR/alignment"
    local stats_dir="$QC_DIR/alignment_stats"
    
    # 检查BAM文件是否存在
    if [ ! -f "$align_dir/H3K4me3_ChIP.bam" ]; then
        log "BAM文件不存在，跳过比对质量统计"
        return
    fi
    
    # 生成详细的比对统计
    log "生成比对统计..."
    
    # SAMtools flagstat
    samtools flagstat $align_dir/H3K4me3_ChIP.bam > $stats_dir/H3K4me3_ChIP.flagstat
    samtools flagstat $align_dir/Input_control.bam > $stats_dir/Input_control.flagstat
    
    # SAMtools stats
    samtools stats $align_dir/H3K4me3_ChIP.bam > $stats_dir/H3K4me3_ChIP.stats
    samtools stats $align_dir/Input_control.bam > $stats_dir/Input_control.stats
    
    # 计算比对率
    chip_total=$(samtools view -c $align_dir/H3K4me3_ChIP.bam)
    chip_mapped=$(samtools view -c -F 4 $align_dir/H3K4me3_ChIP.bam)
    chip_rate=$(echo "scale=2; $chip_mapped * 100 / $chip_total" | bc -l)
    
    input_total=$(samtools view -c $align_dir/Input_control.bam)
    input_mapped=$(samtools view -c -F 4 $align_dir/Input_control.bam)
    input_rate=$(echo "scale=2; $input_mapped * 100 / $input_total" | bc -l)
    
    # 生成比对质量报告
    cat > $stats_dir/alignment_summary.txt << EOF
ChIP-seq比对质量报告
==================

ChIP样本 (H3K4me3):
  总reads数: $chip_total
  比对reads数: $chip_mapped
  比对率: ${chip_rate}%

Input对照:
  总reads数: $input_total
  比对reads数: $input_mapped
  比对率: ${input_rate}%

质量评估:
$([ ${chip_rate%.*} -gt 80 ] && echo "✓ ChIP样本比对率良好" || echo "⚠ ChIP样本比对率偏低")
$([ ${input_rate%.*} -gt 80 ] && echo "✓ Input对照比对率良好" || echo "⚠ Input对照比对率偏低")
EOF
    
    log "✓ 比对质量统计完成"
}

# 3. Peak质量评估
peak_quality() {
    log "=== 步骤3：Peak质量评估 ==="
    
    local peaks_dir="$RESULTS_DIR/peaks"
    local peak_qc_dir="$QC_DIR/peak_qc"
    
    # 检查Peak文件是否存在
    if [ ! -f "$peaks_dir/${OUTPUT_PREFIX}_peaks.narrowPeak" ]; then
        log "Peak文件不存在，跳过Peak质量评估"
        return
    fi
    
    local peak_file="$peaks_dir/${OUTPUT_PREFIX}_peaks.narrowPeak"
    local bam_file="$RESULTS_DIR/alignment/H3K4me3_ChIP.bam"
    
    # Peak基本统计
    log "计算Peak基本统计..."
    local peak_count=$(wc -l < $peak_file)
    local total_peak_length=$(awk '{sum += $3-$2} END {print sum}' $peak_file)
    local avg_peak_length=$(echo "scale=0; $total_peak_length / $peak_count" | bc -l)
    
    # 计算FRiP score
    log "计算FRiP score..."
    local total_reads=$(samtools view -c $bam_file)
    local peak_reads=$(bedtools intersect -a $bam_file -b $peak_file -u | samtools view -c)
    local frip=$(echo "scale=4; $peak_reads / $total_reads" | bc -l)
    
    # Peak质量分布
    log "分析Peak质量分布..."
    awk '{print $5}' $peak_file | sort -nr > $peak_qc_dir/peak_scores.txt
    awk '{print $3-$2}' $peak_file | sort -n > $peak_qc_dir/peak_lengths.txt
    
    # 计算质量分位数
    local q25_score=$(awk 'NR==int(0.25*NR_total)+1' NR_total=$(wc -l < $peak_qc_dir/peak_scores.txt) $peak_qc_dir/peak_scores.txt)
    local q50_score=$(awk 'NR==int(0.50*NR_total)+1' NR_total=$(wc -l < $peak_qc_dir/peak_scores.txt) $peak_qc_dir/peak_scores.txt)
    local q75_score=$(awk 'NR==int(0.75*NR_total)+1' NR_total=$(wc -l < $peak_qc_dir/peak_scores.txt) $peak_qc_dir/peak_scores.txt)
    
    # 生成Peak质量报告
    cat > $peak_qc_dir/peak_quality_report.txt << EOF
ChIP-seq Peak质量报告
===================

基本统计:
  Peak数量: $peak_count
  总Peak长度: $total_peak_length bp
  平均Peak长度: $avg_peak_length bp
  FRiP Score: $frip

Peak得分分布:
  25%分位数: $q25_score
  50%分位数(中位数): $q50_score
  75%分位数: $q75_score

质量评估:
$([ $peak_count -gt 1000 ] && [ $peak_count -lt 100000 ] && echo "✓ Peak数量在合理范围内" || echo "⚠ Peak数量异常")
$([ $(echo "$frip > 0.01" | bc -l) -eq 1 ] && echo "✓ FRiP score达到最低标准" || echo "⚠ FRiP score过低")
$([ $(echo "$frip > 0.05" | bc -l) -eq 1 ] && echo "✓ FRiP score良好" || echo "⚠ FRiP score偏低")
$([ $avg_peak_length -gt 100 ] && [ $avg_peak_length -lt 2000 ] && echo "✓ 平均Peak长度合理" || echo "⚠ 平均Peak长度异常")

建议:
$([ $peak_count -lt 1000 ] && echo "- Peak数量过少，建议检查实验条件或降低阈值")
$([ $peak_count -gt 100000 ] && echo "- Peak数量过多，建议提高阈值或检查背景噪音")
$([ $(echo "$frip < 0.01" | bc -l) -eq 1 ] && echo "- FRiP score过低，建议检查ChIP实验效率")
EOF
    
    log "✓ Peak质量评估完成"
}

# 4. 生成质量控制图表
generate_qc_plots() {
    log "=== 步骤4：生成质量控制图表 ==="
    
    # 使用R生成QC图表
    R --slave --no-restore --no-save << 'EOF'
# 设置工作目录
setwd(Sys.getenv("WORK_DIR", "~/ngs-analysis/lesson-06"))

# 创建图表目录
dir.create("figures/qc", recursive = TRUE, showWarnings = FALSE)

# 读取Peak质量数据
peak_qc_dir <- "results/quality_control/peak_qc"

if (file.exists(file.path(peak_qc_dir, "peak_scores.txt")) && 
    file.exists(file.path(peak_qc_dir, "peak_lengths.txt"))) {
    
    # 读取数据
    peak_scores <- read.table(file.path(peak_qc_dir, "peak_scores.txt"), header = FALSE)$V1
    peak_lengths <- read.table(file.path(peak_qc_dir, "peak_lengths.txt"), header = FALSE)$V1
    
    # Peak得分分布图
    pdf("figures/qc/peak_score_distribution.pdf", width = 8, height = 6)
    hist(peak_scores, breaks = 50, main = "Peak Score Distribution", 
         xlab = "Peak Score", ylab = "Frequency", col = "lightblue")
    abline(v = median(peak_scores), col = "red", lwd = 2, lty = 2)
    legend("topright", legend = paste("Median =", round(median(peak_scores), 2)), 
           col = "red", lty = 2, lwd = 2)
    dev.off()
    
    # Peak长度分布图
    pdf("figures/qc/peak_length_distribution.pdf", width = 8, height = 6)
    hist(peak_lengths, breaks = 50, main = "Peak Length Distribution", 
         xlab = "Peak Length (bp)", ylab = "Frequency", col = "lightgreen")
    abline(v = median(peak_lengths), col = "red", lwd = 2, lty = 2)
    legend("topright", legend = paste("Median =", round(median(peak_lengths)), "bp"), 
           col = "red", lty = 2, lwd = 2)
    dev.off()
    
    # Peak得分vs长度散点图
    pdf("figures/qc/peak_score_vs_length.pdf", width = 8, height = 6)
    plot(peak_lengths, peak_scores, pch = 16, alpha = 0.6, 
         main = "Peak Score vs Length", 
         xlab = "Peak Length (bp)", ylab = "Peak Score")
    # 添加趋势线
    fit <- lm(peak_scores ~ peak_lengths)
    abline(fit, col = "red", lwd = 2)
    dev.off()
    
    cat("QC图表已生成到 figures/qc/ 目录\n")
} else {
    cat("Peak质量数据不存在，跳过图表生成\n")
}
EOF
    
    log "✓ 质量控制图表生成完成"
}

# 5. 生成综合质量报告
generate_comprehensive_report() {
    log "=== 步骤5：生成综合质量报告 ==="
    
    local report_file="$QC_DIR/comprehensive_qc_report.html"
    
    # 收集统计信息
    local chip_total_reads=$(samtools view -c $RESULTS_DIR/alignment/H3K4me3_ChIP.bam 2>/dev/null || echo "N/A")
    local chip_mapped_reads=$(samtools view -c -F 4 $RESULTS_DIR/alignment/H3K4me3_ChIP.bam 2>/dev/null || echo "N/A")
    
    if [ "$chip_total_reads" != "N/A" ] && [ "$chip_mapped_reads" != "N/A" ]; then
        local mapping_rate=$(echo "scale=2; $chip_mapped_reads * 100 / $chip_total_reads" | bc -l)
    else
        local mapping_rate="N/A"
    fi
    
    local peak_count=$(wc -l < $RESULTS_DIR/peaks/${OUTPUT_PREFIX}_peaks.narrowPeak 2>/dev/null || echo "N/A")
    local frip_score=$(grep "FRiP Score:" $QC_DIR/peak_qc/peak_quality_report.txt 2>/dev/null | cut -d' ' -f3 || echo "N/A")
    
    # 生成HTML报告
    cat > $report_file << EOF
<!DOCTYPE html>
<html>
<head>
    <title>ChIP-seq质量控制报告</title>
    <meta charset="UTF-8">
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; }
        h1 { color: #1e3a8a; border-bottom: 3px solid #1e3a8a; }
        h2 { color: #374151; border-bottom: 2px solid #e5e7eb; }
        h3 { color: #4b5563; }
        .metric { background: #f3f4f6; padding: 15px; margin: 10px 0; border-radius: 8px; border-left: 4px solid #3b82f6; }
        .good { color: #059669; font-weight: bold; }
        .warning { color: #f59e0b; font-weight: bold; }
        .error { color: #dc2626; font-weight: bold; }
        table { border-collapse: collapse; width: 100%; margin: 20px 0; }
        th, td { border: 1px solid #d1d5db; padding: 12px; text-align: left; }
        th { background-color: #f9fafb; font-weight: bold; }
        .summary-box { background: #eff6ff; padding: 20px; border-radius: 8px; margin: 20px 0; }
        .file-list { background: #f0fdf4; padding: 15px; border-radius: 8px; }
        ul { padding-left: 20px; }
        li { margin: 5px 0; }
    </style>
</head>
<body>
    <h1>ChIP-seq质量控制综合报告</h1>
    
    <div class="summary-box">
        <h2>分析概要</h2>
        <p><strong>样本：</strong> H3K4me3 ChIP-seq</p>
        <p><strong>分析时间：</strong> $(date)</p>
        <p><strong>分析流程：</strong> 原始数据质量评估 → 比对质量统计 → Peak质量评估</p>
    </div>
    
    <h2>1. 原始数据质量</h2>
    <div class="metric">
        <h3>FastQC结果</h3>
        <p>详细的FastQC报告请查看：</p>
        <ul>
            <li>ChIP样本：<code>results/quality_control/fastqc/H3K4me3_ChIP_fastqc.html</code></li>
            <li>Input对照：<code>results/quality_control/fastqc/Input_control_fastqc.html</code></li>
        </ul>
        <p><strong>关键指标检查：</strong></p>
        <ul>
            <li>Per base sequence quality: 应该大部分在绿色区域</li>
            <li>Per sequence quality scores: 峰值应该在高质量区域</li>
            <li>Sequence duplication levels: 适度重复是正常的</li>
            <li>Adapter content: 应该很少或没有adapter污染</li>
        </ul>
    </div>
    
    <h2>2. 比对质量统计</h2>
    <table>
        <tr><th>样本</th><th>总reads数</th><th>比对reads数</th><th>比对率</th><th>评价</th></tr>
        <tr>
            <td>H3K4me3 ChIP</td>
            <td>$chip_total_reads</td>
            <td>$chip_mapped_reads</td>
            <td>${mapping_rate}%</td>
            <td>$([ "$mapping_rate" != "N/A" ] && [ ${mapping_rate%.*} -gt 80 ] && echo '<span class="good">良好</span>' || echo '<span class="warning">需关注</span>')</td>
        </tr>
    </table>
    
    <div class="metric">
        <h3>比对质量标准</h3>
        <ul>
            <li><span class="good">优秀：</span> 比对率 > 90%</li>
            <li><span class="good">良好：</span> 比对率 80-90%</li>
            <li><span class="warning">可接受：</span> 比对率 70-80%</li>
            <li><span class="error">需改进：</span> 比对率 < 70%</li>
        </ul>
    </div>
    
    <h2>3. Peak质量评估</h2>
    <table>
        <tr><th>指标</th><th>数值</th><th>标准</th><th>评价</th></tr>
        <tr>
            <td>Peak数量</td>
            <td>$peak_count</td>
            <td>H3K4me3: 10,000-50,000</td>
            <td>$([ "$peak_count" != "N/A" ] && [ $peak_count -gt 10000 ] && [ $peak_count -lt 50000 ] && echo '<span class="good">正常</span>' || echo '<span class="warning">异常</span>')</td>
        </tr>
        <tr>
            <td>FRiP Score</td>
            <td>$frip_score</td>
            <td>> 0.05 (理想 > 0.1)</td>
            <td>$([ "$frip_score" != "N/A" ] && [ $(echo "$frip_score > 0.05" | bc -l 2>/dev/null || echo 0) -eq 1 ] && echo '<span class="good">良好</span>' || echo '<span class="warning">偏低</span>')</td>
        </tr>
    </table>
    
    <div class="metric">
        <h3>Peak质量标准</h3>
        <ul>
            <li><strong>FRiP Score (Fraction of Reads in Peaks)：</strong></li>
            <ul>
                <li><span class="good">优秀：</span> > 0.1 (10%)</li>
                <li><span class="good">良好：</span> 0.05-0.1 (5-10%)</li>
                <li><span class="warning">可接受：</span> 0.01-0.05 (1-5%)</li>
                <li><span class="error">需改进：</span> < 0.01 (1%)</li>
            </ul>
            <li><strong>Peak数量：</strong> 因蛋白质类型而异，H3K4me3通常10,000-50,000个</li>
        </ul>
    </div>
    
    <h2>4. 质量控制图表</h2>
    <div class="file-list">
        <h3>生成的图表文件：</h3>
        <ul>
            <li><code>figures/qc/peak_score_distribution.pdf</code> - Peak得分分布</li>
            <li><code>figures/qc/peak_length_distribution.pdf</code> - Peak长度分布</li>
            <li><code>figures/qc/peak_score_vs_length.pdf</code> - Peak得分与长度关系</li>
        </ul>
    </div>
    
    <h2>5. 问题诊断和建议</h2>
    <div class="metric">
        <h3>常见问题及解决方案：</h3>
        <ul>
            <li><strong>比对率低 (<70%)：</strong>
                <ul>
                    <li>检查参考基因组版本是否正确</li>
                    <li>检查测序数据质量</li>
                    <li>考虑更换比对软件或参数</li>
                </ul>
            </li>
            <li><strong>FRiP score低 (<1%)：</strong>
                <ul>
                    <li>检查ChIP实验效率</li>
                    <li>检查抗体特异性</li>
                    <li>调整MACS2参数（降低q-value阈值）</li>
                </ul>
            </li>
            <li><strong>Peak数量异常：</strong>
                <ul>
                    <li>过多：提高q-value阈值，检查背景噪音</li>
                    <li>过少：降低q-value阈值，检查实验条件</li>
                </ul>
            </li>
        </ul>
    </div>
    
    <h2>6. 下一步分析建议</h2>
    <div class="summary-box">
        <p>基于质量控制结果，建议的下一步分析：</p>
        <ul>
            <li>如果所有质量指标都良好，可以进行Peak注释和功能分析</li>
            <li>如果存在质量问题，建议先优化实验条件或分析参数</li>
            <li>可以考虑与其他ChIP-seq数据进行比较分析</li>
            <li>结合RNA-seq数据分析转录调控关系</li>
        </ul>
    </div>
    
    <h2>7. 输出文件清单</h2>
    <div class="file-list">
        <h3>质量控制相关文件：</h3>
        <ul>
            <li><code>results/quality_control/fastqc/</code> - FastQC报告</li>
            <li><code>results/quality_control/alignment_stats/</code> - 比对统计</li>
            <li><code>results/quality_control/peak_qc/</code> - Peak质量评估</li>
            <li><code>figures/qc/</code> - 质量控制图表</li>
        </ul>
    </div>
    
    <hr>
    <p><small>报告生成时间：$(date) | 联系方式：wangys@hunau.edu.cn</small></p>
</body>
</html>
EOF
    
    log "✓ 综合质量报告已生成：$report_file"
}

# 主函数
main() {
    log "开始ChIP-seq质量控制分析"
    
    # 执行质量控制步骤
    fastqc_analysis
    alignment_stats
    peak_quality
    generate_qc_plots
    generate_comprehensive_report
    
    log "=== ChIP-seq质量控制分析完成 ==="
    log "综合报告：$QC_DIR/comprehensive_qc_report.html"
    
    echo ""
    echo "质量控制分析完成！"
    echo "主要输出："
    echo "- FastQC报告：$QC_DIR/fastqc/"
    echo "- 比对统计：$QC_DIR/alignment_stats/"
    echo "- Peak质量：$QC_DIR/peak_qc/"
    echo "- 综合报告：$QC_DIR/comprehensive_qc_report.html"
    echo ""
    echo "建议："
    echo "1. 在浏览器中查看综合质量报告"
    echo "2. 检查FastQC报告中的关键指标"
    echo "3. 根据质量评估结果决定是否需要优化参数"
}

# 运行主函数
main "$@"