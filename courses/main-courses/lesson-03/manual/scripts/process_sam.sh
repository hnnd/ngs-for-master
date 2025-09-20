#!/bin/bash
# SAM/BAM文件处理脚本
# 作者：王运生
# 日期：2025-01-20
# 用法：bash process_sam.sh [input_bam] [output_prefix]

set -e  # 遇到错误立即退出
set -u  # 使用未定义变量时报错

echo "=========================================="
echo "SAM/BAM文件处理和分析"
echo "=========================================="

# 参数设置
INPUT_BAM=${1:-"results/bwa_default_sorted.bam"}
OUTPUT_PREFIX=${2:-"processed"}

# 检查输入文件
if [ ! -f "$INPUT_BAM" ]; then
    echo "错误：输入BAM文件不存在: $INPUT_BAM"
    echo "用法: $0 [input_bam] [output_prefix]"
    exit 1
fi

echo "输入BAM文件: $INPUT_BAM"
echo "输出前缀: $OUTPUT_PREFIX"
echo "文件大小: $(du -h "$INPUT_BAM" | cut -f1)"

# 创建输出目录
mkdir -p results logs

# 记录开始时间
START_TIME=$(date +%s)
echo "开始时间: $(date)"

# 1. 基本文件信息
echo ""
echo "1. 分析BAM文件基本信息..."
echo "----------------------------------------"

echo "BAM文件头部信息:"
samtools view -H "$INPUT_BAM" | head -10

echo ""
echo "BAM文件统计信息:"
samtools idxstats "$INPUT_BAM" > results/${OUTPUT_PREFIX}_idxstats.txt
cat results/${OUTPUT_PREFIX}_idxstats.txt

# 2. 比对质量过滤
echo ""
echo "2. 进行比对质量过滤..."
echo "----------------------------------------"

# 过滤MAPQ < 20的比对
echo "过滤MAPQ < 20的比对..."
samtools view -q 20 -b "$INPUT_BAM" > results/${OUTPUT_PREFIX}_mapq20.bam
samtools index results/${OUTPUT_PREFIX}_mapq20.bam

# 过滤MAPQ < 30的比对
echo "过滤MAPQ < 30的比对..."
samtools view -q 30 -b "$INPUT_BAM" > results/${OUTPUT_PREFIX}_mapq30.bam
samtools index results/${OUTPUT_PREFIX}_mapq30.bam

# 统计过滤结果
echo "过滤统计:"
echo "原始BAM: $(samtools view -c "$INPUT_BAM") reads"
echo "MAPQ≥20: $(samtools view -c results/${OUTPUT_PREFIX}_mapq20.bam) reads"
echo "MAPQ≥30: $(samtools view -c results/${OUTPUT_PREFIX}_mapq30.bam) reads"

# 3. 配对状态过滤
echo ""
echo "3. 进行配对状态过滤..."
echo "----------------------------------------"

# 只保留正确配对的reads
echo "提取正确配对的reads..."
samtools view -f 2 -b "$INPUT_BAM" > results/${OUTPUT_PREFIX}_proper_pairs.bam
samtools index results/${OUTPUT_PREFIX}_proper_pairs.bam

# 提取单端比对的reads
echo "提取单端比对的reads..."
samtools view -F 2 -f 1 -b "$INPUT_BAM" > results/${OUTPUT_PREFIX}_singleton.bam
samtools index results/${OUTPUT_PREFIX}_singleton.bam

# 提取未比对的reads
echo "提取未比对的reads..."
samtools view -f 4 -b "$INPUT_BAM" > results/${OUTPUT_PREFIX}_unmapped.bam

# 统计配对结果
echo "配对统计:"
echo "正确配对: $(samtools view -c results/${OUTPUT_PREFIX}_proper_pairs.bam) reads"
echo "单端比对: $(samtools view -c results/${OUTPUT_PREFIX}_singleton.bam) reads"
echo "未比对: $(samtools view -c results/${OUTPUT_PREFIX}_unmapped.bam) reads"

# 4. 重复序列标记
echo ""
echo "4. 标记重复序列..."
echo "----------------------------------------"

# 使用samtools markdup标记重复
echo "使用samtools markdup标记重复序列..."
samtools markdup "$INPUT_BAM" results/${OUTPUT_PREFIX}_markdup.bam 2> logs/${OUTPUT_PREFIX}_markdup.log
samtools index results/${OUTPUT_PREFIX}_markdup.bam

# 统计重复率
echo "重复序列统计:"
grep "DUPLICATE" logs/${OUTPUT_PREFIX}_markdup.log || echo "未找到重复统计信息"

# 5. 插入片段大小分析
echo ""
echo "5. 分析插入片段大小..."
echo "----------------------------------------"

# 提取插入片段大小信息
echo "提取插入片段大小..."
samtools view -f 2 "$INPUT_BAM" | \
awk '$9 > 0 {print $9}' | \
sort -n > results/${OUTPUT_PREFIX}_insert_sizes.txt

# 计算插入片段统计
echo "计算插入片段统计..."
awk '{
    sum += $1; 
    sumsq += $1*$1; 
    count++; 
    sizes[count] = $1
} END {
    mean = sum/count
    variance = (sumsq - sum*sum/count)/(count-1)
    stddev = sqrt(variance)
    
    # 计算中位数
    n = asort(sizes)
    if (n % 2 == 1) {
        median = sizes[(n+1)/2]
    } else {
        median = (sizes[n/2] + sizes[n/2+1])/2
    }
    
    printf "插入片段统计:\n"
    printf "  总数: %d\n", count
    printf "  平均值: %.2f\n", mean
    printf "  中位数: %.2f\n", median
    printf "  标准差: %.2f\n", stddev
    printf "  最小值: %d\n", sizes[1]
    printf "  最大值: %d\n", sizes[n]
}' results/${OUTPUT_PREFIX}_insert_sizes.txt > results/${OUTPUT_PREFIX}_insert_stats.txt

cat results/${OUTPUT_PREFIX}_insert_stats.txt

# 6. 比对位置分析
echo ""
echo "6. 分析比对位置分布..."
echo "----------------------------------------"

# 统计每个染色体的比对数量
echo "统计染色体比对分布..."
samtools view "$INPUT_BAM" | \
cut -f3 | \
sort | \
uniq -c | \
sort -nr > results/${OUTPUT_PREFIX}_chr_distribution.txt

echo "染色体比对分布:"
head -10 results/${OUTPUT_PREFIX}_chr_distribution.txt

# 7. 软剪切分析
echo ""
echo "7. 分析软剪切情况..."
echo "----------------------------------------"

# 提取包含软剪切的reads
echo "统计软剪切reads..."
samtools view "$INPUT_BAM" | \
awk '$6 ~ /S/ {print $6}' | \
sort | \
uniq -c | \
sort -nr > results/${OUTPUT_PREFIX}_soft_clipping.txt

echo "软剪切统计:"
head -10 results/${OUTPUT_PREFIX}_soft_clipping.txt

# 计算软剪切比例
total_reads=$(samtools view -c "$INPUT_BAM")
soft_clipped=$(samtools view "$INPUT_BAM" | awk '$6 ~ /S/ {count++} END {print count+0}')
if [ "$total_reads" -gt 0 ]; then
    soft_clip_rate=$(echo "scale=4; $soft_clipped * 100 / $total_reads" | bc -l)
    echo "软剪切比例: ${soft_clip_rate}%"
fi

# 8. 错配分析
echo ""
echo "8. 分析错配情况..."
echo "----------------------------------------"

# 提取NM标签（错配数）
echo "统计错配分布..."
samtools view "$INPUT_BAM" | \
grep -o "NM:i:[0-9]*" | \
cut -d: -f3 | \
sort -n | \
uniq -c > results/${OUTPUT_PREFIX}_mismatches.txt

echo "错配分布:"
head -10 results/${OUTPUT_PREFIX}_mismatches.txt

# 9. 生成综合报告
echo ""
echo "9. 生成综合处理报告..."
echo "----------------------------------------"

cat > results/${OUTPUT_PREFIX}_processing_report.txt << EOF
SAM/BAM文件处理报告
生成时间: $(date)
输入文件: $INPUT_BAM

文件基本信息:
- 文件大小: $(du -h "$INPUT_BAM" | cut -f1)
- 总reads数: $(samtools view -c "$INPUT_BAM")
- 索引状态: $([ -f "${INPUT_BAM}.bai" ] && echo "已建立" || echo "未建立")

质量过滤结果:
- MAPQ≥20: $(samtools view -c results/${OUTPUT_PREFIX}_mapq20.bam) reads
- MAPQ≥30: $(samtools view -c results/${OUTPUT_PREFIX}_mapq30.bam) reads

配对状态统计:
- 正确配对: $(samtools view -c results/${OUTPUT_PREFIX}_proper_pairs.bam) reads
- 单端比对: $(samtools view -c results/${OUTPUT_PREFIX}_singleton.bam) reads  
- 未比对: $(samtools view -c results/${OUTPUT_PREFIX}_unmapped.bam) reads

插入片段信息:
EOF

cat results/${OUTPUT_PREFIX}_insert_stats.txt >> results/${OUTPUT_PREFIX}_processing_report.txt

cat >> results/${OUTPUT_PREFIX}_processing_report.txt << EOF

软剪切统计:
- 包含软剪切的reads: $soft_clipped
- 软剪切比例: ${soft_clip_rate:-"未计算"}%

生成的文件:
- 质量过滤: ${OUTPUT_PREFIX}_mapq20.bam, ${OUTPUT_PREFIX}_mapq30.bam
- 配对过滤: ${OUTPUT_PREFIX}_proper_pairs.bam, ${OUTPUT_PREFIX}_singleton.bam
- 重复标记: ${OUTPUT_PREFIX}_markdup.bam
- 统计文件: ${OUTPUT_PREFIX}_*_stats.txt, ${OUTPUT_PREFIX}_*.txt
EOF

# 10. 创建处理后文件的快速统计
echo ""
echo "10. 生成处理后文件统计..."
echo "----------------------------------------"

cat > results/${OUTPUT_PREFIX}_file_summary.txt << EOF
处理后文件统计摘要
生成时间: $(date)

文件名                           大小      reads数    描述
================================================================
EOF

for bam_file in results/${OUTPUT_PREFIX}_*.bam; do
    if [ -f "$bam_file" ]; then
        filename=$(basename "$bam_file")
        size=$(du -h "$bam_file" | cut -f1)
        reads=$(samtools view -c "$bam_file" 2>/dev/null || echo "0")
        
        case "$filename" in
            *mapq20*) desc="MAPQ≥20过滤" ;;
            *mapq30*) desc="MAPQ≥30过滤" ;;
            *proper_pairs*) desc="正确配对" ;;
            *singleton*) desc="单端比对" ;;
            *unmapped*) desc="未比对reads" ;;
            *markdup*) desc="重复标记" ;;
            *) desc="其他处理" ;;
        esac
        
        printf "%-30s %8s %10s    %s\n" "$filename" "$size" "$reads" "$desc" >> results/${OUTPUT_PREFIX}_file_summary.txt
    fi
done

cat results/${OUTPUT_PREFIX}_file_summary.txt

# 计算总用时
END_TIME=$(date +%s)
TOTAL_TIME=$((END_TIME - START_TIME))

echo ""
echo "=========================================="
echo "SAM/BAM文件处理完成！"
echo "=========================================="
echo "总用时: ${TOTAL_TIME}秒"
echo ""
echo "生成的主要文件:"
echo "- 处理报告: results/${OUTPUT_PREFIX}_processing_report.txt"
echo "- 文件摘要: results/${OUTPUT_PREFIX}_file_summary.txt"
echo "- 质量过滤: results/${OUTPUT_PREFIX}_mapq*.bam"
echo "- 配对过滤: results/${OUTPUT_PREFIX}_proper_pairs.bam"
echo "- 重复标记: results/${OUTPUT_PREFIX}_markdup.bam"
echo "- 统计文件: results/${OUTPUT_PREFIX}_*_stats.txt"
echo ""
echo "可以使用以下命令查看结果:"
echo "  cat results/${OUTPUT_PREFIX}_processing_report.txt"
echo "  cat results/${OUTPUT_PREFIX}_file_summary.txt"
echo "  samtools flagstat results/${OUTPUT_PREFIX}_mapq30.bam"

# 记录到日志
cat >> logs/experiment.log << EOF

SAM/BAM处理完成时间: $(date)
输入文件: $INPUT_BAM
输出前缀: $OUTPUT_PREFIX
处理用时: ${TOTAL_TIME}秒
EOF

echo ""
echo "实验日志已更新: logs/experiment.log"