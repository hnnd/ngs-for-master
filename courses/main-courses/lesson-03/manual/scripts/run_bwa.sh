#!/bin/bash
# BWA序列比对脚本
# 作者：王运生
# 日期：2025-01-20
# 用法：bash run_bwa.sh [threads] [output_prefix]

set -e  # 遇到错误立即退出
set -u  # 使用未定义变量时报错

echo "=========================================="
echo "BWA序列比对分析"
echo "=========================================="

# 参数设置
THREADS=${1:-4}
OUTPUT_PREFIX=${2:-"bwa"}
REFERENCE="reference/chr22_fragment.fa"
READ1="data/sample_R1.fastq"
READ2="data/sample_R2.fastq"

# 检查输入文件
echo "检查输入文件..."
for file in "$REFERENCE" "$READ1" "$READ2"; do
    if [ ! -f "$file" ]; then
        echo "错误：文件不存在: $file"
        exit 1
    fi
    echo "✓ $file ($(du -h "$file" | cut -f1))"
done

# 检查BWA索引
echo ""
echo "检查BWA索引..."
if [ ! -f "${REFERENCE}.bwt" ]; then
    echo "错误：BWA索引不存在，请先运行 build_index.sh"
    exit 1
fi
echo "✓ BWA索引文件存在"

# 创建输出目录
mkdir -p results logs

# 记录开始时间
START_TIME=$(date +%s)
echo ""
echo "开始BWA比对分析..."
echo "参考基因组: $REFERENCE"
echo "测序数据: $READ1, $READ2"
echo "线程数: $THREADS"
echo "输出前缀: $OUTPUT_PREFIX"
echo "开始时间: $(date)"

# 1. BWA-MEM默认参数比对
echo ""
echo "1. BWA-MEM默认参数比对..."
echo "----------------------------------------"

BWA_DEFAULT_START=$(date +%s)

bwa mem \
    -t $THREADS \
    -M \
    -R "@RG\tID:${OUTPUT_PREFIX}_default\tSM:sample\tPL:illumina\tLB:lib1\tPU:unit1" \
    "$REFERENCE" \
    "$READ1" \
    "$READ2" \
    2> logs/${OUTPUT_PREFIX}_default.log \
    > results/${OUTPUT_PREFIX}_default.sam

BWA_DEFAULT_END=$(date +%s)
BWA_DEFAULT_TIME=$((BWA_DEFAULT_END - BWA_DEFAULT_START))

echo "默认参数比对完成，用时: ${BWA_DEFAULT_TIME}秒"
echo "输出文件: results/${OUTPUT_PREFIX}_default.sam"

# 2. BWA-MEM优化参数比对
echo ""
echo "2. BWA-MEM优化参数比对..."
echo "----------------------------------------"

BWA_OPT_START=$(date +%s)

bwa mem \
    -t $THREADS \
    -k 19 \
    -w 100 \
    -d 100 \
    -r 1.5 \
    -A 1 -B 4 -O 6 -E 1 \
    -L 5 \
    -M \
    -R "@RG\tID:${OUTPUT_PREFIX}_optimized\tSM:sample\tPL:illumina\tLB:lib1\tPU:unit1" \
    "$REFERENCE" \
    "$READ1" \
    "$READ2" \
    2> logs/${OUTPUT_PREFIX}_optimized.log \
    > results/${OUTPUT_PREFIX}_optimized.sam

BWA_OPT_END=$(date +%s)
BWA_OPT_TIME=$((BWA_OPT_END - BWA_OPT_START))

echo "优化参数比对完成，用时: ${BWA_OPT_TIME}秒"
echo "输出文件: results/${OUTPUT_PREFIX}_optimized.sam"

# 3. BWA-MEM严格参数比对（用于高精度分析）
echo ""
echo "3. BWA-MEM严格参数比对..."
echo "----------------------------------------"

BWA_STRICT_START=$(date +%s)

bwa mem \
    -t $THREADS \
    -k 25 \
    -w 200 \
    -d 200 \
    -r 2.0 \
    -A 2 -B 8 -O 12 -E 2 \
    -L 10 \
    -M \
    -R "@RG\tID:${OUTPUT_PREFIX}_strict\tSM:sample\tPL:illumina\tLB:lib1\tPU:unit1" \
    "$REFERENCE" \
    "$READ1" \
    "$READ2" \
    2> logs/${OUTPUT_PREFIX}_strict.log \
    > results/${OUTPUT_PREFIX}_strict.sam

BWA_STRICT_END=$(date +%s)
BWA_STRICT_TIME=$((BWA_STRICT_END - BWA_STRICT_START))

echo "严格参数比对完成，用时: ${BWA_STRICT_TIME}秒"
echo "输出文件: results/${OUTPUT_PREFIX}_strict.sam"

# 4. 转换SAM为BAM并排序
echo ""
echo "4. 转换SAM为BAM格式并排序..."
echo "----------------------------------------"

for mode in default optimized strict; do
    echo "处理 ${mode} 模式结果..."
    
    SAM_FILE="results/${OUTPUT_PREFIX}_${mode}.sam"
    BAM_FILE="results/${OUTPUT_PREFIX}_${mode}_sorted.bam"
    
    # 转换为BAM并排序
    samtools view -bS "$SAM_FILE" | \
    samtools sort -@ $((THREADS-1)) -o "$BAM_FILE"
    
    # 建立索引
    samtools index "$BAM_FILE"
    
    # 删除SAM文件以节省空间（可选）
    # rm "$SAM_FILE"
    
    echo "  ✓ $BAM_FILE"
done

# 5. 生成比对统计
echo ""
echo "5. 生成比对统计..."
echo "----------------------------------------"

for mode in default optimized strict; do
    BAM_FILE="results/${OUTPUT_PREFIX}_${mode}_sorted.bam"
    STATS_FILE="results/${OUTPUT_PREFIX}_${mode}_stats.txt"
    
    echo "统计 ${mode} 模式结果..."
    
    # 基本统计
    samtools flagstat "$BAM_FILE" > "$STATS_FILE"
    
    # 详细统计
    samtools stats "$BAM_FILE" >> "$STATS_FILE"
    
    echo "  ✓ $STATS_FILE"
done

# 6. 提取关键统计信息
echo ""
echo "6. 提取关键统计信息..."
echo "----------------------------------------"

cat > results/${OUTPUT_PREFIX}_summary.txt << EOF
BWA比对结果摘要
生成时间: $(date)

比对参数:
- 线程数: $THREADS
- 参考基因组: $REFERENCE
- 测序数据: $READ1, $READ2

运行时间:
- 默认参数: ${BWA_DEFAULT_TIME}秒
- 优化参数: ${BWA_OPT_TIME}秒  
- 严格参数: ${BWA_STRICT_TIME}秒

比对统计:
EOF

for mode in default optimized strict; do
    STATS_FILE="results/${OUTPUT_PREFIX}_${mode}_stats.txt"
    
    echo "" >> results/${OUTPUT_PREFIX}_summary.txt
    echo "${mode} 模式:" >> results/${OUTPUT_PREFIX}_summary.txt
    
    # 提取关键统计信息
    total_reads=$(grep "in total" "$STATS_FILE" | cut -d' ' -f1)
    mapped_reads=$(grep "mapped (" "$STATS_FILE" | head -1 | cut -d' ' -f1)
    paired_reads=$(grep "properly paired" "$STATS_FILE" | cut -d' ' -f1)
    
    if [ -n "$total_reads" ] && [ -n "$mapped_reads" ]; then
        mapping_rate=$(echo "scale=2; $mapped_reads * 100 / $total_reads" | bc -l 2>/dev/null || echo "计算错误")
        echo "  总reads数: $total_reads" >> results/${OUTPUT_PREFIX}_summary.txt
        echo "  比对reads数: $mapped_reads" >> results/${OUTPUT_PREFIX}_summary.txt
        echo "  比对率: ${mapping_rate}%" >> results/${OUTPUT_PREFIX}_summary.txt
        
        if [ -n "$paired_reads" ]; then
            pairing_rate=$(echo "scale=2; $paired_reads * 100 / $total_reads" | bc -l 2>/dev/null || echo "计算错误")
            echo "  正确配对数: $paired_reads" >> results/${OUTPUT_PREFIX}_summary.txt
            echo "  配对率: ${pairing_rate}%" >> results/${OUTPUT_PREFIX}_summary.txt
        fi
    fi
done

# 7. 生成MAPQ分布统计
echo ""
echo "7. 生成MAPQ分布统计..."
echo "----------------------------------------"

for mode in default optimized strict; do
    BAM_FILE="results/${OUTPUT_PREFIX}_${mode}_sorted.bam"
    MAPQ_FILE="results/${OUTPUT_PREFIX}_${mode}_mapq.txt"
    
    echo "分析 ${mode} 模式MAPQ分布..."
    
    # 提取MAPQ值并统计分布
    samtools view "$BAM_FILE" | cut -f5 | sort -n | uniq -c > "$MAPQ_FILE"
    
    echo "  ✓ $MAPQ_FILE"
done

# 8. 计算覆盖度
echo ""
echo "8. 计算基因组覆盖度..."
echo "----------------------------------------"

for mode in default optimized strict; do
    BAM_FILE="results/${OUTPUT_PREFIX}_${mode}_sorted.bam"
    COV_FILE="results/${OUTPUT_PREFIX}_${mode}_coverage.txt"
    
    echo "计算 ${mode} 模式覆盖度..."
    
    # 计算每个位点的覆盖度
    samtools depth "$BAM_FILE" > "$COV_FILE"
    
    # 统计覆盖度信息
    awk '{
        sum += $3; 
        count++; 
        if($3 > 0) covered++
    } END {
        printf "平均覆盖度: %.2f\n", sum/count
        printf "覆盖率: %.2f%%\n", covered/count*100
        printf "总位点数: %d\n", count
        printf "覆盖位点数: %d\n", covered
    }' "$COV_FILE" > "results/${OUTPUT_PREFIX}_${mode}_coverage_summary.txt"
    
    echo "  ✓ $COV_FILE"
    echo "  ✓ results/${OUTPUT_PREFIX}_${mode}_coverage_summary.txt"
done

# 计算总用时
END_TIME=$(date +%s)
TOTAL_TIME=$((END_TIME - START_TIME))

echo ""
echo "=========================================="
echo "BWA比对分析完成！"
echo "=========================================="
echo "总用时: ${TOTAL_TIME}秒"
echo ""
echo "生成的文件:"
echo "- SAM文件: results/${OUTPUT_PREFIX}_*.sam"
echo "- BAM文件: results/${OUTPUT_PREFIX}_*_sorted.bam"
echo "- 统计文件: results/${OUTPUT_PREFIX}_*_stats.txt"
echo "- 摘要文件: results/${OUTPUT_PREFIX}_summary.txt"
echo "- MAPQ分布: results/${OUTPUT_PREFIX}_*_mapq.txt"
echo "- 覆盖度文件: results/${OUTPUT_PREFIX}_*_coverage.txt"
echo ""
echo "可以使用以下命令查看结果:"
echo "  cat results/${OUTPUT_PREFIX}_summary.txt"
echo "  samtools view results/${OUTPUT_PREFIX}_default_sorted.bam | head"
echo "  samtools tview results/${OUTPUT_PREFIX}_default_sorted.bam $REFERENCE"

# 记录到日志
cat >> logs/experiment.log << EOF

BWA比对完成时间: $(date)
总用时: ${TOTAL_TIME}秒
默认参数用时: ${BWA_DEFAULT_TIME}秒
优化参数用时: ${BWA_OPT_TIME}秒
严格参数用时: ${BWA_STRICT_TIME}秒
EOF

echo ""
echo "实验日志已更新: logs/experiment.log"