#!/bin/bash
# Bowtie2序列比对脚本
# 作者：王运生
# 日期：2025-01-20
# 用法：bash run_bowtie2.sh [threads] [output_prefix]

set -e  # 遇到错误立即退出
set -u  # 使用未定义变量时报错

echo "=========================================="
echo "Bowtie2序列比对分析"
echo "=========================================="

# 参数设置
THREADS=${1:-4}
OUTPUT_PREFIX=${2:-"bowtie2"}
REFERENCE_PREFIX="reference/chr22_fragment_bt2"
READ1="data/sample_R1.fastq"
READ2="data/sample_R2.fastq"

# 检查输入文件
echo "检查输入文件..."
for file in "$READ1" "$READ2"; do
    if [ ! -f "$file" ]; then
        echo "错误：文件不存在: $file"
        exit 1
    fi
    echo "✓ $file ($(du -h "$file" | cut -f1))"
done

# 检查Bowtie2索引
echo ""
echo "检查Bowtie2索引..."
if [ ! -f "${REFERENCE_PREFIX}.1.bt2" ]; then
    echo "错误：Bowtie2索引不存在，请先运行 build_index.sh"
    exit 1
fi
echo "✓ Bowtie2索引文件存在"

# 创建输出目录
mkdir -p results logs

# 记录开始时间
START_TIME=$(date +%s)
echo ""
echo "开始Bowtie2比对分析..."
echo "索引前缀: $REFERENCE_PREFIX"
echo "测序数据: $READ1, $READ2"
echo "线程数: $THREADS"
echo "输出前缀: $OUTPUT_PREFIX"
echo "开始时间: $(date)"

# 1. Bowtie2快速模式
echo ""
echo "1. Bowtie2快速模式比对..."
echo "----------------------------------------"

BT2_FAST_START=$(date +%s)

bowtie2 \
    --fast \
    -p $THREADS \
    --rg-id "${OUTPUT_PREFIX}_fast" \
    --rg "SM:sample" \
    --rg "PL:illumina" \
    --rg "LB:lib1" \
    -x "$REFERENCE_PREFIX" \
    -1 "$READ1" \
    -2 "$READ2" \
    -S results/${OUTPUT_PREFIX}_fast.sam \
    2> logs/${OUTPUT_PREFIX}_fast.log

BT2_FAST_END=$(date +%s)
BT2_FAST_TIME=$((BT2_FAST_END - BT2_FAST_START))

echo "快速模式比对完成，用时: ${BT2_FAST_TIME}秒"
echo "输出文件: results/${OUTPUT_PREFIX}_fast.sam"

# 2. Bowtie2敏感模式（默认）
echo ""
echo "2. Bowtie2敏感模式比对..."
echo "----------------------------------------"

BT2_SENS_START=$(date +%s)

bowtie2 \
    --sensitive \
    -p $THREADS \
    --rg-id "${OUTPUT_PREFIX}_sensitive" \
    --rg "SM:sample" \
    --rg "PL:illumina" \
    --rg "LB:lib1" \
    -x "$REFERENCE_PREFIX" \
    -1 "$READ1" \
    -2 "$READ2" \
    -S results/${OUTPUT_PREFIX}_sensitive.sam \
    2> logs/${OUTPUT_PREFIX}_sensitive.log

BT2_SENS_END=$(date +%s)
BT2_SENS_TIME=$((BT2_SENS_END - BT2_SENS_START))

echo "敏感模式比对完成，用时: ${BT2_SENS_TIME}秒"
echo "输出文件: results/${OUTPUT_PREFIX}_sensitive.sam"

# 3. Bowtie2高敏感模式
echo ""
echo "3. Bowtie2高敏感模式比对..."
echo "----------------------------------------"

BT2_VSENS_START=$(date +%s)

bowtie2 \
    --very-sensitive \
    -p $THREADS \
    --rg-id "${OUTPUT_PREFIX}_very_sensitive" \
    --rg "SM:sample" \
    --rg "PL:illumina" \
    --rg "LB:lib1" \
    -x "$REFERENCE_PREFIX" \
    -1 "$READ1" \
    -2 "$READ2" \
    -S results/${OUTPUT_PREFIX}_very_sensitive.sam \
    2> logs/${OUTPUT_PREFIX}_very_sensitive.log

BT2_VSENS_END=$(date +%s)
BT2_VSENS_TIME=$((BT2_VSENS_END - BT2_VSENS_START))

echo "高敏感模式比对完成，用时: ${BT2_VSENS_TIME}秒"
echo "输出文件: results/${OUTPUT_PREFIX}_very_sensitive.sam"

# 4. Bowtie2局部比对模式
echo ""
echo "4. Bowtie2局部比对模式..."
echo "----------------------------------------"

BT2_LOCAL_START=$(date +%s)

bowtie2 \
    --local \
    --very-sensitive-local \
    -p $THREADS \
    --rg-id "${OUTPUT_PREFIX}_local" \
    --rg "SM:sample" \
    --rg "PL:illumina" \
    --rg "LB:lib1" \
    -x "$REFERENCE_PREFIX" \
    -1 "$READ1" \
    -2 "$READ2" \
    -S results/${OUTPUT_PREFIX}_local.sam \
    2> logs/${OUTPUT_PREFIX}_local.log

BT2_LOCAL_END=$(date +%s)
BT2_LOCAL_TIME=$((BT2_LOCAL_END - BT2_LOCAL_START))

echo "局部比对模式完成，用时: ${BT2_LOCAL_TIME}秒"
echo "输出文件: results/${OUTPUT_PREFIX}_local.sam"

# 5. Bowtie2自定义参数模式
echo ""
echo "5. Bowtie2自定义参数模式..."
echo "----------------------------------------"

BT2_CUSTOM_START=$(date +%s)

bowtie2 \
    -p $THREADS \
    -L 20 \
    -i S,1,0.50 \
    --mp 6,2 \
    --np 1 \
    --rdg 5,3 \
    --rfg 5,3 \
    -I 50 \
    -X 800 \
    --no-mixed \
    --no-discordant \
    --rg-id "${OUTPUT_PREFIX}_custom" \
    --rg "SM:sample" \
    --rg "PL:illumina" \
    --rg "LB:lib1" \
    -x "$REFERENCE_PREFIX" \
    -1 "$READ1" \
    -2 "$READ2" \
    -S results/${OUTPUT_PREFIX}_custom.sam \
    2> logs/${OUTPUT_PREFIX}_custom.log

BT2_CUSTOM_END=$(date +%s)
BT2_CUSTOM_TIME=$((BT2_CUSTOM_END - BT2_CUSTOM_START))

echo "自定义参数模式完成，用时: ${BT2_CUSTOM_TIME}秒"
echo "输出文件: results/${OUTPUT_PREFIX}_custom.sam"

# 6. 转换SAM为BAM并排序
echo ""
echo "6. 转换SAM为BAM格式并排序..."
echo "----------------------------------------"

for mode in fast sensitive very_sensitive local custom; do
    echo "处理 ${mode} 模式结果..."
    
    SAM_FILE="results/${OUTPUT_PREFIX}_${mode}.sam"
    BAM_FILE="results/${OUTPUT_PREFIX}_${mode}_sorted.bam"
    
    # 转换为BAM并排序
    samtools view -bS "$SAM_FILE" | \
    samtools sort -@ $((THREADS-1)) -o "$BAM_FILE"
    
    # 建立索引
    samtools index "$BAM_FILE"
    
    echo "  ✓ $BAM_FILE"
done

# 7. 生成比对统计
echo ""
echo "7. 生成比对统计..."
echo "----------------------------------------"

for mode in fast sensitive very_sensitive local custom; do
    BAM_FILE="results/${OUTPUT_PREFIX}_${mode}_sorted.bam"
    STATS_FILE="results/${OUTPUT_PREFIX}_${mode}_stats.txt"
    
    echo "统计 ${mode} 模式结果..."
    
    # 基本统计
    samtools flagstat "$BAM_FILE" > "$STATS_FILE"
    
    # 详细统计
    samtools stats "$BAM_FILE" >> "$STATS_FILE"
    
    echo "  ✓ $STATS_FILE"
done

# 8. 提取关键统计信息
echo ""
echo "8. 提取关键统计信息..."
echo "----------------------------------------"

cat > results/${OUTPUT_PREFIX}_summary.txt << EOF
Bowtie2比对结果摘要
生成时间: $(date)

比对参数:
- 线程数: $THREADS
- 索引前缀: $REFERENCE_PREFIX
- 测序数据: $READ1, $READ2

运行时间:
- 快速模式: ${BT2_FAST_TIME}秒
- 敏感模式: ${BT2_SENS_TIME}秒
- 高敏感模式: ${BT2_VSENS_TIME}秒
- 局部比对模式: ${BT2_LOCAL_TIME}秒
- 自定义参数模式: ${BT2_CUSTOM_TIME}秒

比对统计:
EOF

for mode in fast sensitive very_sensitive local custom; do
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

# 9. 分析比对日志中的详细信息
echo ""
echo "9. 分析比对日志信息..."
echo "----------------------------------------"

cat > results/${OUTPUT_PREFIX}_alignment_details.txt << EOF
Bowtie2比对详细信息
生成时间: $(date)

EOF

for mode in fast sensitive very_sensitive local custom; do
    LOG_FILE="logs/${OUTPUT_PREFIX}_${mode}.log"
    
    echo "" >> results/${OUTPUT_PREFIX}_alignment_details.txt
    echo "${mode} 模式详细信息:" >> results/${OUTPUT_PREFIX}_alignment_details.txt
    
    if [ -f "$LOG_FILE" ]; then
        # 提取比对率信息
        grep "overall alignment rate" "$LOG_FILE" >> results/${OUTPUT_PREFIX}_alignment_details.txt 2>/dev/null || echo "  未找到比对率信息" >> results/${OUTPUT_PREFIX}_alignment_details.txt
        
        # 提取配对信息
        grep "aligned concordantly" "$LOG_FILE" >> results/${OUTPUT_PREFIX}_alignment_details.txt 2>/dev/null || true
        grep "aligned discordantly" "$LOG_FILE" >> results/${OUTPUT_PREFIX}_alignment_details.txt 2>/dev/null || true
    else
        echo "  日志文件不存在: $LOG_FILE" >> results/${OUTPUT_PREFIX}_alignment_details.txt
    fi
done

# 10. 生成MAPQ分布统计
echo ""
echo "10. 生成MAPQ分布统计..."
echo "----------------------------------------"

for mode in fast sensitive very_sensitive local custom; do
    BAM_FILE="results/${OUTPUT_PREFIX}_${mode}_sorted.bam"
    MAPQ_FILE="results/${OUTPUT_PREFIX}_${mode}_mapq.txt"
    
    echo "分析 ${mode} 模式MAPQ分布..."
    
    # 提取MAPQ值并统计分布
    samtools view "$BAM_FILE" | cut -f5 | sort -n | uniq -c > "$MAPQ_FILE"
    
    echo "  ✓ $MAPQ_FILE"
done

# 11. 计算覆盖度
echo ""
echo "11. 计算基因组覆盖度..."
echo "----------------------------------------"

for mode in fast sensitive very_sensitive local custom; do
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
echo "Bowtie2比对分析完成！"
echo "=========================================="
echo "总用时: ${TOTAL_TIME}秒"
echo ""
echo "生成的文件:"
echo "- SAM文件: results/${OUTPUT_PREFIX}_*.sam"
echo "- BAM文件: results/${OUTPUT_PREFIX}_*_sorted.bam"
echo "- 统计文件: results/${OUTPUT_PREFIX}_*_stats.txt"
echo "- 摘要文件: results/${OUTPUT_PREFIX}_summary.txt"
echo "- 详细信息: results/${OUTPUT_PREFIX}_alignment_details.txt"
echo "- MAPQ分布: results/${OUTPUT_PREFIX}_*_mapq.txt"
echo "- 覆盖度文件: results/${OUTPUT_PREFIX}_*_coverage.txt"
echo ""
echo "可以使用以下命令查看结果:"
echo "  cat results/${OUTPUT_PREFIX}_summary.txt"
echo "  cat results/${OUTPUT_PREFIX}_alignment_details.txt"
echo "  samtools view results/${OUTPUT_PREFIX}_sensitive_sorted.bam | head"

# 记录到日志
cat >> logs/experiment.log << EOF

Bowtie2比对完成时间: $(date)
总用时: ${TOTAL_TIME}秒
快速模式用时: ${BT2_FAST_TIME}秒
敏感模式用时: ${BT2_SENS_TIME}秒
高敏感模式用时: ${BT2_VSENS_TIME}秒
局部比对模式用时: ${BT2_LOCAL_TIME}秒
自定义参数模式用时: ${BT2_CUSTOM_TIME}秒
EOF

echo ""
echo "实验日志已更新: logs/experiment.log"