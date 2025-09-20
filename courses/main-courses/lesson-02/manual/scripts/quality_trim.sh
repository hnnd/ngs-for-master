#!/bin/bash

# 数据清洗脚本 - Trimmomatic质量修剪
# 课程：高通量测序数据分析 - 测序数据质量控制与预处理
# 作者：王运生
# 日期：2025-01-20
# 用法：bash quality_trim.sh [input_dir] [output_dir] [adapter_file] [threads]

set -e  # 遇到错误立即退出

# 默认参数设置
DEFAULT_INPUT_DIR="data"
DEFAULT_OUTPUT_DIR="results/trimmed"
DEFAULT_ADAPTER_FILE="data/TruSeq3-PE.fa"
DEFAULT_THREADS=4

# 参数解析
INPUT_DIR=${1:-$DEFAULT_INPUT_DIR}
OUTPUT_DIR=${2:-$DEFAULT_OUTPUT_DIR}
ADAPTER_FILE=${3:-$DEFAULT_ADAPTER_FILE}
THREADS=${4:-$DEFAULT_THREADS}

echo "=========================================="
echo "Trimmomatic数据清洗脚本"
echo "=========================================="
echo "输入目录: $INPUT_DIR"
echo "输出目录: $OUTPUT_DIR"
echo "接头文件: $ADAPTER_FILE"
echo "线程数: $THREADS"
echo "=========================================="

# 检查输入目录
if [ ! -d "$INPUT_DIR" ]; then
    echo "错误：输入目录不存在: $INPUT_DIR"
    exit 1
fi

# 创建输出目录
mkdir -p "$OUTPUT_DIR"
mkdir -p logs

# 检查接头文件
if [ ! -f "$ADAPTER_FILE" ]; then
    echo "警告：接头文件不存在: $ADAPTER_FILE"
    echo "创建默认的TruSeq3-PE接头文件..."
    
    mkdir -p "$(dirname "$ADAPTER_FILE")"
    cat > "$ADAPTER_FILE" << 'EOF'
>PrefixPE/1
TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PrefixPE/2
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
>PE1
TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PE1_rc
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
>PE2
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
>PE2_rc
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
EOF
    echo "✓ 默认接头文件已创建"
fi

# 检查Trimmomatic是否可用
check_trimmomatic() {
    if command -v trimmomatic &> /dev/null; then
        echo "✓ 使用系统安装的Trimmomatic"
        TRIMMOMATIC_CMD="trimmomatic"
        return 0
    elif [ -f "trimmomatic-0.39.jar" ]; then
        echo "✓ 使用本地Trimmomatic JAR文件"
        TRIMMOMATIC_CMD="java -jar trimmomatic-0.39.jar"
        return 0
    elif [ -f "~/ngs-analysis/lesson-02/trimmomatic-0.39.jar" ]; then
        echo "✓ 使用工作目录中的Trimmomatic"
        TRIMMOMATIC_CMD="java -jar ~/ngs-analysis/lesson-02/trimmomatic-0.39.jar"
        return 0
    else
        echo "错误：找不到Trimmomatic"
        echo "请运行setup.sh安装必要软件"
        return 1
    fi
}

if ! check_trimmomatic; then
    exit 1
fi

# 查找配对的FASTQ文件
echo "查找配对的FASTQ文件..."
declare -A R1_FILES
declare -A R2_FILES

# 查找R1文件
for file in "$INPUT_DIR"/*_R1.fastq* "$INPUT_DIR"/*_1.fastq*; do
    if [ -f "$file" ]; then
        basename_file=$(basename "$file")
        # 提取样本名称
        if [[ "$basename_file" =~ (.+)_R1\.(fastq|fq)(\.gz)?$ ]]; then
            sample_name="${BASH_REMATCH[1]}"
            R1_FILES["$sample_name"]="$file"
        elif [[ "$basename_file" =~ (.+)_1\.(fastq|fq)(\.gz)?$ ]]; then
            sample_name="${BASH_REMATCH[1]}"
            R1_FILES["$sample_name"]="$file"
        fi
    fi
done

# 查找对应的R2文件
for sample_name in "${!R1_FILES[@]}"; do
    r1_file="${R1_FILES[$sample_name]}"
    r1_basename=$(basename "$r1_file")
    
    # 尝试不同的R2文件命名模式
    if [[ "$r1_basename" =~ (.+)_R1\.(fastq|fq)(\.gz)?$ ]]; then
        r2_pattern="${BASH_REMATCH[1]}_R2.${BASH_REMATCH[2]}${BASH_REMATCH[3]}"
    elif [[ "$r1_basename" =~ (.+)_1\.(fastq|fq)(\.gz)?$ ]]; then
        r2_pattern="${BASH_REMATCH[1]}_2.${BASH_REMATCH[2]}${BASH_REMATCH[3]}"
    else
        continue
    fi
    
    r2_file="$INPUT_DIR/$r2_pattern"
    if [ -f "$r2_file" ]; then
        R2_FILES["$sample_name"]="$r2_file"
    else
        echo "警告：找不到 $sample_name 的R2文件: $r2_file"
        unset R1_FILES["$sample_name"]
    fi
done

# 检查是否找到配对文件
if [ ${#R1_FILES[@]} -eq 0 ]; then
    echo "错误：未找到配对的FASTQ文件"
    echo "支持的命名格式："
    echo "  - sample_R1.fastq.gz 和 sample_R2.fastq.gz"
    echo "  - sample_1.fastq.gz 和 sample_2.fastq.gz"
    exit 1
fi

echo "找到 ${#R1_FILES[@]} 对配对文件："
for sample_name in "${!R1_FILES[@]}"; do
    echo "  样本: $sample_name"
    echo "    R1: $(basename "${R1_FILES[$sample_name]}")"
    echo "    R2: $(basename "${R2_FILES[$sample_name]}")"
done

# 定义Trimmomatic参数
TRIMMOMATIC_PARAMS="ILLUMINACLIP:${ADAPTER_FILE}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

echo "Trimmomatic参数: $TRIMMOMATIC_PARAMS"

# 创建参数说明文件
cat > "$OUTPUT_DIR/trimming_parameters.txt" << EOF
Trimmomatic参数说明
==================

使用的参数: $TRIMMOMATIC_PARAMS

参数详解:
- ILLUMINACLIP:${ADAPTER_FILE}:2:30:10
  * 接头文件: ${ADAPTER_FILE}
  * 种子错配数: 2 (允许的错配数)
  * 回文阈值: 30 (回文匹配的阈值)
  * 简单clip阈值: 10 (简单匹配的阈值)

- LEADING:3
  * 去除5'端质量低于3的碱基

- TRAILING:3
  * 去除3'端质量低于3的碱基

- SLIDINGWINDOW:4:15
  * 滑动窗口大小: 4bp
  * 平均质量阈值: 15 (Phred分数)

- MINLEN:36
  * 保留长度至少36bp的序列

质量阈值说明:
- Q15: 错误率 ~3.2% (96.8%准确率)
- Q20: 错误率 1% (99%准确率)
- Q30: 错误率 0.1% (99.9%准确率)

生成时间: $(date)
EOF

# 开始处理文件
echo "开始数据清洗..."
START_TIME=$(date +%s)

TOTAL_SAMPLES=${#R1_FILES[@]}
CURRENT_SAMPLE=0
SUCCESSFUL_SAMPLES=0
FAILED_SAMPLES=0

for sample_name in "${!R1_FILES[@]}"; do
    CURRENT_SAMPLE=$((CURRENT_SAMPLE + 1))
    echo ""
    echo "处理样本 $CURRENT_SAMPLE/$TOTAL_SAMPLES: $sample_name"
    echo "----------------------------------------"
    
    R1_INPUT="${R1_FILES[$sample_name]}"
    R2_INPUT="${R2_FILES[$sample_name]}"
    
    # 输出文件名
    R1_PAIRED="$OUTPUT_DIR/${sample_name}_R1_paired.fastq.gz"
    R1_UNPAIRED="$OUTPUT_DIR/${sample_name}_R1_unpaired.fastq.gz"
    R2_PAIRED="$OUTPUT_DIR/${sample_name}_R2_paired.fastq.gz"
    R2_UNPAIRED="$OUTPUT_DIR/${sample_name}_R2_unpaired.fastq.gz"
    
    # 日志文件
    LOG_FILE="logs/trimmomatic_${sample_name}.log"
    
    echo "输入文件:"
    echo "  R1: $(basename "$R1_INPUT")"
    echo "  R2: $(basename "$R2_INPUT")"
    echo "输出文件:"
    echo "  R1_paired: $(basename "$R1_PAIRED")"
    echo "  R2_paired: $(basename "$R2_PAIRED")"
    
    # 构建Trimmomatic命令
    TRIM_CMD="$TRIMMOMATIC_CMD PE -threads $THREADS -phred33 \
        \"$R1_INPUT\" \"$R2_INPUT\" \
        \"$R1_PAIRED\" \"$R1_UNPAIRED\" \
        \"$R2_PAIRED\" \"$R2_UNPAIRED\" \
        $TRIMMOMATIC_PARAMS"
    
    echo "执行命令:"
    echo "$TRIM_CMD"
    
    # 记录开始时间
    SAMPLE_START=$(date +%s)
    
    # 执行Trimmomatic
    if eval "$TRIM_CMD" 2> "$LOG_FILE"; then
        SAMPLE_END=$(date +%s)
        SAMPLE_DURATION=$((SAMPLE_END - SAMPLE_START))
        
        echo "✓ 样本 $sample_name 处理成功 (耗时: ${SAMPLE_DURATION}秒)"
        SUCCESSFUL_SAMPLES=$((SUCCESSFUL_SAMPLES + 1))
        
        # 统计处理结果
        if [ -f "$R1_PAIRED" ] && [ -f "$R2_PAIRED" ]; then
            # 计算保留率
            if command -v zcat &> /dev/null; then
                INPUT_READS=$(zcat "$R1_INPUT" | wc -l)
                INPUT_READS=$((INPUT_READS / 4))
                
                OUTPUT_READS=$(zcat "$R1_PAIRED" | wc -l)
                OUTPUT_READS=$((OUTPUT_READS / 4))
                
                if [ $INPUT_READS -gt 0 ]; then
                    RETENTION_RATE=$(echo "scale=2; $OUTPUT_READS * 100 / $INPUT_READS" | bc 2>/dev/null || echo "N/A")
                    echo "  输入序列: $INPUT_READS"
                    echo "  输出序列: $OUTPUT_READS"
                    echo "  保留率: ${RETENTION_RATE}%"
                fi
            fi
            
            # 文件大小比较
            INPUT_SIZE=$(ls -lh "$R1_INPUT" | awk '{print $5}')
            OUTPUT_SIZE=$(ls -lh "$R1_PAIRED" | awk '{print $5}')
            echo "  文件大小: $INPUT_SIZE -> $OUTPUT_SIZE"
        fi
        
    else
        echo "✗ 样本 $sample_name 处理失败"
        echo "  错误日志: $LOG_FILE"
        FAILED_SAMPLES=$((FAILED_SAMPLES + 1))
        
        # 显示错误信息
        if [ -f "$LOG_FILE" ]; then
            echo "  错误详情:"
            tail -5 "$LOG_FILE" | sed 's/^/    /'
        fi
    fi
done

# 计算总耗时
END_TIME=$(date +%s)
TOTAL_DURATION=$((END_TIME - START_TIME))

echo ""
echo "=========================================="
echo "数据清洗完成！"
echo "=========================================="
echo "处理统计:"
echo "  总样本数: $TOTAL_SAMPLES"
echo "  成功处理: $SUCCESSFUL_SAMPLES"
echo "  处理失败: $FAILED_SAMPLES"
echo "  总耗时: ${TOTAL_DURATION}秒"
echo ""
echo "输出目录: $OUTPUT_DIR"
echo "日志目录: logs/"

# 生成处理摘要
SUMMARY_FILE="$OUTPUT_DIR/trimming_summary.txt"
cat > "$SUMMARY_FILE" << EOF
Trimmomatic数据清洗摘要
======================

处理时间: $(date)
输入目录: $INPUT_DIR
输出目录: $OUTPUT_DIR
接头文件: $ADAPTER_FILE
使用线程: $THREADS

处理统计:
- 总样本数: $TOTAL_SAMPLES
- 成功处理: $SUCCESSFUL_SAMPLES
- 处理失败: $FAILED_SAMPLES
- 总耗时: ${TOTAL_DURATION}秒

使用参数:
$TRIMMOMATIC_PARAMS

输出文件说明:
- *_R1_paired.fastq.gz: 清洗后的R1配对序列
- *_R2_paired.fastq.gz: 清洗后的R2配对序列
- *_R1_unpaired.fastq.gz: 清洗后的R1非配对序列
- *_R2_unpaired.fastq.gz: 清洗后的R2非配对序列

样本处理详情:
EOF

for sample_name in "${!R1_FILES[@]}"; do
    echo "样本: $sample_name" >> "$SUMMARY_FILE"
    LOG_FILE="logs/trimmomatic_${sample_name}.log"
    if [ -f "$LOG_FILE" ]; then
        # 提取Trimmomatic的统计信息
        if grep -q "Input Read Pairs" "$LOG_FILE"; then
            grep "Input Read Pairs\|Both Surviving\|Forward Only Surviving\|Reverse Only Surviving\|Dropped" "$LOG_FILE" | sed 's/^/  /' >> "$SUMMARY_FILE"
        else
            echo "  日志文件存在但无统计信息" >> "$SUMMARY_FILE"
        fi
    else
        echo "  处理失败或日志文件缺失" >> "$SUMMARY_FILE"
    fi
    echo "" >> "$SUMMARY_FILE"
done

echo "处理摘要: $SUMMARY_FILE"

# 创建后续分析脚本
cat > "$OUTPUT_DIR/run_fastqc_clean.sh" << 'EOF'
#!/bin/bash
# 对清洗后的数据运行FastQC分析

echo "对清洗后的数据运行FastQC分析..."

# 创建输出目录
mkdir -p ../fastqc_clean

# 只分析配对的序列文件
fastqc -t 4 -o ../fastqc_clean/ *_paired.fastq.gz

echo "FastQC分析完成，结果保存在 ../fastqc_clean/"
echo "可以运行MultiQC进行整合分析："
echo "multiqc ../fastqc_clean/ -o ../multiqc_clean/"
EOF

chmod +x "$OUTPUT_DIR/run_fastqc_clean.sh"

echo "后续分析脚本: $OUTPUT_DIR/run_fastqc_clean.sh"

if [ $FAILED_SAMPLES -eq 0 ]; then
    echo ""
    echo "✓ 所有样本处理成功！"
    echo "建议下一步:"
    echo "1. 运行 $OUTPUT_DIR/run_fastqc_clean.sh 分析清洗后数据"
    echo "2. 比较清洗前后的质量差异"
    echo "3. 检查保留率是否在合理范围内"
else
    echo ""
    echo "⚠ 有 $FAILED_SAMPLES 个样本处理失败"
    echo "请检查相应的日志文件排查问题"
fi

echo "=========================================="