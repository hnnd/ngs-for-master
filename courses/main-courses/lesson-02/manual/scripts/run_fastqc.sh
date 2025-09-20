#!/bin/bash

# FastQC批量运行脚本
# 课程：高通量测序数据分析 - 测序数据质量控制与预处理
# 作者：王运生
# 日期：2025-01-20
# 用法：bash run_fastqc.sh [input_dir] [output_dir] [threads]

set -e  # 遇到错误立即退出

# 默认参数设置
DEFAULT_INPUT_DIR="data"
DEFAULT_OUTPUT_DIR="results/fastqc_raw"
DEFAULT_THREADS=4

# 参数解析
INPUT_DIR=${1:-$DEFAULT_INPUT_DIR}
OUTPUT_DIR=${2:-$DEFAULT_OUTPUT_DIR}
THREADS=${3:-$DEFAULT_THREADS}

echo "=========================================="
echo "FastQC批量质量分析脚本"
echo "=========================================="
echo "输入目录: $INPUT_DIR"
echo "输出目录: $OUTPUT_DIR"
echo "线程数: $THREADS"
echo "=========================================="

# 检查输入目录
if [ ! -d "$INPUT_DIR" ]; then
    echo "错误：输入目录不存在: $INPUT_DIR"
    exit 1
fi

# 创建输出目录
mkdir -p "$OUTPUT_DIR"

# 检查FastQC是否安装
if ! command -v fastqc &> /dev/null; then
    echo "错误：FastQC未安装或不在PATH中"
    echo "请先运行setup.sh安装必要软件"
    exit 1
fi

# 查找FASTQ文件
echo "查找FASTQ文件..."
FASTQ_FILES=($(find "$INPUT_DIR" -name "*.fastq" -o -name "*.fastq.gz" -o -name "*.fq" -o -name "*.fq.gz"))

if [ ${#FASTQ_FILES[@]} -eq 0 ]; then
    echo "错误：在 $INPUT_DIR 中未找到FASTQ文件"
    echo "支持的文件格式：*.fastq, *.fastq.gz, *.fq, *.fq.gz"
    exit 1
fi

echo "找到 ${#FASTQ_FILES[@]} 个FASTQ文件："
for file in "${FASTQ_FILES[@]}"; do
    echo "  - $(basename "$file")"
done

# 运行FastQC分析
echo "开始FastQC分析..."
START_TIME=$(date +%s)

# 创建临时日志文件
LOG_FILE="logs/fastqc_$(date +%Y%m%d_%H%M%S).log"
mkdir -p logs

echo "FastQC分析日志 - $(date)" > "$LOG_FILE"
echo "输入目录: $INPUT_DIR" >> "$LOG_FILE"
echo "输出目录: $OUTPUT_DIR" >> "$LOG_FILE"
echo "线程数: $THREADS" >> "$LOG_FILE"
echo "文件列表:" >> "$LOG_FILE"
for file in "${FASTQ_FILES[@]}"; do
    echo "  - $file" >> "$LOG_FILE"
done
echo "========================================" >> "$LOG_FILE"

# 执行FastQC
echo "运行FastQC命令..."
fastqc -t "$THREADS" -o "$OUTPUT_DIR" "${FASTQ_FILES[@]}" 2>&1 | tee -a "$LOG_FILE"

# 检查FastQC是否成功完成
if [ $? -eq 0 ]; then
    echo "✓ FastQC分析成功完成"
else
    echo "✗ FastQC分析失败，请检查日志文件: $LOG_FILE"
    exit 1
fi

# 计算运行时间
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
echo "分析耗时: ${DURATION}秒"

# 检查输出文件
echo "检查输出文件..."
HTML_FILES=($(find "$OUTPUT_DIR" -name "*.html"))
ZIP_FILES=($(find "$OUTPUT_DIR" -name "*_fastqc.zip"))

echo "生成的报告文件:"
echo "HTML报告: ${#HTML_FILES[@]} 个"
for file in "${HTML_FILES[@]}"; do
    echo "  - $(basename "$file")"
done

echo "ZIP数据包: ${#ZIP_FILES[@]} 个"
for file in "${ZIP_FILES[@]}"; do
    echo "  - $(basename "$file")"
done

# 生成分析摘要
echo "生成分析摘要..."
SUMMARY_FILE="$OUTPUT_DIR/fastqc_summary.txt"

cat > "$SUMMARY_FILE" << EOF
FastQC分析摘要报告
==================

分析时间: $(date)
输入目录: $INPUT_DIR
输出目录: $OUTPUT_DIR
使用线程: $THREADS
分析耗时: ${DURATION}秒

文件统计:
- 输入FASTQ文件: ${#FASTQ_FILES[@]} 个
- 生成HTML报告: ${#HTML_FILES[@]} 个
- 生成ZIP数据包: ${#ZIP_FILES[@]} 个

输入文件列表:
EOF

for file in "${FASTQ_FILES[@]}"; do
    FILE_SIZE=$(ls -lh "$file" | awk '{print $5}')
    echo "- $(basename "$file") ($FILE_SIZE)" >> "$SUMMARY_FILE"
done

cat >> "$SUMMARY_FILE" << EOF

输出文件列表:
EOF

for file in "${HTML_FILES[@]}"; do
    echo "- $(basename "$file")" >> "$SUMMARY_FILE"
done

# 提取关键质量指标
echo "提取质量指标摘要..."
METRICS_FILE="$OUTPUT_DIR/quality_metrics.txt"

cat > "$METRICS_FILE" << EOF
质量指标摘要
============

EOF

for zip_file in "${ZIP_FILES[@]}"; do
    SAMPLE_NAME=$(basename "$zip_file" _fastqc.zip)
    echo "样本: $SAMPLE_NAME" >> "$METRICS_FILE"
    
    # 提取基本统计信息
    if command -v unzip &> /dev/null; then
        FASTQC_DATA=$(unzip -p "$zip_file" "*/fastqc_data.txt" 2>/dev/null || echo "")
        
        if [ -n "$FASTQC_DATA" ]; then
            # 提取总序列数
            TOTAL_SEQS=$(echo "$FASTQC_DATA" | grep "Total Sequences" | cut -f2)
            # 提取GC含量
            GC_CONTENT=$(echo "$FASTQC_DATA" | grep "%GC" | cut -f2)
            # 提取序列长度
            SEQ_LENGTH=$(echo "$FASTQC_DATA" | grep "Sequence length" | cut -f2)
            
            echo "  总序列数: $TOTAL_SEQS" >> "$METRICS_FILE"
            echo "  GC含量: $GC_CONTENT%" >> "$METRICS_FILE"
            echo "  序列长度: $SEQ_LENGTH" >> "$METRICS_FILE"
            
            # 检查模块状态
            echo "  模块状态:" >> "$METRICS_FILE"
            echo "$FASTQC_DATA" | grep -E "(PASS|WARN|FAIL)" | while read line; do
                MODULE=$(echo "$line" | cut -f1)
                STATUS=$(echo "$line" | cut -f2)
                echo "    $MODULE: $STATUS" >> "$METRICS_FILE"
            done
        else
            echo "  无法提取详细信息" >> "$METRICS_FILE"
        fi
    else
        echo "  需要unzip命令来提取详细信息" >> "$METRICS_FILE"
    fi
    
    echo "" >> "$METRICS_FILE"
done

# 创建快速查看脚本
cat > "$OUTPUT_DIR/view_reports.sh" << 'EOF'
#!/bin/bash
# 快速查看FastQC报告的脚本

echo "可用的FastQC报告:"
ls -1 *.html 2>/dev/null | nl

echo ""
read -p "请选择要查看的报告编号 (或按Enter查看所有): " choice

if [ -z "$choice" ]; then
    echo "在浏览器中打开所有报告..."
    for html in *.html; do
        if command -v firefox &> /dev/null; then
            firefox "$html" &
        elif command -v google-chrome &> /dev/null; then
            google-chrome "$html" &
        elif command -v chromium-browser &> /dev/null; then
            chromium-browser "$html" &
        else
            echo "请手动打开: $html"
        fi
    done
else
    HTML_FILE=$(ls -1 *.html 2>/dev/null | sed -n "${choice}p")
    if [ -n "$HTML_FILE" ]; then
        echo "打开报告: $HTML_FILE"
        if command -v firefox &> /dev/null; then
            firefox "$HTML_FILE" &
        elif command -v google-chrome &> /dev/null; then
            google-chrome "$HTML_FILE" &
        elif command -v chromium-browser &> /dev/null; then
            chromium-browser "$HTML_FILE" &
        else
            echo "请手动打开: $HTML_FILE"
        fi
    else
        echo "无效的选择"
    fi
fi
EOF

chmod +x "$OUTPUT_DIR/view_reports.sh"

# 最终报告
echo "=========================================="
echo "FastQC分析完成！"
echo "=========================================="
echo "输出目录: $OUTPUT_DIR"
echo "分析摘要: $SUMMARY_FILE"
echo "质量指标: $METRICS_FILE"
echo "查看报告: $OUTPUT_DIR/view_reports.sh"
echo ""
echo "下一步建议："
echo "1. 查看HTML报告了解数据质量"
echo "2. 运行MultiQC整合多个样本的报告"
echo "3. 根据质量问题制定数据清洗策略"
echo "=========================================="

# 记录完成信息到日志
echo "========================================" >> "$LOG_FILE"
echo "FastQC分析完成 - $(date)" >> "$LOG_FILE"
echo "分析耗时: ${DURATION}秒" >> "$LOG_FILE"
echo "生成文件: ${#HTML_FILES[@]} 个HTML报告, ${#ZIP_FILES[@]} 个ZIP数据包" >> "$LOG_FILE"