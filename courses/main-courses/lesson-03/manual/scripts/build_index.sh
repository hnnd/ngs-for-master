#!/bin/bash
# 构建比对工具索引脚本
# 作者：王运生
# 日期：2025-01-20
# 用法：bash build_index.sh [reference.fa]

set -e  # 遇到错误立即退出
set -u  # 使用未定义变量时报错

echo "=========================================="
echo "构建序列比对索引"
echo "=========================================="

# 参数检查
REFERENCE_FILE=${1:-"reference/chr22_fragment.fa"}

if [ ! -f "$REFERENCE_FILE" ]; then
    echo "错误：参考基因组文件不存在: $REFERENCE_FILE"
    echo "用法: $0 [reference.fa]"
    exit 1
fi

# 获取文件信息
REF_DIR=$(dirname "$REFERENCE_FILE")
REF_BASE=$(basename "$REFERENCE_FILE")
REF_NAME="${REF_BASE%.*}"

echo "参考基因组文件: $REFERENCE_FILE"
echo "文件大小: $(du -h "$REFERENCE_FILE" | cut -f1)"
echo "序列数量: $(grep -c '^>' "$REFERENCE_FILE")"

# 记录开始时间
START_TIME=$(date +%s)
echo "开始时间: $(date)"

# 构建BWA索引
echo ""
echo "1. 构建BWA索引..."
echo "----------------------------------------"

if [ -f "${REFERENCE_FILE}.bwt" ]; then
    echo "BWA索引已存在，跳过构建"
else
    echo "开始构建BWA索引..."
    BWA_START=$(date +%s)
    
    bwa index "$REFERENCE_FILE" 2>&1 | tee logs/bwa_index.log
    
    BWA_END=$(date +%s)
    BWA_TIME=$((BWA_END - BWA_START))
    echo "BWA索引构建完成，用时: ${BWA_TIME}秒"
fi

# 检查BWA索引文件
echo ""
echo "BWA索引文件："
for ext in amb ann bwt pac sa; do
    index_file="${REFERENCE_FILE}.${ext}"
    if [ -f "$index_file" ]; then
        size=$(du -h "$index_file" | cut -f1)
        echo "  ✓ ${index_file} (${size})"
    else
        echo "  ✗ ${index_file} (缺失)"
    fi
done

# 构建Bowtie2索引
echo ""
echo "2. 构建Bowtie2索引..."
echo "----------------------------------------"

BOWTIE2_PREFIX="${REF_DIR}/${REF_NAME}_bt2"

if [ -f "${BOWTIE2_PREFIX}.1.bt2" ]; then
    echo "Bowtie2索引已存在，跳过构建"
else
    echo "开始构建Bowtie2索引..."
    BT2_START=$(date +%s)
    
    bowtie2-build "$REFERENCE_FILE" "$BOWTIE2_PREFIX" 2>&1 | tee logs/bowtie2_index.log
    
    BT2_END=$(date +%s)
    BT2_TIME=$((BT2_END - BT2_START))
    echo "Bowtie2索引构建完成，用时: ${BT2_TIME}秒"
fi

# 检查Bowtie2索引文件
echo ""
echo "Bowtie2索引文件："
for ext in 1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2; do
    index_file="${BOWTIE2_PREFIX}.${ext}"
    if [ -f "$index_file" ]; then
        size=$(du -h "$index_file" | cut -f1)
        echo "  ✓ ${index_file} (${size})"
    else
        echo "  ✗ ${index_file} (缺失)"
    fi
done

# 计算总用时
END_TIME=$(date +%s)
TOTAL_TIME=$((END_TIME - START_TIME))

echo ""
echo "=========================================="
echo "索引构建完成！"
echo "=========================================="
echo "总用时: ${TOTAL_TIME}秒"

# 生成索引信息摘要
echo ""
echo "3. 生成索引信息摘要..."
echo "----------------------------------------"

cat > results/index_summary.txt << EOF
序列比对索引构建摘要
生成时间: $(date)

参考基因组信息:
- 文件: $REFERENCE_FILE
- 大小: $(du -h "$REFERENCE_FILE" | cut -f1)
- 序列数: $(grep -c '^>' "$REFERENCE_FILE")

BWA索引:
- 构建时间: ${BWA_TIME:-"已存在"}秒
- 索引文件大小: $(du -ch "${REFERENCE_FILE}".* 2>/dev/null | tail -1 | cut -f1 || echo "未知")

Bowtie2索引:
- 构建时间: ${BT2_TIME:-"已存在"}秒  
- 索引前缀: $BOWTIE2_PREFIX
- 索引文件大小: $(du -ch "${BOWTIE2_PREFIX}".* 2>/dev/null | tail -1 | cut -f1 || echo "未知")

总构建时间: ${TOTAL_TIME}秒
EOF

echo "索引信息已保存到: results/index_summary.txt"

# 验证索引完整性
echo ""
echo "4. 验证索引完整性..."
echo "----------------------------------------"

# 验证BWA索引
echo "验证BWA索引..."
if bwa 2>&1 | grep -q "index"; then
    echo "✓ BWA索引格式正确"
else
    echo "✗ BWA索引可能有问题"
fi

# 验证Bowtie2索引
echo "验证Bowtie2索引..."
if bowtie2-inspect -s "$BOWTIE2_PREFIX" > /dev/null 2>&1; then
    echo "✓ Bowtie2索引格式正确"
    
    # 获取索引统计信息
    echo "Bowtie2索引统计信息:"
    bowtie2-inspect -s "$BOWTIE2_PREFIX" | head -5
else
    echo "✗ Bowtie2索引可能有问题"
fi

echo ""
echo "索引构建和验证完成！"
echo "现在可以使用这些索引进行序列比对了。"

# 更新配置文件
if [ -f "config.txt" ]; then
    echo ""
    echo "更新配置文件..."
    
    # 添加索引路径到配置文件
    cat >> config.txt << EOF

# 索引文件路径
BWA_INDEX=$REFERENCE_FILE
BOWTIE2_INDEX=$BOWTIE2_PREFIX
EOF
    
    echo "配置文件已更新"
fi