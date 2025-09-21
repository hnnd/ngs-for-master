#!/bin/bash
# 单细胞RNA测序数据下载脚本
# 课程：高通量测序数据分析 - 第7次课
# 作者：王运生
# 日期：2025-01-21
# 用法：bash scripts/download_data.sh

set -e  # 遇到错误立即退出
set -u  # 使用未定义变量时报错

echo "=== 单细胞RNA测序数据下载脚本 ==="

# 检查必要的工具
check_command() {
    if ! command -v $1 &> /dev/null; then
        echo "错误: $1 命令未找到，请先安装"
        exit 1
    fi
}

echo "检查必要工具..."
check_command wget
check_command tar
echo "工具检查完成"

# 创建数据目录
echo "创建数据目录..."
mkdir -p data
cd data

# PBMC 3K数据集信息
DATASET_NAME="pbmc3k"
DATASET_URL="https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
DATASET_FILE="pbmc3k_filtered_gene_bc_matrices.tar.gz"
EXPECTED_SIZE="20971520"  # 约20MB

echo "准备下载PBMC 3K数据集..."
echo "数据集: $DATASET_NAME"
echo "URL: $DATASET_URL"

# 检查文件是否已存在
if [ -f "$DATASET_FILE" ]; then
    echo "数据文件已存在: $DATASET_FILE"
    
    # 检查文件大小
    FILE_SIZE=$(stat -f%z "$DATASET_FILE" 2>/dev/null || stat -c%s "$DATASET_FILE" 2>/dev/null || echo "0")
    
    if [ "$FILE_SIZE" -gt 10000000 ]; then  # 大于10MB认为下载完整
        echo "文件大小正常: $(($FILE_SIZE / 1024 / 1024)) MB"
        echo "跳过下载步骤"
    else
        echo "文件大小异常，重新下载..."
        rm -f "$DATASET_FILE"
    fi
fi

# 下载数据
if [ ! -f "$DATASET_FILE" ]; then
    echo "开始下载数据..."
    echo "这可能需要几分钟时间，请耐心等待..."
    
    # 使用wget下载，显示进度条
    wget --progress=bar:force:noscroll -O "$DATASET_FILE" "$DATASET_URL"
    
    if [ $? -eq 0 ]; then
        echo "数据下载完成"
    else
        echo "数据下载失败"
        rm -f "$DATASET_FILE"
        exit 1
    fi
    
    # 验证下载的文件
    FILE_SIZE=$(stat -f%z "$DATASET_FILE" 2>/dev/null || stat -c%s "$DATASET_FILE" 2>/dev/null || echo "0")
    echo "下载文件大小: $(($FILE_SIZE / 1024 / 1024)) MB"
    
    if [ "$FILE_SIZE" -lt 10000000 ]; then
        echo "警告: 文件大小可能不正确，请检查网络连接"
    fi
fi

# 解压数据
echo "解压数据文件..."
if [ -d "filtered_gene_bc_matrices" ]; then
    echo "数据已解压，跳过解压步骤"
else
    tar -xzf "$DATASET_FILE"
    
    if [ $? -eq 0 ]; then
        echo "数据解压完成"
    else
        echo "数据解压失败"
        exit 1
    fi
fi

# 验证解压结果
echo "验证数据文件..."
DATA_DIR="filtered_gene_bc_matrices/hg19"

if [ -d "$DATA_DIR" ]; then
    echo "数据目录存在: $DATA_DIR"
    
    # 检查必要文件
    REQUIRED_FILES=("matrix.mtx" "barcodes.tsv" "features.tsv")
    
    # 兼容旧版本文件名
    if [ ! -f "$DATA_DIR/features.tsv" ] && [ -f "$DATA_DIR/genes.tsv" ]; then
        echo "发现旧版本基因文件，创建符号链接..."
        ln -sf genes.tsv "$DATA_DIR/features.tsv"
    fi
    
    for file in "${REQUIRED_FILES[@]}"; do
        if [ -f "$DATA_DIR/$file" ]; then
            FILE_SIZE=$(stat -f%z "$DATA_DIR/$file" 2>/dev/null || stat -c%s "$DATA_DIR/$file" 2>/dev/null || echo "0")
            echo "✓ $file ($(($FILE_SIZE / 1024)) KB)"
        else
            echo "✗ 缺少文件: $file"
            exit 1
        fi
    done
    
    echo "所有必要文件都存在"
else
    echo "错误: 数据目录不存在"
    exit 1
fi

# 显示数据统计信息
echo ""
echo "=== 数据集信息 ==="
echo "数据集名称: PBMC 3K"
echo "数据类型: 单细胞RNA测序"
echo "样本来源: 外周血单核细胞(PBMC)"
echo "细胞数量: ~2700"
echo "基因数量: ~32000"
echo "测序平台: 10X Genomics Chromium"
echo "参考基因组: GRCh37 (hg19)"

# 显示文件信息
echo ""
echo "=== 文件说明 ==="
echo "matrix.mtx    - 基因表达矩阵（稀疏矩阵格式）"
echo "barcodes.tsv  - 细胞条码列表"
echo "features.tsv  - 基因信息列表"

# 快速统计
if command -v wc &> /dev/null; then
    echo ""
    echo "=== 快速统计 ==="
    BARCODE_COUNT=$(wc -l < "$DATA_DIR/barcodes.tsv")
    FEATURE_COUNT=$(wc -l < "$DATA_DIR/features.tsv")
    echo "细胞数量: $BARCODE_COUNT"
    echo "基因数量: $FEATURE_COUNT"
fi

# 创建数据描述文件
cat > data_description.txt << EOF
PBMC 3K Dataset Description
===========================

Dataset: Peripheral Blood Mononuclear Cells (PBMC) 3K
Source: 10X Genomics
URL: https://www.10xgenomics.com/resources/datasets/3-k-pbm-cs-from-a-healthy-donor-1-standard-1-1-0

Sample Information:
- Sample type: Fresh frozen human PBMCs
- Donor: Healthy adult
- Cell count: ~2,700 cells
- Sequencing depth: ~69,000 reads per cell
- Chemistry: Single Cell 3' v1

Files:
- matrix.mtx: Gene expression count matrix (Market Matrix format)
- barcodes.tsv: Cell barcodes (one per line)
- features.tsv: Gene information (gene_id, gene_symbol, gene_type)

Reference:
- Genome: GRCh37 (hg19)
- Annotation: Ensembl

Download date: $(date)
Downloaded by: Single-cell RNA-seq Analysis Course
EOF

echo ""
echo "数据描述文件已创建: data_description.txt"

# 返回原目录
cd ..

echo ""
echo "=== 数据下载完成 ==="
echo "数据位置: data/filtered_gene_bc_matrices/hg19/"
echo "现在可以开始R分析了"
echo ""
echo "下一步："
echo "1. 启动RStudio或R控制台"
echo "2. 运行: source('scripts/setup.R')"
echo "3. 运行: source('scripts/complete_analysis.R')"

echo ""
echo "数据下载脚本执行完成！"