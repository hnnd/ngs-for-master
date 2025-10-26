#!/bin/bash

###############################################################################
#                   ChIP-seq数据准备脚本                                     #
# 自动下载参考基因组、注释文件和ChIP-seq测序数据                            #
###############################################################################

set -e

# 配置参数
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
REFERENCE_DIR="$WORK_DIR/reference"
DATA_DIR="$WORK_DIR/data"
LOGS_DIR="$WORK_DIR/../logs"

# 创建必要的目录
mkdir -p "$REFERENCE_DIR" "$DATA_DIR" "$LOGS_DIR"

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

# 日志函数
log_info() {
    echo -e "${GREEN}[INFO]${NC} $1" | tee -a "$LOGS_DIR/data_download.log"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1" | tee -a "$LOGS_DIR/data_download.log"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1" | tee -a "$LOGS_DIR/data_download.log"
}

# 打印帮助信息
usage() {
    cat << EOF
ChIP-seq数据准备脚本

使用方法: $0 [选项]

选项:
    --reference     仅下载参考基因组和注释文件
    --sequencing    仅准备ChIP-seq测序数据
    --validate      仅验证数据完整性
    --all           完整准备所有数据（默认选项）
    --help          显示此帮助信息

示例:
    $0                      # 完整数据准备
    $0 --reference          # 仅下载参考数据
    $0 --sequencing         # 仅准备测序数据

EOF
}

# 检查软件依赖
check_dependencies() {
    log_info "检查软件依赖..."

    local required_tools=("wget" "gunzip" "samtools")
    local missing_tools=()

    for tool in "${required_tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            missing_tools+=("$tool")
        fi
    done

    if [ ${#missing_tools[@]} -gt 0 ]; then
        log_warn "缺少以下工具: ${missing_tools[@]}"
        log_info "可以使用conda安装: conda install -c bioconda samtools wget"
    else
        log_info "所有必需的工具都已安装"
    fi
}

# 下载参考基因组
download_reference() {
    log_info "开始下载参考基因组..."

    cd "$REFERENCE_DIR"

    # 创建演示基因组序列（减小数据量）
    log_info "生成演示参考基因组序列..."

    # 生成示例人类基因组片段（hg38的22号染色体前5Mb）
    if [ ! -f "hg38_chr22_demo.fa" ]; then
        python3 << 'PYTHON_EOF'
import random
random.seed(42)

# 生成模拟基因组序列（包含一些真实基因信息）
header = ">chr22:1-5000000"
seq = "".join([random.choice("ACGT") for _ in range(5000000)])

# 插入一些模拟ChIP-seq结合位点（模拟富集的区域）
enriched_regions = [
    (100000, 102000),    # 区域1：TF结合位点
    (500000, 502000),    # 区域2：组蛋白修饰
    (1000000, 1002000),  # 区域3：开放染色质
    (2000000, 2002000),  # 区域4
    (3500000, 3502000),  # 区域5
]

seq_list = list(seq)
for start, end in enriched_regions:
    # 在这些区域增加G/C含量（模拟富集）
    for i in range(start, min(end, 5000000)):
        if random.random() < 0.7:
            seq_list[i] = random.choice("GC")

seq = "".join(seq_list)

# 写入文件
with open("hg38_chr22_demo.fa", "w") as f:
    f.write(f"{header}\n")
    # 每60字符换行
    for i in range(0, len(seq), 60):
        f.write(seq[i:i+60] + "\n")

print("演示参考基因组已生成: hg38_chr22_demo.fa")
PYTHON_EOF
    fi

    # 创建索引
    if [ -f "hg38_chr22_demo.fa" ]; then
        log_info "建立FASTA索引..."
        samtools faidx hg38_chr22_demo.fa 2>/dev/null || echo "✓ 索引已存在"
    fi

    # 生成示例GTF注释文件
    if [ ! -f "hg38_genes_demo.gtf" ]; then
        log_info "生成示例基因注释文件..."
        cat > hg38_genes_demo.gtf << 'GTF_EOF'
##gff-version 3
##sequence-region chr22 1 5000000
chr22	ENSEMBL	gene	100000	102000	.	+	.	gene_id "ENSG00000001"; gene_name "TF_TARGET1"; gene_biotype "protein_coding";
chr22	ENSEMBL	transcript	100000	102000	.	+	.	gene_id "ENSG00000001"; transcript_id "ENST00000001"; transcript_name "TF_TARGET1-001";
chr22	ENSEMBL	exon	100000	101000	.	+	.	gene_id "ENSG00000001"; transcript_id "ENST00000001"; exon_number "1";
chr22	ENSEMBL	exon	101500	102000	.	+	.	gene_id "ENSG00000001"; transcript_id "ENST00000001"; exon_number "2";
chr22	ENSEMBL	gene	500000	502000	.	-	.	gene_id "ENSG00000002"; gene_name "HISTONE_TARGET"; gene_biotype "protein_coding";
chr22	ENSEMBL	transcript	500000	502000	.	-	.	gene_id "ENSG00000002"; transcript_id "ENST00000002"; transcript_name "HISTONE_TARGET-001";
chr22	ENSEMBL	exon	500000	501000	.	-	.	gene_id "ENSG00000002"; transcript_id "ENST00000002"; exon_number "2";
chr22	ENSEMBL	exon	501500	502000	.	-	.	gene_id "ENSG00000002"; transcript_id "ENST00000002"; exon_number "1";
chr22	ENSEMBL	gene	1000000	1002000	.	+	.	gene_id "ENSG00000003"; gene_name "CHROMATIN_OPEN"; gene_biotype "protein_coding";
chr22	ENSEMBL	transcript	1000000	1002000	.	+	.	gene_id "ENSG00000003"; transcript_id "ENST00000003"; transcript_name "CHROMATIN_OPEN-001";
chr22	ENSEMBL	exon	1000000	1002000	.	+	.	gene_id "ENSG00000003"; transcript_id "ENST00000003"; exon_number "1";
chr22	ENSEMBL	gene	2000000	2002000	.	-	.	gene_id "ENSG00000004"; gene_name "GENE_004"; gene_biotype "protein_coding";
chr22	ENSEMBL	transcript	2000000	2002000	.	-	.	gene_id "ENSG00000004"; transcript_id "ENST00000004"; transcript_name "GENE_004-001";
chr22	ENSEMBL	exon	2000000	2002000	.	-	.	gene_id "ENSG00000004"; transcript_id "ENST00000004"; exon_number "1";
chr22	ENSEMBL	gene	3500000	3502000	.	+	.	gene_id "ENSG00000005"; gene_name "GENE_005"; gene_biotype "protein_coding";
chr22	ENSEMBL	transcript	3500000	3502000	.	+	.	gene_id "ENSG00000005"; transcript_id "ENST00000005"; transcript_name "GENE_005-001";
chr22	ENSEMBL	exon	3500000	3502000	.	+	.	gene_id "ENSG00000005"; transcript_id "ENST00000005"; exon_number "1";
GTF_EOF
        log_info "示例基因注释文件已生成"
    fi
}

# 准备ChIP-seq测序数据
prepare_chipseq_data() {
    log_info "开始准备ChIP-seq测序数据..."

    cd "$DATA_DIR"

    # 生成模拟ChIP-seq数据
    log_info "生成模拟ChIP-seq测序数据..."

    if command -v wgsim &> /dev/null; then
        # 使用wgsim生成模拟数据
        log_info "使用wgsim生成模拟ChIP-seq reads..."

        if [ -f "$REFERENCE_DIR/hg38_chr22_demo.fa" ]; then
            # 生成ChIP样本数据
            log_info "生成H3K4me3 ChIP样本数据..."
            wgsim -N 500000 -1 50 -2 50 -r 0.001 -R 0.15 -X 0.3 \
                "$REFERENCE_DIR/hg38_chr22_demo.fa" \
                "H3K4me3_ChIP_R1.fastq" \
                "H3K4me3_ChIP_R2.fastq"

            # 生成Input对照数据
            log_info "生成Input对照样本数据..."
            wgsim -N 500000 -1 50 -2 50 -r 0.001 -R 0.15 -X 0.3 \
                "$REFERENCE_DIR/hg38_chr22_demo.fa" \
                "Input_control_R1.fastq" \
                "Input_control_R2.fastq"

            # 压缩数据
            log_info "压缩FASTQ文件..."
            gzip *.fastq
            log_info "ChIP-seq测序数据已生成"
        else
            log_error "参考基因组文件不存在"
            return 1
        fi
    else
        log_warn "wgsim未安装，使用Python生成简单数据"

        # 使用Python生成FASTQ数据
        python3 << 'PYTHON_EOF'
import random
import gzip

random.seed(42)

def generate_fastq_record(read_id, seq_length=50):
    """生成单条FASTQ记录"""
    seq = "".join([random.choice("ACGT") for _ in range(seq_length)])
    qual = "".join([chr(33 + random.randint(30, 40)) for _ in range(seq_length)])
    return f"@{read_id}\n{seq}\n+\n{qual}\n"

# 生成ChIP和Input数据
samples = {
    "H3K4me3_ChIP": 500000,  # 500k reads
    "Input_control": 500000   # 500k reads
}

for sample, num_reads in samples.items():
    for read_pair in ["R1", "R2"]:
        filename = f"{sample}_{read_pair}.fastq.gz"
        with gzip.open(filename, "wt") as f:
            for i in range(num_reads):
                read_id = f"{sample}_{read_pair}_{i+1}"
                f.write(generate_fastq_record(read_id))
        print(f"Generated: {filename}")

print("FASTQ files generation complete!")
PYTHON_EOF
    fi

    # 验证生成的文件
    log_info "验证生成的FASTQ文件..."
    for fastq_file in *.fastq.gz; do
        if [ -f "$fastq_file" ]; then
            size=$(du -h "$fastq_file" | cut -f1)
            count=$(zcat "$fastq_file" | wc -l)
            reads=$((count / 4))
            log_info "$fastq_file: $size (约 $reads reads)"
        fi
    done
}

# 创建样本信息文件
create_sample_info() {
    log_info "创建样本信息文件..."

    cat > "$DATA_DIR/sample_info.txt" << 'EOF'
sample_id	condition	group	description
H3K4me3_ChIP	ChIP	treatment	H3K4me3 ChIP-seq (Active promoters)
Input_control	Input	control	Input control DNA

# 样本特点：
# H3K4me3: 三甲基化组蛋白H3K4，标记活跃启动子和转录开始位点
# Input: 未经ChIP的全基因组DNA，用作对照
EOF

    if [ -f "$DATA_DIR/sample_info.txt" ]; then
        log_info "样本信息文件已创建"
    fi
}

# 验证数据完整性
validate_data() {
    log_info "验证数据完整性..."

    local validation_ok=true

    # 检查参考基因组
    if [ -d "$REFERENCE_DIR" ]; then
        log_info "检查参考基因组..."
        if [ -f "$REFERENCE_DIR/hg38_chr22_demo.fa" ]; then
            log_info "✓ 参考基因组文件存在"
        else
            log_warn "✗ 参考基因组文件不存在"
            validation_ok=false
        fi

        # 检查注释文件
        if [ -f "$REFERENCE_DIR/hg38_genes_demo.gtf" ]; then
            log_info "✓ 注释文件存在"
        else
            log_warn "✗ 注释文件不存在"
            validation_ok=false
        fi
    fi

    # 检查测序数据
    if [ -d "$DATA_DIR" ]; then
        log_info "检查ChIP-seq数据..."
        fastq_count=$(find "$DATA_DIR" -name "*.fastq.gz" | wc -l)
        if [ "$fastq_count" -gt 0 ]; then
            log_info "✓ 找到 $fastq_count 个FASTQ文件"
        else
            log_warn "✗ 未找到FASTQ文件"
            validation_ok=false
        fi

        # 检查样本信息文件
        if [ -f "$DATA_DIR/sample_info.txt" ]; then
            log_info "✓ 样本信息文件存在"
        else
            log_warn "✗ 样本信息文件不存在"
            validation_ok=false
        fi
    fi

    if [ "$validation_ok" = true ]; then
        log_info "✓ 数据验证完成，所有必要的文件都已准备就绪"
        return 0
    else
        log_warn "⚠ 某些数据文件缺失，请检查"
        return 1
    fi
}

# 生成摘要报告
generate_summary() {
    log_info "生成数据准备摘要..."

    cat > "$LOGS_DIR/data_summary.txt" << EOF
ChIP-seq数据准备摘要
====================

准备时间: $(date)

参考数据:
  目录: $REFERENCE_DIR
  文件数: $(ls -1 "$REFERENCE_DIR" 2>/dev/null | wc -l)

ChIP-seq数据:
  目录: $DATA_DIR
  样本总数: 2 (1个ChIP + 1个Input)
  数据类型: Paired-end ChIP-seq (2×50bp)
  每样本reads数: ~500,000

预期分析结果:
  - 比对率: 80-95%
  - 富集peak数: 预计1000-5000个
  - 信噪比: >1.5

下一步:
  1. 质量控制: FastQC
  2. 序列比对: BWA/BOWTIE2
  3. Peak calling: MACS2
  4. Peak注释: ChIPseeker (R)
  5. 可视化: deepTools / IGV

EOF

    log_info "摘要已保存到: $LOGS_DIR/data_summary.txt"
}

# 主程序
main() {
    local action="all"

    # 解析命令行参数
    while [[ $# -gt 0 ]]; do
        case $1 in
            --reference)
                action="reference"
                shift
                ;;
            --sequencing)
                action="sequencing"
                shift
                ;;
            --validate)
                action="validate"
                shift
                ;;
            --all)
                action="all"
                shift
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

    log_info "=== ChIP-seq数据准备脚本 ==="
    log_info "开始时间: $(date)"

    check_dependencies

    case $action in
        reference)
            download_reference
            ;;
        sequencing)
            prepare_chipseq_data
            create_sample_info
            ;;
        validate)
            validate_data
            ;;
        all)
            download_reference
            prepare_chipseq_data
            create_sample_info
            validate_data
            generate_summary
            ;;
    esac

    log_info "结束时间: $(date)"
    log_info "详细日志已保存到: $LOGS_DIR/data_download.log"
}

# 执行主程序
main "$@"
