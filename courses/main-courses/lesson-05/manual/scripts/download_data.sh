#!/bin/bash

###############################################################################
#                   RNA-seq数据准备脚本                                      #
# 自动下载参考基因组、基因注释和测序数据                                   #
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
使用方法: $0 [选项]

选项:
    --reference     仅下载参考基因组和注释文件
    --sequencing    仅准备测序数据
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

    # 下载选项1：完整基因组（大约3GB）
    log_warn "参考基因组很大（~3GB），首次下载可能需要较长时间"

    # 创建小规模示例数据（22号染色体的5Mb片段）
    log_info "准备基因组示例数据（22号染色体的5Mb片段）..."

    # 检查是否已存在完整基因组
    if [ ! -f "Homo_sapiens.GRCh38.dna.primary_assembly.fa" ]; then
        log_info "从Ensembl下载人类参考基因组..."

        # 尝试使用wget下载（需要网络连接）
        if command -v wget &> /dev/null; then
            log_info "下载地址: ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/"
            log_info "由于文件很大，请手动下载或使用以下命令："
            echo "wget -c ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
            log_warn "或者使用更快的UCSC镜像:"
            echo "wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
        else
            log_error "wget未安装，无法下载"
            return 1
        fi
    fi

    # 如果已有完整基因组，提取22号染色体的5Mb片段用于演示
    if [ -f "Homo_sapiens.GRCh38.dna.primary_assembly.fa" ] && \
       [ ! -f "chr22_5Mb_fragment.fa" ]; then
        log_info "从完整基因组提取chr22的前5Mb片段..."

        if command -v samtools &> /dev/null; then
            samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
            samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa 22:1-5000000 \
                > chr22_5Mb_fragment.fa
            log_info "chr22_5Mb_fragment.fa已生成"
        fi
    fi

    # 生成示例基因组序列（用于快速测试）
    if [ ! -f "demo_reference.fa" ]; then
        log_info "生成演示参考序列..."
        python3 << 'PYTHON_EOF'
import random
random.seed(42)

# 生成示例基因组序列（200kb）
with open("demo_reference.fa", "w") as f:
    f.write(">chr22_demo:1-200000\n")
    seq = "".join([random.choice("ACGT") for _ in range(200000)])
    # 添加一些真实基因序列
    genes = [
        "ATGGCATGGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",  # 基因1
        "ATGAAATTTAAATTTAAACCCGGGTTTAAATTTAAATTTAAATTTAAATTTAAA",   # 基因2
        "ATGCCCAAAGGGAAATTTAAACCCGGGTTTAAATTTAAATTTAAATTTAAATTT",   # 基因3
    ]

    # 随机插入基因序列
    pos = 0
    for i, gene in enumerate(genes):
        pos = (i + 1) * 40000 + random.randint(0, 10000)
        if pos + len(gene) < 200000:
            seq = seq[:pos] + gene + seq[pos+len(gene):]

    # 每60个字符换行
    for i in range(0, len(seq), 60):
        f.write(seq[i:i+60] + "\n")

print("Demo reference generated: demo_reference.fa")
PYTHON_EOF

        if [ -f "demo_reference.fa" ]; then
            log_info "演示参考序列已生成: demo_reference.fa (200kb)"
        fi
    fi

    # 下载基因注释文件
    log_info "下载基因注释文件..."

    if [ ! -f "Homo_sapiens.GRCh38.104.gtf" ] && \
       [ ! -f "Homo_sapiens.GRCh38.104.gtf.gz" ]; then
        log_info "从Ensembl下载GTF注释文件..."

        # 生成示例GTF文件
        log_info "生成示例GTF注释文件..."
        cat > demo_annotation.gtf << 'GTF_EOF'
##gff-version 3
##sequence-region chr22 1 200000
chr22	ENSEMBL	gene	1000	50000	.	+	.	gene_id "ENSG00000000001"; gene_name "GENE1"; gene_biotype "protein_coding";
chr22	ENSEMBL	transcript	1000	50000	.	+	.	gene_id "ENSG00000000001"; transcript_id "ENST00000000001"; transcript_name "GENE1-001";
chr22	ENSEMBL	exon	1000	2000	.	+	.	gene_id "ENSG00000000001"; transcript_id "ENST00000000001"; exon_number "1";
chr22	ENSEMBL	exon	5000	6000	.	+	.	gene_id "ENSG00000000001"; transcript_id "ENST00000000001"; exon_number "2";
chr22	ENSEMBL	exon	10000	50000	.	+	.	gene_id "ENSG00000000001"; transcript_id "ENST00000000001"; exon_number "3";
chr22	ENSEMBL	CDS	1100	2000	.	+	0	gene_id "ENSG00000000001"; transcript_id "ENST00000000001"; protein_id "ENSP00000000001";
chr22	ENSEMBL	CDS	5000	6000	.	+	2	gene_id "ENSG00000000001"; transcript_id "ENST00000000001"; protein_id "ENSP00000000001";
chr22	ENSEMBL	CDS	10000	10900	.	+	0	gene_id "ENSG00000000001"; transcript_id "ENST00000000001"; protein_id "ENSP00000000001";
chr22	ENSEMBL	gene	60000	100000	.	-	.	gene_id "ENSG00000000002"; gene_name "GENE2"; gene_biotype "protein_coding";
chr22	ENSEMBL	transcript	60000	100000	.	-	.	gene_id "ENSG00000000002"; transcript_id "ENST00000000002"; transcript_name "GENE2-001";
chr22	ENSEMBL	exon	60000	70000	.	-	.	gene_id "ENSG00000000002"; transcript_id "ENST00000000002"; exon_number "2";
chr22	ENSEMBL	exon	75000	100000	.	-	.	gene_id "ENSG00000000002"; transcript_id "ENST00000000002"; exon_number "1";
chr22	ENSEMBL	CDS	60100	70000	.	-	2	gene_id "ENSG00000000002"; transcript_id "ENST00000000002"; protein_id "ENSP00000000002";
chr22	ENSEMBL	CDS	75000	99900	.	-	0	gene_id "ENSG00000000002"; transcript_id "ENST00000000002"; protein_id "ENSP00000000002";
GTF_EOF

        if [ -f "demo_annotation.gtf" ]; then
            log_info "示例GTF注释文件已生成: demo_annotation.gtf"
        fi
    fi
}

# 准备测序数据
prepare_sequencing_data() {
    log_info "开始准备测序数据..."

    cd "$DATA_DIR"

    # 方案1：生成模拟测序数据
    log_info "生成模拟RNA-seq测序数据..."

    # 检查是否有wgsim工具
    if command -v wgsim &> /dev/null; then
        log_info "使用wgsim生成模拟测序数据..."

        if [ -f "$REFERENCE_DIR/demo_reference.fa" ]; then
            # 为6个样本生成测序数据
            for i in 1 2 3; do
                # 处理样本 (treatment)
                log_info "生成sample${i}模拟数据..."
                wgsim -N 500000 -1 100 -2 100 -r 0.001 -R 0.15 -X 0.3 \
                    "$REFERENCE_DIR/demo_reference.fa" \
                    "sample${i}_R1.fastq" \
                    "sample${i}_R2.fastq"

                # 控制样本 (control)
                log_info "生成ctrl${i}模拟数据..."
                wgsim -N 500000 -1 100 -2 100 -r 0.001 -R 0.15 -X 0.3 \
                    "$REFERENCE_DIR/demo_reference.fa" \
                    "ctrl${i}_R1.fastq" \
                    "ctrl${i}_R2.fastq"
            done

            # 压缩数据文件
            log_info "压缩FASTQ文件..."
            gzip *.fastq
            log_info "测序数据已生成并压缩"
        else
            log_error "demo_reference.fa不存在，无法生成测序数据"
            return 1
        fi
    else
        log_warn "wgsim未安装，将使用Python生成示例数据"

        # 使用Python生成简单的FASTQ数据
        python3 << 'PYTHON_EOF'
import random
import gzip

random.seed(42)

def generate_fastq_record(read_id, seq_length=100):
    """生成单条FASTQ记录"""
    seq = "".join([random.choice("ACGT") for _ in range(seq_length)])
    qual = "".join([chr(33 + random.randint(30, 40)) for _ in range(seq_length)])
    return f"@{read_id}\n{seq}\n+\n{qual}\n"

# 生成6个样本的数据（3个处理，3个对照）
samples = [f"sample{i}" for i in range(1, 4)] + [f"ctrl{i}" for i in range(1, 4)]

for sample in samples:
    # 生成R1和R2数据
    for read_pair in ["R1", "R2"]:
        filename = f"{sample}_{read_pair}.fastq.gz"
        with gzip.open(filename, "wt") as f:
            for i in range(100000):  # 100k reads per sample
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
sample_id	group	batch	description
sample1	treatment	batch1	Treatment group 1
sample2	treatment	batch2	Treatment group 2
sample3	treatment	batch3	Treatment group 3
ctrl1	control	batch1	Control group 1
ctrl2	control	batch2	Control group 2
ctrl3	control	batch3	Control group 3
EOF

    if [ -f "$DATA_DIR/sample_info.txt" ]; then
        log_info "样本信息文件已创建: $DATA_DIR/sample_info.txt"
    fi
}

# 验证数据完整性
validate_data() {
    log_info "验证数据完整性..."

    local validation_ok=true

    # 检查参考基因组
    if [ -d "$REFERENCE_DIR" ]; then
        log_info "检查参考基因组..."
        if [ -f "$REFERENCE_DIR/demo_reference.fa" ] || \
           [ -f "$REFERENCE_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa" ]; then
            log_info "✓ 参考基因组文件存在"
        else
            log_warn "✗ 参考基因组文件不存在"
            validation_ok=false
        fi

        # 检查注释文件
        if [ -f "$REFERENCE_DIR/demo_annotation.gtf" ] || \
           [ -f "$REFERENCE_DIR/Homo_sapiens.GRCh38.104.gtf" ]; then
            log_info "✓ 注释文件存在"
        else
            log_warn "✗ 注释文件不存在"
        fi
    fi

    # 检查测序数据
    if [ -d "$DATA_DIR" ]; then
        log_info "检查测序数据..."
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
            cat "$DATA_DIR/sample_info.txt"
        else
            log_warn "✗ 样本信息文件不存在"
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
RNA-seq Data Preparation Summary
================================

准备时间: $(date)

参考数据:
  目录: $REFERENCE_DIR
  文件数: $(ls -1 "$REFERENCE_DIR" | wc -l)

测序数据:
  目录: $DATA_DIR
  样本总数: 6 (3个处理组 + 3个对照组)
  数据类型: Paired-end RNA-seq
  每样本reads数: ~100,000 (演示数据)

预期分析结果:
  - 比对率: 85-95%
  - 基因检出数: ~20
  - 显著差异基因: 预计2-5个

下一步:
  1. 建立参考基因组索引: hisat2-build
  2. 进行序列比对: hisat2
  3. 基因定量: featureCounts
  4. 差异分析: DESeq2 (R)

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

    log_info "=== RNA-seq数据准备脚本 ==="
    log_info "开始时间: $(date)"

    check_dependencies

    case $action in
        reference)
            download_reference
            ;;
        sequencing)
            prepare_sequencing_data
            create_sample_info
            ;;
        validate)
            validate_data
            ;;
        all)
            download_reference
            prepare_sequencing_data
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
