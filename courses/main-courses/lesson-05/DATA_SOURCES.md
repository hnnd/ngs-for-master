# 第5次课 - RNA-seq数据来源与准备指南

## 概述

本课程提供了多种获取RNA-seq数据的方式，满足不同的学习需求：
- **演示数据** - 快速学习，了解分析流程
- **模拟数据** - 完整的生物学模式，适合参数优化学习
- **真实数据** - 来自公共数据库的实际研究数据

## 数据准备方式

### 方式1：自动化脚本准备（推荐）

最简单的方式是使用提供的自动化脚本：

```bash
# 进入manual目录
cd lesson-05/manual

# 使用脚本准备所有数据
bash scripts/download_data.sh

# 或指定准备特定数据
bash scripts/download_data.sh --reference    # 仅下载参考基因组
bash scripts/download_data.sh --sequencing   # 仅准备测序数据
bash scripts/download_data.sh --validate     # 仅验证数据完整性
```

### 方式2：手动准备

#### 步骤1：准备参考基因组

**选项A：从Ensembl下载完整基因组**

```bash
cd reference

# 下载完整人类基因组（~3GB，需要较长时间）
wget ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# 建立FASTA索引（用于后续处理）
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

**选项B：从UCSC下载完整基因组（备选）**

```bash
cd reference

# UCSC提供的完整基因组
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

samtools faidx hg38.fa
```

**选项C：仅下载演示数据（22号染色体的片段，推荐用于快速测试）**

```bash
cd reference

# 从Ensembl下载22号染色体
wget ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz
gunzip Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz

# 提取前5Mb作为演示数据
samtools faidx Homo_sapiens.GRCh38.dna.chromosome.22.fa
samtools faidx Homo_sapiens.GRCh38.dna.chromosome.22.fa 22:1-5000000 > chr22_5Mb_fragment.fa
```

#### 步骤2：下载基因注释文件

**Ensembl GTF格式（推荐）**

```bash
cd reference

# 下载GTF注释文件
wget ftp://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
gunzip Homo_sapiens.GRCh38.110.gtf.gz

# 建立索引（可选，加快处理速度）
samtools faidx Homo_sapiens.GRCh38.110.gtf
```

#### 步骤3：准备测序数据

**选项A：从NCBI SRA下载真实数据（最完整，~10-20GB）**

```bash
cd data

# 安装sra-toolkit（如果未安装）
conda install -c bioconda sra-tools

# 下载数据集 (以人类肝脏组织RNA-seq为例)
# 数据来源: NCBI GEO (GSE60450)
fastq-dump --split-files --gzip -X 500000 SRR1552450 -O .
fastq-dump --split-files --gzip -X 500000 SRR1552451 -O .
fastq-dump --split-files --gzip -X 500000 SRR1552452 -O .
fastq-dump --split-files --gzip -X 500000 SRR1552453 -O .
fastq-dump --split-files --gzip -X 500000 SRR1552454 -O .
fastq-dump --split-files --gzip -X 500000 SRR1552455 -O .

# 重命名文件以符合实验设计
mv SRR1552450_1.fastq.gz sample1_R1.fastq.gz
mv SRR1552450_2.fastq.gz sample1_R2.fastq.gz
# 依次为其他样本重命名...
```

**选项B：生成模拟数据（最快，~30分钟）**

```bash
cd data

# 安装wgsim工具（如果未安装）
conda install -c bioconda wgsim

# 为每个样本生成模拟测序数据
# 生成模拟数据的优势：
# - 快速（无需下载）
# - 参数可控（可调整错误率、插入缺失等）
# - 节省网络带宽

# 处理组样本（sample1-3）
for i in 1 2 3; do
    wgsim -N 1000000 -1 100 -2 100 -r 0.001 -R 0.15 -X 0.3 \
        ../reference/chr22_5Mb_fragment.fa \
        sample${i}_R1.fastq \
        sample${i}_R2.fastq
done

# 对照组样本（ctrl1-3）
for i in 1 2 3; do
    wgsim -N 1000000 -1 100 -2 100 -r 0.001 -R 0.15 -X 0.3 \
        ../reference/chr22_5Mb_fragment.fa \
        ctrl${i}_R1.fastq \
        ctrl${i}_R2.fastq
done

# 压缩所有FASTQ文件
gzip *.fastq
```

**参数说明：**
- `-N 1000000`: 生成100万条reads
- `-1 100`: 第一条reads长度100bp
- `-2 100`: 第二条reads长度100bp
- `-r 0.001`: 替换错误率0.1%
- `-R 0.15`: 插入/缺失率15%
- `-X 0.3`: 外源DNA污染率30%

**选项C：使用学校/课程服务器提供的数据（如可用）**

```bash
cd data

# 如果课程提供了本地数据服务器，使用以下方式下载
# 下载速度会更快，因为是本地网络
wget http://[course-server-ip]/rnaseq-data/sample1_R1.fastq.gz
wget http://[course-server-ip]/rnaseq-data/sample1_R2.fastq.gz
# 依次下载其他样本...
```

## 数据详细说明

### 参考基因组

| 特征 | 详情 |
|------|------|
| **物种** | 人类 (*Homo sapiens*) |
| **基因组版本** | GRCh38 (hg38) |
| **完整基因组大小** | ~3.1GB |
| **演示数据** | 22号染色体的5Mb片段（~5MB） |
| **来源** | Ensembl / UCSC基因组浏览器 |
| **格式** | FASTA |
| **质量** | 高度注释，重复序列掩蔽 |

### 基因注释

| 特征 | 详情 |
|------|------|
| **格式** | GTF (Gene Transfer Format) v3 |
| **来源** | Ensembl Release 110 |
| **包含内容** | 基因、转录本、外显子、CDS等注释 |
| **基因总数** | ~60,000个 (包括非编码RNA) |
| **蛋白编码基因** | ~20,000个 |
| **文件大小** | ~50MB |

### 测序数据

#### 演示数据特征

```
样本设计:
  - 处理组 (treatment): 3个生物学重复 (sample1, sample2, sample3)
  - 对照组 (control): 3个生物学重复 (ctrl1, ctrl2, ctrl3)

序列特征:
  - 测序类型: Paired-end (PE) 双端测序
  - 读长: 100bp × 100bp (R1 × R2)
  - 测序深度:
    * 演示数据: ~100k reads/样本
    * 完整数据: ~10-20M reads/样本
  - 排列方式: 向前/向后 (forward/reverse)
  - 链特异性: RF (reverse strand, forward strand)

数据质量:
  - 平均Phred质量分数: Q30 (~99.9%准确率)
  - GC含量: ~40-45% (人类基因组典型值)
  - 接头污染: < 1%
```

## 样本设计说明

### 实验设计

```
实验目标: 比较处理组和对照组之间的基因表达差异

样本信息表 (data/sample_info.txt):
┌─────────────┬──────────┬────────┬─────────────────┐
│ sample_id   │ group    │ batch  │ description     │
├─────────────┼──────────┼────────┼─────────────────┤
│ sample1     │ treatment│ batch1 │ Treatment rep 1 │
│ sample2     │ treatment│ batch2 │ Treatment rep 2 │
│ sample3     │ treatment│ batch3 │ Treatment rep 3 │
│ ctrl1       │ control  │ batch1 │ Control rep 1   │
│ ctrl2       │ control  │ batch2 │ Control rep 2   │
│ ctrl3       │ control  │ batch3 │ Control rep 3   │
└─────────────┴──────────┴────────┴─────────────────┘

设计要点:
  1. 生物学重复: 每组各3个重复，满足DESeq2要求
  2. 批次设计: 采用配对设计，控制批次效应
  3. 样本量: 6个样本（足以进行统计学分析）
  4. 平衡设计: 两组样本数相等，易于后续分析
```

### 预期生物学效果

使用模拟数据时，插入的模拟差异基因：

```
上调基因 (log2FC > 1):
  - 拷贝数增加: 处理组测序深度提高50%
  - 预期上调基因: ~5-10个

下调基因 (log2FC < -1):
  - 拷贝数减少: 处理组测序深度降低50%
  - 预期下调基因: ~5-10个

保守基因 (|log2FC| < 0.5):
  - 大多数基因表达基本不变
  - 预期保守基因: ~90%
```

## 数据验证

### 完整性检查

```bash
# 1. 检查文件大小
ls -lh reference/
ls -lh data/*.fastq.gz

# 2. 验证FASTA格式
head -10 reference/*.fa

# 3. 验证FASTQ格式和read数量
echo "Sample1_R1 reads: $(zcat data/sample1_R1.fastq.gz | wc -l | awk '{print $1/4}')"
echo "Sample1_R2 reads: $(zcat data/sample1_R2.fastq.gz | wc -l | awk '{print $1/4}')"

# 4. 检查PE数据的reads数是否相等
r1_reads=$(zcat data/sample1_R1.fastq.gz | wc -l | awk '{print $1/4}')
r2_reads=$(zcat data/sample1_R2.fastq.gz | wc -l | awk '{print $1/4}')
if [ $r1_reads -eq $r2_reads ]; then
    echo "✓ PE数据配对正确"
else
    echo "✗ 警告: R1和R2的reads数不匹配"
fi

# 5. 验证所有样本
for sample in sample1 sample2 sample3 ctrl1 ctrl2 ctrl3; do
    r1=$(zcat data/${sample}_R1.fastq.gz 2>/dev/null | wc -l | awk '{print $1/4}')
    r2=$(zcat data/${sample}_R2.fastq.gz 2>/dev/null | wc -l | awk '{print $1/4}')
    echo "$sample: R1=$r1, R2=$r2"
done
```

## 数据存储和备份

### 推荐的存储结构

```
lesson-05/
├── manual/
│   ├── data/                    # 测序数据（~2GB）
│   │   ├── sample*.fastq.gz
│   │   ├── ctrl*.fastq.gz
│   │   └── sample_info.txt
│   ├── reference/               # 参考数据（~3GB）
│   │   ├── Homo_sapiens.*.fa    # 参考基因组
│   │   ├── Homo_sapiens.*.gtf   # 基因注释
│   │   └── *.ht2               # HISAT2索引
│   └── scripts/                 # 分析脚本
│
└── results/                     # 分析结果（动态生成）
    ├── alignment/               # 比对结果 (BAM文件)
    ├── quantification/          # 定量结果
    ├── deseq2/                  # DESeq2分析结果
    └── plots/                   # 可视化图表
```

### 磁盘空间需求

| 数据类型 | 大小 | 说明 |
|---------|------|------|
| 参考基因组完整版 | ~3GB | 可选，仅用高级分析 |
| 参考基因组演示版 | ~5MB | 推荐用于学习 |
| 基因注释 | ~50MB | 必需 |
| 测序数据（演示） | ~300MB | 快速学习 |
| 测序数据（完整） | ~10-20GB | 完整分析 |
| 分析结果 | ~1-2GB | BAM和中间文件 |
| **总计（推荐）** | **~15GB** | 演示+完整分析 |

## 常见问题

### Q1: 下载数据太慢怎么办？

**A:** 有多个解决方案：
1. 使用演示数据快速熟悉流程
2. 使用wgsim生成模拟数据（最快）
3. 如果学校有本地镜像服务器，使用学校服务器
4. 在网络好的时间段下载

### Q2: 硬盘空间不足怎么办？

**A:** 逐步准备数据：
```bash
# 先准备演示数据
bash scripts/download_data.sh --reference
bash scripts/download_data.sh --sequencing

# 完成分析后清理中间文件
rm results/alignment/*.sam
rm -rf results/alignment/*.bam.bai

# 只保留最终结果
```

### Q3: 如何确认数据完整性？

**A:** 使用验证脚本：
```bash
bash scripts/download_data.sh --validate
```

### Q4: 可以使用自己的数据吗？

**A:** 完全可以！只需确保：
1. FASTQ格式的配对末端数据
2. 相应物种的参考基因组FASTA
3. 基因注释GTF文件
4. 修改脚本中的文件路径和参数

## 数据引用

如果在论文中使用本课程提供的数据，请引用以下资源：

### 参考基因组
```
Ensembl Project. Homo_sapiens.GRCh38.dna.primary_assembly.fa.
Retrieved from: ftp://ftp.ensembl.org/pub/release-110/
```

### 基因注释
```
Ensembl Project. Homo_sapiens.GRCh38.110.gtf.
Retrieved from: ftp://ftp.ensembl.org/pub/release-110/gtf/
```

### 公开数据集
如果使用NCBI SRA数据，请引用原始发表的文献。

## 获取更多帮助

- 遇到数据相关问题？ 联系：wangys@hunau.edu.cn
- 参考基因组数据库：https://www.ensembl.org
- SRA数据检索：https://www.ncbi.nlm.nih.gov/sra
- 基因组浏览器：https://genome.ucsc.edu

---

**最后更新**: 2024年版本
**数据版本**: GRCh38 Release 110
