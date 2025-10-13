# 第4次课数据来源说明

## 概述

本课程实验使用经过优化的小规模数据集,旨在让学生在有限的计算资源和时间内完成完整的变异检测流程。

## 数据文件详情

### 1. 参考基因组

**文件名:** `chr22_fragment.fa`
**大小:** ~5 MB (约5,000,000 bp)
**来源:**
- 主要来源: [Ensembl](http://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/)
- 备用来源: [UCSC Genome Browser](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/)

**描述:**
- 人类基因组GRCh38/hg38版本
- 22号染色体的10,000,000-15,000,000区域片段
- 选择该区域的原因:
  - 22号染色体是最小的常染色体
  - 该区域包含多个基因,有足够的变异位点
  - 大小适中,适合教学演示

**获取方式:**
```bash
# 自动下载(推荐)
bash scripts/download_data.sh

# 手动下载
wget -O chr22.fa.gz \
  "http://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz"
gunzip chr22.fa.gz
samtools faidx chr22.fa
samtools faidx chr22.fa 22:10000000-15000000 > chr22_fragment.fa
```

### 2. 测序数据

**文件名:** `sample.bam`
**大小:** ~50-100 MB
**来源:**
- 优先: 使用lesson-03生成的比对数据
- 备选: 模拟生成的双端测序数据

**数据特征:**
- 测序平台: 模拟Illumina HiSeq
- 读长: 2×150 bp (双端测序)
- Reads数量: ~100,000 pairs (~200,000 reads)
- 覆盖深度: ~30X
- 质量分数范围: Phred 30-40
- 插入片段大小: ~350 bp (±50 bp)

**数据生成方式:**
如果从头生成模拟数据,使用以下步骤:
1. 从参考基因组随机采样序列
2. 引入随机测序错误 (错误率 ~0.1%)
3. 使用BWA进行比对
4. 标记重复序列和排序

**质量特点:**
- 比对率: ~95%
- 正确配对率: ~90%
- 重复率: ~2-5%

### 3. 已知变异数据库 (dbSNP)

**文件名:** `dbsnp.vcf`
**大小:** ~50 KB
**来源:** 模拟生成(教学用)

**注意事项:**
⚠️ **这是教学用的简化版本,真实项目应使用官方数据库!**

**真实数据来源:**
- 官方: [NCBI dbSNP](https://www.ncbi.nlm.nih.gov/snp/)
- GATK资源包: [GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811)

**模拟数据特征:**
- 变异数量: ~1,000个SNP位点
- 覆盖区域: chr22:10000000-15000000
- 用途: 用于BQSR校正和变异质量评估

**获取真实dbSNP数据:**
```bash
# 从GATK Resource Bundle下载 (需要Google Cloud)
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz .

# 提取22号染色体区域
bcftools view -r chr22:10000000-15000000 \
  Homo_sapiens_assembly38.dbsnp138.vcf.gz \
  > dbsnp_chr22_fragment.vcf
```

### 4. HapMap高质量SNP

**文件名:** `hapmap.vcf`
**大小:** ~25 KB
**来源:** 模拟生成(教学用)

**真实数据来源:**
- HapMap项目: [International HapMap Project](https://www.genome.gov/10001688/international-hapmap-project)
- GATK资源包: [GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811)

**模拟数据特征:**
- 变异数量: ~500个高质量SNP
- 等位基因频率: 0.1-0.5
- 用途: 用于VQSR训练(如果样本量足够)

## 数据准备工作流程

### 自动化准备(推荐)

```bash
# 在课程目录下运行
cd ~/ngs-analysis/lesson-04
bash scripts/download_data.sh
```

该脚本将自动:
1. ✓ 检查必要的软件环境
2. ✓ 下载参考基因组并提取片段
3. ✓ 生成或复用测序数据
4. ✓ 进行序列比对和排序
5. ✓ 创建资源文件
6. ✓ 建立所有必要的索引

### 手动准备步骤

如果自动脚本无法运行,可以按以下步骤手动准备:

1. **准备工作目录**
```bash
mkdir -p ~/ngs-analysis/lesson-04/data
cd ~/ngs-analysis/lesson-04/data
```

2. **下载参考基因组**
```bash
# 见"参考基因组"部分
```

3. **准备测序数据**
```bash
# 如果有lesson-03的数据
cp ~/ngs-analysis/lesson-03/results/bwa_default_sorted.bam sample.bam
samtools index sample.bam

# 或生成新的模拟数据(见download_data.sh脚本)
```

4. **创建资源文件**
```bash
# 运行download_data.sh中的Python脚本生成VCF文件
```

## 数据验证

运行以下命令验证数据完整性:

```bash
cd ~/ngs-analysis/lesson-04/data

# 检查参考基因组
echo "参考基因组:"
ls -lh chr22_fragment.fa*
head -2 chr22_fragment.fa

# 检查BAM文件
echo -e "\nBAM文件:"
ls -lh sample.bam*
samtools flagstat sample.bam | head -5

# 检查资源文件
echo -e "\n资源文件:"
ls -lh *.vcf
echo "dbSNP变异数: $(grep -v '^#' dbsnp.vcf | wc -l)"
echo "HapMap变异数: $(grep -v '^#' hapmap.vcf | wc -l)"
```

**预期输出:**
```
参考基因组:
-rw-r--r-- 1 user user 5.0M chr22_fragment.fa
-rw-r--r-- 1 user user  157 chr22_fragment.fa.fai

BAM文件:
-rw-r--r-- 1 user user  80M sample.bam
-rw-r--r-- 1 user user 128K sample.bam.bai
100000 + 0 in total (QC-passed reads + QC-failed reads)
95000 + 0 properly paired

资源文件:
-rw-r--r-- 1 user user  50K dbsnp.vcf
-rw-r--r-- 1 user user  25K hapmap.vcf
dbSNP变异数: 1000
HapMap变异数: 500
```

## 数据使用注意事项

### 教学使用
✅ **适合用于:**
- 学习GATK变异检测流程
- 理解VCF格式和质量控制
- 练习变异注释和解读
- 掌握生物信息学工作流程

⚠️ **不适合用于:**
- 发表研究论文
- 临床诊断决策
- 真实样本分析比较

### 与真实数据的差异

| 特征 | 教学数据 | 真实全基因组数据 |
|-----|---------|----------------|
| 基因组大小 | 5 Mb | 3,000 Mb |
| Reads数量 | 10万 | 3-10亿 |
| 文件大小 | ~100 MB | ~100-300 GB |
| 变异数量 | 500-1000 | 400-500万 |
| 处理时间 | 10-30分钟 | 数小时到数天 |
| 内存需求 | 4-8 GB | 32-128 GB |
| dbSNP重叠率 | 30-50% | >95% |

### 扩展到真实项目

如果要分析真实样本,需要:

1. **使用完整参考基因组**
   - GRCh38完整版本 (~3 GB)
   - 包含所有染色体和未定位序列

2. **下载官方资源文件**
   - GATK Resource Bundle (dbSNP, 1000G, HapMap, Mills等)
   - 大小: 数GB到数十GB

3. **准备充足的计算资源**
   - 内存: 至少32 GB, 推荐64 GB+
   - 存储: 500 GB+ 可用空间
   - CPU: 16核+ 推荐

4. **增加处理时间预期**
   - 变异检测: 数小时
   - VQSR过滤: 需要大样本量(通常>30个样本)
   - 注释: 可能需要数小时

## 故障排除

### 下载失败
**问题:** 无法从Ensembl/UCSC下载数据
**解决:**
- 检查网络连接
- 尝试备用数据源
- 使用提供的本地数据(如果可用)

### 空间不足
**问题:** 磁盘空间不足
**解决:**
```bash
# 检查可用空间
df -h .

# 清理临时文件
rm -rf ~/ngs-analysis/lesson-*/temp/*
```

### 文件损坏
**问题:** 下载的文件损坏或不完整
**解决:**
```bash
# 删除并重新下载
rm chr22_fragment.fa
bash scripts/download_data.sh
```

## 参考资源

### 数据库
- [Ensembl](https://www.ensembl.org/) - 基因组数据库
- [UCSC Genome Browser](https://genome.ucsc.edu/) - 基因组浏览器
- [NCBI dbSNP](https://www.ncbi.nlm.nih.gov/snp/) - 变异数据库
- [1000 Genomes Project](https://www.internationalgenome.org/) - 人群变异数据

### GATK资源
- [GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811)
- [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651)

### 工具文档
- [samtools](http://www.htslib.org/doc/samtools.html)
- [bcftools](http://samtools.github.io/bcftools/)
- [IGV](https://software.broadinstitute.org/software/igv/)

---

**最后更新:** 2025年1月
**维护者:** 王运生 (wangys@hunau.edu.cn)
