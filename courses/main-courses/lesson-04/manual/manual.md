# 第4次课实践操作手册

## 课程信息
- **课程名称**：高通量测序数据分析
- **主讲教师**：王运生
- **联系邮箱**：wangys@hunau.edu.cn
- **办公室**：16教420室
- **上课地点**：105机房
- **课程时间**：第4次课（2学时实践）

## 实验目标

### 主要目标
- 掌握GATK变异检测的完整流程
- 学会VCF文件的处理和质量评估
- 了解变异注释和功能预测方法
- 理解变异质量控制的重要性

### 预期成果
- 生成高质量的变异检测结果
- 掌握VCF格式文件的分析技能
- 能够进行变异功能注释和解读
- 建立变异检测质量控制体系

## 环境要求

### 软件环境
| 软件名称 | 版本要求 | 安装方式 | 说明 |
|---------|---------|---------|------|
| GATK | 4.2+ | conda install gatk4 | 变异检测主要工具 |
| samtools | 1.15+ | conda install samtools | BAM文件处理 |
| bcftools | 1.15+ | conda install bcftools | VCF文件处理 |
| VEP | 104+ | conda install ensembl-vep | 变异注释工具 |
| IGV | 2.12+ | 官网下载 | 变异可视化 |
| Python | 3.8+ | 系统自带 | 脚本运行环境 |

### 硬件要求
- **内存**：至少 8 GB RAM
- **存储空间**：至少 20 GB 可用空间
- **CPU**：4核心以上推荐
- **网络**：稳定的互联网连接（用于下载参考数据）

### 数据准备
本次实验数据通过自动化脚本准备，包括：

| 数据文件 | 大小 | 来源 | 说明 |
|---------|------|------|------|
| chr22_fragment.fa | ~5MB | Ensembl/UCSC | 人类22号染色体5Mb片段 |
| sample.bam | ~50-100MB | 模拟生成或lesson-03 | 已比对的测序数据(~10万reads) |
| dbsnp.vcf | ~50KB | 模拟生成 | 已知变异数据库(~1000个位点) |
| hapmap.vcf | ~25KB | 模拟生成 | 高质量SNP集合(~500个位点) |

**数据获取方式：**
- 运行 `bash scripts/download_data.sh` 自动下载/生成所有数据
- 数据将保存在 `~/ngs-analysis/lesson-04/data/` 目录
- 参考[DATA_SOURCES.md](../DATA_SOURCES.md)了解数据详细说明

## 操作步骤

### 步骤1：环境设置和准备工作

#### 1.1 创建工作目录
```bash
# 创建本次实验的工作目录
mkdir -p ~/ngs-analysis/lesson-04
cd ~/ngs-analysis/lesson-04

# 创建子目录结构
mkdir -p {data,scripts,results,logs,plots}
```

#### 1.2 检查软件环境
```bash
# 检查GATK安装
gatk --version

# 检查samtools安装
samtools --version

# 检查bcftools安装
bcftools --version

# 检查VEP安装
vep --help | head -5
```

**预期输出：**
```
The Genome Analysis Toolkit (GATK) v4.2.6.1
samtools 1.15.1
bcftools 1.15.1
VEP version 104
```

#### 1.3 下载和准备数据

**自动化数据准备（推荐）：**
```bash
# 运行数据准备脚本
bash scripts/download_data.sh

# 该脚本将自动完成:
# 1. 下载人类22号染色体片段作为参考基因组
# 2. 生成模拟测序数据或使用lesson-03的结果
# 3. 创建变异资源文件(dbSNP和HapMap)
# 4. 建立必要的索引文件
```

**手动数据链接（如果数据已准备好）：**
```bash
# 如果已有准备好的数据，创建符号链接
ln -s ~/ngs-analysis/lesson-04/data/* data/

# 验证数据完整性
ls -lh data/
```

**预期文件列表：**
```
chr22_fragment.fa       参考基因组
chr22_fragment.fa.fai   参考基因组索引
sample.bam              比对数据
sample.bam.bai          BAM索引
dbsnp.vcf               dbSNP数据库
hapmap.vcf              HapMap数据
```

**检查点：** 运行 `ls -lh ~/ngs-analysis/lesson-04/data/` 确认所有数据文件已生成。

---

### 步骤2：数据预处理

#### 2.1 标记重复序列

**操作说明：**
使用GATK MarkDuplicates工具标记PCR重复序列。这些重复序列可能影响变异检测的准确性。

**执行命令：**
```bash
# 标记重复序列
gatk MarkDuplicates \
    -I ~/ngs-analysis/lesson-04/data/sample.bam \
    -O results/sample.marked.bam \
    -M results/duplicate_metrics.txt \
    --CREATE_INDEX true

# 查看重复序列统计
cat results/duplicate_metrics.txt | head -10
```

**参数解释：**
- `-I`：输入BAM文件
- `-O`：输出标记重复后的BAM文件
- `-M`：重复序列统计报告
- `--CREATE_INDEX`：自动创建索引文件

**预期输出示例：**
```
LIBRARY	UNPAIRED_READS_EXAMINED	READ_PAIRS_EXAMINED	UNMAPPED_READS	PERCENT_DUPLICATION
sample	0	100000	5234	0.02-0.05
```

**说明：**
- READ_PAIRS_EXAMINED: ~100,000 (与生成的数据量一致)
- UNMAPPED_READS: ~5% (正常范围)
- PERCENT_DUPLICATION: 2-5% (低重复率,因为是模拟数据)

**结果验证：**
```bash
# 验证输出文件
samtools flagstat results/sample.marked.bam
```

#### 2.2 碱基质量分数校正 (BQSR)

**操作说明：**
BQSR通过已知变异位点校正系统性的测序错误，提高变异检测的准确性。

**执行命令：**
```bash
# 第一步：生成校正表
gatk BaseRecalibrator \
    -I results/sample.marked.bam \
    -R ~/ngs-analysis/lesson-04/data/chr22_fragment.fa \
    --known-sites ~/ngs-analysis/lesson-04/data/dbsnp.vcf \
    -O results/recal_data.table

# 第二步：应用校正
gatk ApplyBQSR \
    -I results/sample.marked.bam \
    -R ~/ngs-analysis/lesson-04/data/chr22_fragment.fa \
    --bqsr-recal-file results/recal_data.table \
    -O results/sample.recal.bam

# 生成校正后的统计表（可选）
gatk BaseRecalibrator \
    -I results/sample.recal.bam \
    -R ~/ngs-analysis/lesson-04/data/chr22_fragment.fa \
    --known-sites ~/ngs-analysis/lesson-04/data/dbsnp.vcf \
    -O results/post_recal_data.table
```

**参数解释：**
- `--known-sites`：已知变异位点，用于校正模型训练
- `--bqsr-recal-file`：校正表文件
- `-R`：参考基因组文件

**预期输出：**
校正表文件包含质量分数校正信息，用于后续变异检测。

**检查点：** 确认生成了 `sample.recal.bam` 文件和相应的索引文件。

---

### 步骤3：变异检测

#### 3.1 使用HaplotypeCaller检测变异

**操作说明：**
HaplotypeCaller是GATK的核心变异检测工具，基于单倍型重组装算法，能够同时检测SNP和InDel。

**执行命令：**
```bash
# 变异检测 - 生成GVCF文件
gatk HaplotypeCaller \
    -I results/sample.recal.bam \
    -R ~/ngs-analysis/lesson-04/data/chr22_fragment.fa \
    -O results/sample.g.vcf \
    -ERC GVCF \
    --native-pair-hmm-threads 4

# 基因分型 - 从GVCF生成最终VCF
gatk GenotypeGVCFs \
    -R ~/ngs-analysis/lesson-04/data/chr22_fragment.fa \
    -V results/sample.g.vcf \
    -O results/sample.vcf
```

**参数解释：**
- `-ERC GVCF`：生成GVCF格式输出，包含所有位点信息
- `--native-pair-hmm-threads`：使用多线程加速
- `-V`：输入的GVCF文件

**预期输出：**
```bash
# 查看VCF文件基本信息
bcftools stats results/sample.vcf | grep "number of records"
```

**结果验证：**
```bash
# 检查VCF文件格式
bcftools view -h results/sample.vcf | tail -5
bcftools view results/sample.vcf | head -10
```

#### 3.2 变异统计分析

**执行命令：**
```bash
# 生成详细统计报告
bcftools stats results/sample.vcf > results/variant_stats.txt

# 查看主要统计信息
grep "^SN" results/variant_stats.txt

# 统计变异类型
bcftools view -v snps results/sample.vcf | bcftools stats | grep "number of records"
bcftools view -v indels results/sample.vcf | bcftools stats | grep "number of records"
```

**检查点：** 对于5Mb基因组区域,期望检测到约500-1000个SNP和50-100个InDel变异。

**注意：** 变异数量取决于:
- 参考基因组区域大小 (本实验使用5Mb片段)
- 样本类型 (模拟数据的变异密度可能与真实数据不同)
- 测序深度 (~30X)

---

### 步骤4：变异质量过滤

#### 4.1 硬过滤方法

**操作说明：**
对于单样本分析，通常使用硬过滤方法，根据预设的质量标准过滤低质量变异。

**执行命令：**
```bash
# SNP硬过滤
gatk VariantFiltration \
    -R ~/ngs-analysis/lesson-04/data/chr22_fragment.fa \
    -V results/sample.vcf \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "SNP_filter" \
    -O results/sample.filtered.vcf

# 提取通过过滤的变异
bcftools view -f PASS results/sample.filtered.vcf > results/sample.pass.vcf

# 统计过滤结果
echo "原始变异数量:"
bcftools view results/sample.vcf | grep -v "^#" | wc -l

echo "通过过滤的变异数量:"
bcftools view results/sample.pass.vcf | grep -v "^#" | wc -l
```

**参数解释：**
- `QD < 2.0`：质量密度过低
- `FS > 60.0`：Fisher链偏倚过高
- `MQ < 40.0`：比对质量过低
- `MQRankSum < -12.5`：比对质量偏倚
- `ReadPosRankSum < -8.0`：reads位置偏倚

**预期输出示例：**
```
原始变异数量: 800-1200
通过过滤的变异数量: 600-900
过滤通过率: ~75-85%
```

**说明：**
- 本实验使用5Mb参考区域,变异数量远少于全基因组
- 过滤通过率应在75-85%之间
- 如果通过率过低(<60%),可能需要检查数据质量

#### 4.2 质量控制评估

**执行命令：**
```bash
# 计算Ti/Tv比值
python3 scripts/calculate_titv.py results/sample.pass.vcf

# 计算dbSNP重叠率
bcftools isec -p results/dbsnp_overlap results/sample.pass.vcf ~/ngs-analysis/lesson-04/data/dbsnp.vcf
echo "dbSNP重叠率:"
wc -l results/dbsnp_overlap/0002.vcf results/sample.pass.vcf
```

**检查点：**
- Ti/Tv比值: 应在1.8-2.2之间 (全基因组约2.0-2.1)
- dbSNP重叠率: 模拟数据约30-50% (真实数据应>95%)

**注意：** 由于使用模拟的dbSNP数据,重叠率会较低。真实项目中使用官方dbSNP数据库,重叠率应>95%。

---

### 步骤5：变异注释

#### 5.1 使用VEP进行注释

**操作说明：**
VEP (Variant Effect Predictor) 是Ensembl开发的变异注释工具，能够预测变异的功能影响。

**执行命令：**
```bash
# 使用VEP注释变异
vep --input_file results/sample.pass.vcf \
    --output_file results/sample.annotated.txt \
    --format vcf \
    --species homo_sapiens \
    --assembly GRCh38 \
    --cache \
    --offline \
    --everything \
    --force_overwrite \
    --stats_file results/vep_summary.html

# 转换为VCF格式
vep --input_file results/sample.pass.vcf \
    --output_file results/sample.annotated.vcf \
    --format vcf \
    --vcf \
    --species homo_sapiens \
    --assembly GRCh38 \
    --cache \
    --offline \
    --sift b \
    --polyphen b \
    --force_overwrite
```

**参数解释：**
- `--cache`：使用本地缓存数据库
- `--offline`：离线模式
- `--everything`：包含所有可用注释
- `--sift b`：包含SIFT功能预测
- `--polyphen b`：包含PolyPhen功能预测

**预期输出：**
注释文件包含基因名称、转录本信息、变异后果、功能预测等信息。

#### 5.2 注释结果分析

**执行命令：**
```bash
# 统计变异后果类型
grep -v "^#" results/sample.annotated.txt | cut -f7 | sort | uniq -c | sort -nr

# 提取有害变异
python3 scripts/filter_damaging_variants.py results/sample.annotated.txt > results/damaging_variants.txt

# 查看有害变异数量
wc -l results/damaging_variants.txt
```

**检查点：** 确认注释成功完成，并识别出潜在有害变异。

---

### 步骤6：结果可视化和解读

#### 6.1 生成质量控制图表

**执行命令：**
```bash
# 生成质量分布图
python3 scripts/plot_quality_distribution.py results/sample.pass.vcf plots/

# 生成Ti/Tv比值图
python3 scripts/plot_titv_ratio.py results/sample.pass.vcf plots/

# 生成变异密度图
python3 scripts/plot_variant_density.py results/sample.pass.vcf plots/
```

#### 6.2 IGV可视化验证

**操作说明：**
使用IGV (Integrative Genomics Viewer) 可视化验证重要变异位点。

**执行步骤：**
1. 启动IGV软件
2. 加载参考基因组 (hg38)
3. 加载BAM文件：`results/sample.recal.bam`
4. 加载VCF文件：`results/sample.pass.vcf`
5. 导航到感兴趣的变异位点进行验证

**验证要点：**
- 变异位点的reads支持情况
- 比对质量和碱基质量
- 是否存在比对偏倚
- 周围序列的复杂性

---

## 预期结果

### 主要输出文件
1. **变异检测结果**：`results/sample.pass.vcf`
   - 内容：高质量的SNP和InDel变异
   - 用途：后续分析和临床解读

2. **注释结果**：`results/sample.annotated.vcf`
   - 内容：包含功能注释的变异信息
   - 用途：变异功能影响评估

3. **质量控制报告**：`results/variant_stats.txt`
   - 内容：变异检测质量统计
   - 用途：评估数据质量

4. **可视化图表**：`plots/` 目录
   - 内容：质量分布、Ti/Tv比值等图表
   - 用途：质量控制和结果展示

### 关键结果指标
**本实验数据 (5Mb区域)：**
- **变异数量**：SNP ~500-1000个，InDel ~50-100个
- **Ti/Tv比值**：应该在 1.8-2.2 之间
- **dbSNP重叠率**：模拟数据 ~30-50% (真实数据应>95%)
- **过滤通过率**：应该 >75%
- **平均覆盖深度**：~30X

**全基因组规模 (参考)：**
- **变异数量**：SNP ~4-5百万个，InDel ~40-60万个
- **Ti/Tv比值**：~2.0-2.1
- **dbSNP重叠率**：>95%
- **过滤通过率**：>90%

### 成功标准
- [ ] 数据准备脚本成功运行,生成所有必需文件
- [ ] 所有GATK命令执行无错误
- [ ] 检测到预期数量的变异 (500-1000个SNP)
- [ ] Ti/Tv比值在正常范围内 (1.8-2.2)
- [ ] 变异过滤通过率达标 (>75%)
- [ ] (可选) 注释成功完成
- [ ] (可选) 可视化图表正常生成

## 故障排除

### 常见问题1：内存不足错误
**症状：** GATK运行时出现 "OutOfMemoryError"
**原因：** Java堆内存设置过小
**解决方案：**
```bash
# 增加Java内存设置
export JAVA_OPTS="-Xmx8g"
# 或者在GATK命令中添加
gatk --java-options "-Xmx8g" HaplotypeCaller ...
```

### 常见问题2：参考基因组索引缺失
**症状：** "Reference index not found" 错误
**原因：** 参考基因组缺少必要的索引文件
**解决方案：**
```bash
# 创建参考基因组索引
cd ~/ngs-analysis/lesson-04/data
samtools faidx chr22_fragment.fa
gatk CreateSequenceDictionary -R chr22_fragment.fa
```

### 常见问题3：VEP缓存数据缺失
**症状：** VEP提示缓存数据不存在
**原因：** 本地VEP缓存未安装
**解决方案：**
```bash
# 下载VEP缓存数据
vep_install -a cf -s homo_sapiens -y GRCh38 -c ~/.vep
```

### 常见问题4：变异数量异常
**症状：** 检测到的变异数量过多或过少
**原因：** 数据质量问题或参数设置不当
**解决方案：**
```bash
# 1. 检查BAM文件质量
samtools flagstat ~/ngs-analysis/lesson-04/data/sample.bam

# 2. 检查覆盖深度
samtools depth ~/ngs-analysis/lesson-04/data/sample.bam | \
    awk '{sum+=$3; count++} END {print "平均覆盖度:", sum/count}'

# 3. 验证参考基因组
head ~/ngs-analysis/lesson-04/data/chr22_fragment.fa

# 4. 检查VCF文件内容
bcftools view -H results/sample.vcf | head
```

**预期范围：**
- 5Mb区域: 500-1000个SNP
- 如果<100个: 可能覆盖度太低或数据有问题
- 如果>2000个: 可能样本质量差或有污染

### 获取帮助
如果遇到其他问题：
1. 检查错误日志：`cat logs/error.log`
2. 查看GATK帮助：`gatk HaplotypeCaller --help`
3. 联系助教或老师：wangys@hunau.edu.cn

## 扩展练习

### 练习1：多样本联合分析
**目标：** 学习多样本GVCF联合分析流程
**任务：** 
1. 准备多个样本的GVCF文件
2. 使用CombineGVCFs合并GVCF
3. 进行联合基因分型
4. 比较单样本和多样本结果差异
**提示：** 使用GenomicsDBImport可以提高大规模数据处理效率

### 练习2：VQSR软过滤
**目标：** 掌握VQSR机器学习过滤方法
**任务：**
1. 准备训练数据集（HapMap, 1000G等）
2. 使用VariantRecalibrator训练模型
3. 应用ApplyVQSR进行软过滤
4. 比较硬过滤和软过滤结果
**提示：** VQSR适用于大样本量数据，小样本建议使用硬过滤

### 练习3：结构变异检测
**目标：** 了解结构变异检测方法
**任务：**
1. 使用Manta检测结构变异
2. 使用DELLY进行SV分析
3. 比较不同工具的检测结果
4. 可视化验证大片段变异
**提示：** 结构变异检测需要考虑insert size和异常比对模式

### 思考问题
1. 为什么Ti/Tv比值可以作为变异质量的指标？
2. BQSR校正的原理是什么？为什么能提高变异检测准确性？
3. HaplotypeCaller相比传统方法有什么优势？
4. 如何区分真实的罕见变异和技术假象？
5. 在临床应用中，如何平衡敏感性和特异性？

## 参考资料

### 相关文献
1. McKenna A, et al. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res. 2010.
2. Van der Auwera GA, et al. From FastQ data to high confidence variant calls: the Genome Analysis Toolkit best practices pipeline. Curr Protoc Bioinformatics. 2013.
3. Poplin R, et al. A universal SNP and small-indel variant caller using deep neural networks. Nat Biotechnol. 2018.

### 在线资源
- GATK官方文档：https://gatk.broadinstitute.org/
- VEP用户指南：https://www.ensembl.org/info/docs/tools/vep/
- IGV用户手册：https://software.broadinstitute.org/software/igv/

### 软件文档
- GATK最佳实践：https://gatk.broadinstitute.org/hc/en-us/sections/360007226651
- bcftools手册：http://samtools.github.io/bcftools/
- VEP教程：https://www.ensembl.org/info/docs/tools/vep/script/vep_tutorial.html

## 附录

### 附录A：完整脚本文件
参见：`scripts/` 目录中的相关脚本文件
- `calculate_titv.py`：计算Ti/Tv比值
- `filter_damaging_variants.py`：筛选有害变异
- `plot_quality_distribution.py`：绘制质量分布图
- `plot_titv_ratio.py`：绘制Ti/Tv比值图
- `plot_variant_density.py`：绘制变异密度图

### 附录B：GATK参数配置
```bash
# 推荐的GATK内存设置
export GATK_JAVA_OPTIONS="-Xmx8g -XX:+UseParallelGC"

# HaplotypeCaller关键参数
--min-base-quality-score 20
--standard-min-confidence-threshold-for-calling 30
--native-pair-hmm-threads 4
```

### 附录C：VCF格式详解
VCF文件包含以下关键字段：
- CHROM：染色体
- POS：位置
- REF：参考等位基因
- ALT：变异等位基因
- QUAL：质量分数
- FILTER：过滤状态
- INFO：变异信息
- FORMAT：样本格式
- SAMPLE：样本基因型

---

**实验完成时间：** 预计 2 小时  
**难度等级：** 中级  
**最后更新：** 2025年1月