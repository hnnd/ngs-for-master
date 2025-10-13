# 第3次课实践操作手册：高通量测序序列比对算法与工具

## 课程信息
- **课程名称**：高通量测序数据分析
- **主讲教师**：王运生
- **联系邮箱**：wangys@hunau.edu.cn
- **上课地点**：105机房

## 实验目标

通过本次实践操作，学生将能够：

1. **掌握比对工具使用**
   - 学会构建参考基因组索引
   - 熟练使用BWA-MEM进行序列比对
   - 掌握Bowtie2的参数设置和使用方法
   - 了解不同工具的性能特点

2. **理解比对结果处理**
   - 掌握SAM/BAM文件格式和内容
   - 学会使用samtools进行文件处理
   - 能够进行比对质量统计和评估
   - 掌握比对结果的可视化方法

3. **优化比对策略**
   - 了解不同参数对比对结果的影响
   - 学会根据数据特点选择合适的工具和参数
   - 掌握比对质量控制的方法和标准

## 环境要求

### 软件环境
- **BWA** (>= 0.7.17)
- **Bowtie2** (>= 2.4.0)
- **samtools** (>= 1.10)
- **Python 3.8+** (pandas, matplotlib, seaborn)
- **IGV** (Integrative Genomics Viewer)

### 硬件要求
- **内存**：至少8GB RAM
- **存储**：至少10GB可用空间
- **CPU**：多核处理器（推荐4核以上）

### 数据准备
本次实验将使用以下数据：
- **参考基因组**：人类基因组22号染色体片段 (5Mb)
  - 来源：Ensembl/UCSC公共数据库
  - 版本：GRCh38 (hg38)
- **测序数据**：Paired-end reads (2×100bp, 1M reads)
  - 选项1：NCBI SRA真实数据 (SRR622461)
  - 选项2：使用wgsim工具生成的模拟数据
  - 选项3：课程服务器提供的数据（如果可用）

> **📖 数据详细说明**：关于数据来源、格式要求和使用建议的完整信息，请参阅 [DATA_SOURCES.md](DATA_SOURCES.md)

## 操作步骤

### 步骤1：环境设置和数据准备

#### 1.1 创建工作目录
```bash
# 创建实验目录
mkdir -p ~/ngs_alignment_lab
cd ~/ngs_alignment_lab

# 创建子目录
mkdir -p data reference results logs scripts
```

#### 1.2 检查软件安装（重要！）
```bash
# 检查BWA版本
bwa 2>&1 | head -3

# 检查Bowtie2版本
bowtie2 --version | head -1

# 检查samtools版本
samtools --version | head -1

# 如果缺少工具，使用conda安装
# conda install -c bioconda bwa bowtie2 samtools
```

#### 1.3 快速数据准备（推荐）

**使用自动化脚本准备所有数据**
```bash
# 下载数据准备脚本（如果还没有）
# 脚本应该在课程材料的scripts目录中

# 给脚本添加执行权限
chmod +x scripts/prepare_data.sh

# 运行数据准备脚本（自动下载所有数据）
bash scripts/prepare_data.sh

# 该脚本将自动完成：
# 1. 检查必需的软件
# 2. 下载参考基因组（Ensembl或UCSC）
# 3. 下载/生成测序数据
# 4. 验证数据完整性
```

**脚本使用选项**
```bash
# 仅下载参考基因组
bash scripts/prepare_data.sh --reference

# 仅准备测序数据
bash scripts/prepare_data.sh --sequencing

# 仅验证数据
bash scripts/prepare_data.sh --validate

# 查看帮助
bash scripts/prepare_data.sh --help
```

#### 1.4 手动数据准备（可选）

如果自动脚本无法运行，或者您想了解详细步骤，可以手动准备数据：

**方案一：下载真实参考基因组**
```bash
# 从Ensembl下载人类22号染色体序列
echo "下载人类22号染色体参考序列..."
wget -O reference/chr22.fa.gz \
  "http://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz"

# 解压
gunzip reference/chr22.fa.gz

# 提取前5Mb作为实验数据（减少计算时间）
echo "提取前5Mb序列..."
samtools faidx reference/chr22.fa
samtools faidx reference/chr22.fa 22:1-5000000 > reference/chr22_fragment.fa

# 清理临时文件
rm reference/chr22.fa reference/chr22.fa.fai
```

**方案二：使用UCSC数据源（备选）**
```bash
# 从UCSC下载
wget -O reference/chr22.fa.gz \
  "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz"

gunzip reference/chr22.fa.gz
# 提取前5Mb
samtools faidx reference/chr22.fa
samtools faidx reference/chr22.fa chr22:1-5000000 > reference/chr22_fragment.fa
rm reference/chr22.fa reference/chr22.fa.fai
```

**测序数据准备**

由于公共数据库的原始测序数据通常非常大，本实验提供两种方式获取测序数据：

**选项A：下载示例数据集（推荐用于教学）**
```bash
# 从NCBI SRA下载一个小的测序数据集
# 使用SRA Toolkit下载数据（需要先安装sra-toolkit）
# 这里使用一个小规模的人类基因组测序数据

# 安装sra-toolkit（如果未安装）
# conda install -c bioconda sra-tools

# 下载一个小规模数据集（约100万条reads）
echo "下载测序数据..."
fastq-dump --split-files --gzip -X 1000000 SRR622461 -O data/

# 重命名文件
mv data/SRR622461_1.fastq.gz data/sample_R1.fastq.gz
mv data/SRR622461_2.fastq.gz data/sample_R2.fastq.gz

# 解压
gunzip data/sample_R1.fastq.gz
gunzip data/sample_R2.fastq.gz
```

**选项B：生成模拟测序数据（如果网络受限）**
```bash
# 如果无法下载真实数据，可以使用工具生成模拟数据
# 需要先下载参考基因组片段

# 使用wgsim生成模拟reads（需要先安装wgsim）
# conda install -c bioconda wgsim

echo "生成模拟测序数据..."
wgsim -N 1000000 -1 100 -2 100 -r 0.001 -R 0.15 -X 0.3 \
  reference/chr22_fragment.fa \
  data/sample_R1.fastq \
  data/sample_R2.fastq

echo "数据生成完成"
```

**选项C：使用课程提供的数据（如果可用）**
```bash
# 如果老师提供了本地数据服务器
# 请替换下面的URL为实际的课程数据地址
# wget -O data/sample_R1.fastq.gz "http://course-server/data/sample_R1.fastq.gz"
# wget -O data/sample_R2.fastq.gz "http://course-server/data/sample_R2.fastq.gz"
# gunzip data/*.fastq.gz
```

**验证数据完整性**
```bash
# 检查参考基因组
echo "检查参考基因组..."
head -1 reference/chr22_fragment.fa
wc -l reference/chr22_fragment.fa

# 检查测序数据
echo "检查测序数据..."
head -4 data/sample_R1.fastq
echo "R1 reads数量: $(wc -l < data/sample_R1.fastq | awk '{print $1/4}')"
echo "R2 reads数量: $(wc -l < data/sample_R2.fastq | awk '{print $1/4}')"

# 确认数据格式正确
if [ -f data/sample_R1.fastq ] && [ -f data/sample_R2.fastq ]; then
    echo "✓ 测序数据准备完成"
else
    echo "✗ 测序数据文件缺失，请检查"
fi
```

### 步骤2：参考基因组索引构建

#### 2.1 构建BWA索引
```bash
# 进入参考基因组目录
cd reference

# 构建BWA索引
echo "开始构建BWA索引..."
time bwa index chr22_fragment.fa

# 检查生成的索引文件
ls -lh chr22_fragment.fa.*
```

**预期输出**：
```
chr22_fragment.fa.amb
chr22_fragment.fa.ann  
chr22_fragment.fa.bwt
chr22_fragment.fa.pac
chr22_fragment.fa.sa
```

#### 2.2 构建Bowtie2索引
```bash
# 构建Bowtie2索引
echo "开始构建Bowtie2索引..."
time bowtie2-build chr22_fragment.fa chr22_bt2

# 检查生成的索引文件
ls -lh chr22_bt2.*
```

**预期输出**：
```
chr22_bt2.1.bt2
chr22_bt2.2.bt2
chr22_bt2.3.bt2
chr22_bt2.4.bt2
chr22_bt2.rev.1.bt2
chr22_bt2.rev.2.bt2
```

#### 2.3 比较索引文件大小
```bash
# 统计索引文件大小
echo "BWA索引文件大小："
du -sh chr22_fragment.fa.*

echo "Bowtie2索引文件大小："
du -sh chr22_bt2.*

# 返回工作目录
cd ..
```

### 步骤3：使用BWA进行序列比对

#### 3.1 基本BWA-MEM比对
```bash
# 使用默认参数进行比对
echo "开始BWA-MEM比对..."
time bwa mem \
  -t 4 \
  reference/chr22_fragment.fa \
  data/sample_R1.fastq \
  data/sample_R2.fastq \
  > results/bwa_default.sam

# 检查SAM文件
echo "BWA比对完成，检查结果："
wc -l results/bwa_default.sam
head -20 results/bwa_default.sam
```

#### 3.2 BWA参数优化比对
```bash
# 使用优化参数进行比对
echo "使用优化参数进行BWA比对..."
time bwa mem \
  -t 4 \
  -k 19 \
  -w 100 \
  -d 100 \
  -r 1.5 \
  -A 1 -B 4 -O 6 -E 1 \
  -M \
  reference/chr22_fragment.fa \
  data/sample_R1.fastq \
  data/sample_R2.fastq \
  > results/bwa_optimized.sam

echo "优化参数BWA比对完成"
```

#### 3.3 转换为BAM格式并排序
```bash
# 转换SAM为BAM并排序
echo "转换SAM为BAM格式..."
samtools view -bS results/bwa_default.sam | \
samtools sort -o results/bwa_default_sorted.bam

samtools view -bS results/bwa_optimized.sam | \
samtools sort -o results/bwa_optimized_sorted.bam

# 建立索引
samtools index results/bwa_default_sorted.bam
samtools index results/bwa_optimized_sorted.bam

echo "BAM文件转换和索引完成"
```

### 步骤4：使用Bowtie2进行序列比对

#### 4.1 使用不同敏感性模式
```bash
# 快速模式
echo "Bowtie2快速模式比对..."
time bowtie2 \
  --fast \
  -p 4 \
  -x reference/chr22_bt2 \
  -1 data/sample_R1.fastq \
  -2 data/sample_R2.fastq \
  -S results/bowtie2_fast.sam

# 敏感模式（默认）
echo "Bowtie2敏感模式比对..."
time bowtie2 \
  --sensitive \
  -p 4 \
  -x reference/chr22_bt2 \
  -1 data/sample_R1.fastq \
  -2 data/sample_R2.fastq \
  -S results/bowtie2_sensitive.sam

# 高敏感模式
echo "Bowtie2高敏感模式比对..."
time bowtie2 \
  --very-sensitive \
  -p 4 \
  -x reference/chr22_bt2 \
  -1 data/sample_R1.fastq \
  -2 data/sample_R2.fastq \
  -S results/bowtie2_very_sensitive.sam
```

#### 4.2 局部比对模式
```bash
# 局部比对模式
echo "Bowtie2局部比对模式..."
time bowtie2 \
  --local \
  --very-sensitive-local \
  -p 4 \
  -x reference/chr22_bt2 \
  -1 data/sample_R1.fastq \
  -2 data/sample_R2.fastq \
  -S results/bowtie2_local.sam
```

#### 4.3 自定义参数比对
```bash
# 自定义参数
echo "Bowtie2自定义参数比对..."
time bowtie2 \
  -p 4 \
  -L 20 \
  -i S,1,0.50 \
  --mp 6,2 \
  --np 1 \
  --rdg 5,3 \
  --rfg 5,3 \
  -I 50 \
  -X 800 \
  -x reference/chr22_bt2 \
  -1 data/sample_R1.fastq \
  -2 data/sample_R2.fastq \
  -S results/bowtie2_custom.sam
```

#### 4.4 转换Bowtie2结果为BAM
```bash
# 转换所有Bowtie2结果为BAM
echo "转换Bowtie2结果为BAM..."
for sam_file in results/bowtie2_*.sam; do
    base_name=$(basename "$sam_file" .sam)
    samtools view -bS "$sam_file" | \
    samtools sort -o "results/${base_name}_sorted.bam"
    samtools index "results/${base_name}_sorted.bam"
done

echo "Bowtie2 BAM转换完成"
```

### 步骤5：比对结果统计和分析

#### 5.1 基本比对统计
```bash
# 运行比对统计脚本
echo "生成比对统计报告..."
python scripts/alignment_stats.py

# 查看统计结果
cat results/alignment_statistics.txt
```

#### 5.2 详细质量分析
```bash
# 使用samtools统计
echo "详细比对质量分析..."

for bam_file in results/*_sorted.bam; do
    echo "分析文件: $bam_file"
    
    # 基本统计
    samtools flagstat "$bam_file" > "${bam_file%.bam}_flagstat.txt"
    
    # 插入片段统计
    samtools stats "$bam_file" > "${bam_file%.bam}_stats.txt"
    
    # MAPQ分布
    samtools view "$bam_file" | cut -f5 | sort -n | uniq -c > "${bam_file%.bam}_mapq.txt"
    
    echo "完成: $bam_file"
done
```

#### 5.3 覆盖度分析
```bash
# 计算覆盖度
echo "计算基因组覆盖度..."

for bam_file in results/*_sorted.bam; do
    base_name=$(basename "$bam_file" _sorted.bam)
    
    # 计算每个位点的覆盖度
    samtools depth "$bam_file" > "results/${base_name}_coverage.txt"
    
    # 统计覆盖度分布
    awk '{
        sum += $3; 
        count++; 
        if($3 > 0) covered++
    } END {
        print "平均覆盖度:", sum/count
        print "覆盖率:", covered/count*100"%"
    }' "results/${base_name}_coverage.txt" > "results/${base_name}_coverage_summary.txt"
done
```

### 步骤6：结果可视化

#### 6.1 生成比对统计图表
```bash
# 运行可视化脚本
echo "生成可视化图表..."
python scripts/visualize_results.py

# 查看生成的图片
ls -la results/*.png
```

#### 6.2 使用IGV查看比对结果
```bash
# 准备IGV查看
echo "准备IGV可视化文件..."

# 确保所有BAM文件都有索引
for bam_file in results/*_sorted.bam; do
    if [ ! -f "${bam_file}.bai" ]; then
        samtools index "$bam_file"
    fi
done

echo "IGV文件准备完成"
echo "请打开IGV，加载参考基因组和BAM文件进行可视化"
```

### 步骤7：性能比较和分析

#### 7.1 运行时间比较
```bash
# 提取运行时间信息
echo "比对工具性能比较："
echo "========================"

# 从日志文件中提取时间信息（如果有的话）
if [ -f logs/timing.log ]; then
    cat logs/timing.log
else
    echo "请查看之前运行时显示的时间信息"
fi
```

#### 7.2 结果质量比较
```bash
# 比较不同工具的比对质量
echo "比对质量比较："
echo "========================"

for flagstat_file in results/*_flagstat.txt; do
    echo "文件: $(basename $flagstat_file _flagstat.txt)"
    head -5 "$flagstat_file"
    echo "------------------------"
done
```

#### 7.3 生成综合报告
```bash
# 生成最终比较报告
python scripts/generate_report.py > results/final_report.txt

echo "综合分析报告已生成: results/final_report.txt"
cat results/final_report.txt
```

## 预期结果

### 比对统计预期值

#### BWA-MEM结果
- **总reads数**：2,000,000 (1M pairs)
- **比对率**：~98-99%
- **唯一比对率**：~95-97%
- **正确配对率**：~96-98%
- **平均MAPQ**：~35-40

#### Bowtie2结果
- **快速模式**：比对率 ~96-97%，速度最快
- **敏感模式**：比对率 ~97-98%，平衡性能
- **高敏感模式**：比对率 ~98-99%，速度较慢
- **局部模式**：软剪切比例增加

### 文件大小预期
```
原始FASTQ文件: ~400MB (压缩前)
SAM文件: ~600-800MB
BAM文件: ~150-200MB
索引文件: ~50-100MB
```

### 性能预期
```
BWA索引构建: 2-5分钟
Bowtie2索引构建: 3-8分钟
BWA-MEM比对: 5-10分钟
Bowtie2比对: 8-15分钟
```

## 故障排除

### 常见问题及解决方案

#### 1. 内存不足错误
**问题**：比对过程中出现内存不足
**解决方案**：
```bash
# 减少线程数
bwa mem -t 2 ...  # 从4减少到2

# 或者分批处理数据
split -l 1000000 sample_R1.fastq sample_R1_part_
```

#### 2. 索引文件损坏
**问题**：索引文件不完整或损坏
**解决方案**：
```bash
# 删除所有索引文件并重新构建
rm reference/chr22_fragment.fa.*
bwa index reference/chr22_fragment.fa
```

#### 3. SAM文件格式错误
**问题**：SAM文件无法转换为BAM
**解决方案**：
```bash
# 检查SAM文件头部
head -20 results/sample.sam

# 使用samtools检查
samtools view -H results/sample.sam
```

#### 4. 权限问题
**问题**：无法写入结果文件
**解决方案**：
```bash
# 检查目录权限
ls -la results/

# 修改权限
chmod 755 results/
```

### 性能优化建议

#### 1. 多线程设置
```bash
# 根据CPU核心数设置线程
nproc  # 查看可用核心数
bwa mem -t $(nproc) ...  # 使用所有核心
```

#### 2. 内存优化
```bash
# 监控内存使用
free -h
top -p $(pgrep bwa)
```

#### 3. 磁盘I/O优化
```bash
# 使用SSD存储临时文件
export TMPDIR=/path/to/ssd/tmp
```

## 扩展练习

### 练习1：参数敏感性分析
比较不同BWA参数对比对结果的影响：
```bash
# 测试不同种子长度
for k in 15 19 25; do
    bwa mem -k $k -t 4 reference/chr22_fragment.fa \
        data/sample_R1.fastq data/sample_R2.fastq \
        > results/bwa_k${k}.sam
done
```

### 练习2：错误容忍度测试
在测序数据中引入人工错误，测试比对工具的容错能力：
```bash
# 使用脚本引入错误
python scripts/introduce_errors.py \
    --input data/sample_R1.fastq \
    --output data/sample_R1_errors.fastq \
    --error_rate 0.05
```

### 练习3：不同数据类型比对
尝试比对不同类型的测序数据：
- 单端测序数据
- 不同读长的数据
- 不同插入片段大小的数据

### 练习4：比对结果过滤
学习如何过滤低质量的比对结果：
```bash
# 过滤MAPQ < 20的比对
samtools view -q 20 -b results/bwa_default_sorted.bam \
    > results/bwa_filtered.bam

# 过滤未正确配对的reads
samtools view -f 2 -b results/bwa_default_sorted.bam \
    > results/bwa_proper_pairs.bam
```

## 思考问题

1. **算法选择**：在什么情况下应该选择BWA而不是Bowtie2？
2. **参数优化**：如何根据数据特点调整比对参数？
3. **质量控制**：如何判断比对结果的质量是否可接受？
4. **性能平衡**：如何在比对精度和运行速度之间找到平衡？
5. **下游分析**：比对结果的质量如何影响后续的变异检测？

## 参考资料

### 官方文档
1. [BWA Manual](http://bio-bwa.sourceforge.net/bwa.shtml)
2. [Bowtie2 Manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
3. [SAMtools Documentation](http://www.htslib.org/doc/samtools.html)
4. [SAM Format Specification](https://samtools.github.io/hts-specs/SAMv1.pdf)

### 重要文献
1. Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler transform. *Bioinformatics*, 25, 1754-1760.
2. Langmead B. and Salzberg S. (2012) Fast gapped-read alignment with Bowtie 2. *Nature Methods*, 9, 357-359.
3. Li H. et al. (2009) The Sequence Alignment/Map format and SAMtools. *Bioinformatics*, 25, 2078-2079.

### 在线资源
1. [Galaxy Training Materials](https://training.galaxyproject.org/)
2. [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651)
3. [Bioinformatics Workbook](https://bioinformaticsworkbook.org/)

---

**实验完成后，请保存所有结果文件，并准备在下次课程中使用这些比对结果进行变异检测分析。**