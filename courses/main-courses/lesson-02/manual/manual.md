# 第2次课实践操作手册：测序数据质量控制与预处理

## 课程信息
- **课程名称**：高通量测序数据分析
- **主讲教师**：王运生
- **联系邮箱**：wangys@hunau.edu.cn
- **办公室**：16教420室
- **上课地点**：105机房
- **课程时间**：2学时实践操作

## 实验目标

### 主要目标
- 掌握FastQC进行测序数据质量评估的方法
- 学会使用MultiQC整合多样本质量控制报告
- 熟练使用Trimmomatic进行数据清洗和预处理
- 能够解读质量控制报告并制定相应的处理策略

### 预期成果
- 生成完整的质量控制报告（FastQC + MultiQC）
- 完成测序数据的清洗和预处理
- 掌握质量控制前后效果对比的方法
- 建立标准化的质量控制工作流程

## 环境要求

### 软件环境
| 软件名称 | 版本要求 | 安装方式 | 说明 |
|---------|---------|---------|------|
| FastQC | v0.11.9+ | `sudo apt-get install fastqc` | 质量控制分析工具 |
| MultiQC | v1.12+ | `pip install multiqc` | 多样本报告整合 |
| Trimmomatic | v0.39+ | `conda install -c bioconda trimmomatic` | 序列修剪工具 |
| Python | 3.7+ | 系统自带 | 脚本运行环境 |
| Java | 8+ | `sudo apt-get install openjdk-8-jdk` | Trimmomatic依赖 |

### 硬件要求
- **内存**：至少 4 GB RAM
- **存储空间**：至少 10 GB 可用空间
- **CPU**：多核处理器（推荐4核以上）
- **网络**：稳定的互联网连接（用于下载数据）

### 数据准备
| 数据文件 | 大小 | 下载链接/位置 | 说明 |
|---------|------|-------------|------|
| sample_R1.fastq.gz | ~500MB | SRA:SRR1234567 | 双端测序数据第一端 |
| sample_R2.fastq.gz | ~500MB | SRA:SRR1234567 | 双端测序数据第二端 |
| TruSeq3-PE.fa | ~2KB | Trimmomatic安装目录 | Illumina接头序列文件 |

## 操作步骤

### 步骤1：环境设置和准备工作

#### 1.1 创建工作目录
```bash
# 创建本次实验的工作目录
mkdir -p ~/ngs-analysis/lesson-02
cd ~/ngs-analysis/lesson-02

# 创建子目录结构
mkdir -p {data,scripts,results,logs}
mkdir -p results/{fastqc_raw,fastqc_clean,multiqc,trimmed}
```

#### 1.2 检查软件环境
```bash
# 检查FastQC安装
fastqc --version
# 预期输出：FastQC v0.11.9

# 检查MultiQC安装
multiqc --version
# 预期输出：multiqc, version 1.12

# 检查Trimmomatic安装
trimmomatic -version
# 预期输出：0.39

# 检查Java环境
java -version
# 预期输出：openjdk version "1.8.0_xxx"
```

**预期输出：**
```
FastQC v0.11.9
multiqc, version 1.12
0.39
openjdk version "1.8.0_312"
```

#### 1.3 下载和准备数据
```bash
# 运行数据下载脚本
bash scripts/download_data.sh

# 验证数据完整性
ls -lh data/
md5sum data/*.fastq.gz > data/checksums.md5
```

**检查点：** 确认所有数据文件已正确下载并位于 `data/` 目录中，文件大小合理。

---

### 步骤2：原始数据质量评估

#### 2.1 使用FastQC分析原始数据

**操作说明：**
FastQC是最常用的测序数据质量控制工具，能够快速生成包含多个质量指标的HTML报告。我们首先对原始数据进行质量评估，了解数据的基本特征和潜在问题。

**执行命令：**
```bash
# 对原始数据运行FastQC分析
fastqc -t 4 -o results/fastqc_raw/ data/*.fastq.gz

# 检查生成的报告文件
ls -la results/fastqc_raw/
```bash

**参数解释：**
- `-t 4`：使用4个线程并行处理，提高分析速度
- `-o results/fastqc_raw/`：指定输出目录
- `data/*.fastq.gz`：输入的FASTQ文件（支持压缩格式）

**预期输出：**
```
results/fastqc_raw/
├── sample_R1_fastqc.html    # HTML格式报告
├── sample_R1_fastqc.zip     # 压缩的详细数据
├── sample_R2_fastqc.html
└── sample_R2_fastqc.zip
```

**结果验证：**
```bash
# 检查报告文件是否生成成功
find results/fastqc_raw/ -name "*.html" | wc -l
# 应该输出：2（对应R1和R2两个文件）
```

#### 2.2 查看FastQC报告

**操作说明：**
打开生成的HTML报告，重点关注以下几个关键模块的结果。

**查看报告：**
```bash
# 在浏览器中打开报告（如果有图形界面）
firefox results/fastqc_raw/sample_R1_fastqc.html &

# 或者使用文本方式查看关键统计信息
unzip -p results/fastqc_raw/sample_R1_fastqc.zip sample_R1_fastqc/fastqc_data.txt | head -20
```

**关键指标解读：**
- **Basic Statistics**：总序列数、序列长度、GC含量
- **Per base sequence quality**：每个位置的质量分布
- **Per sequence quality scores**：每条序列的质量分布
- **Sequence Duplication Levels**：序列重复水平
- **Overrepresented sequences**：过度表达的序列

**检查点：** 识别数据中存在的质量问题，如质量下降、接头污染、重复序列等。

---

### 步骤3：MultiQC整合分析

#### 3.1 生成MultiQC报告

**操作说明：**
MultiQC能够整合多个FastQC报告，生成一个综合的质量控制报告，便于比较不同样本的质量特征和识别批次效应。

**执行命令：**
```bash
# 运行MultiQC整合分析
multiqc -o results/multiqc/ results/fastqc_raw/

# 检查生成的报告
ls -la results/multiqc/
```bash

**参数解释：**
- `-o results/multiqc/`：指定输出目录
- `results/fastqc_raw/`：输入目录，包含FastQC结果

**预期输出：**
```
results/multiqc/
├── multiqc_report.html      # 主要的HTML报告
├── multiqc_data/           # 原始数据目录
│   ├── multiqc_fastqc.txt  # FastQC数据汇总
│   └── multiqc_general_stats.txt  # 总体统计
└── multiqc_plots/          # 图表文件
```

#### 3.2 解读MultiQC报告

**查看报告：**
```bash
# 打开MultiQC报告
firefox results/multiqc/multiqc_report.html &

# 查看汇总统计
cat results/multiqc/multiqc_data/multiqc_general_stats.txt
```

**重点关注：**
- **General Statistics**：样本概览和关键指标
- **FastQC: Sequence Quality Histograms**：质量分布对比
- **FastQC: Per Sequence GC Content**：GC含量分布
- **FastQC: Sequence Duplication Levels**：重复序列水平

**检查点：** 确认所有样本的质量特征，识别需要特殊处理的样本。

---

### 步骤4：数据清洗和预处理

#### 4.1 准备Trimmomatic接头文件

**操作说明：**
Trimmomatic需要接头序列文件来识别和去除接头污染。我们需要找到正确的接头文件或创建自定义的接头文件。

**查找接头文件：**
```bash
# 查找Trimmomatic安装目录中的接头文件
find /usr/share/trimmomatic/ -name "*.fa" 2>/dev/null
# 或者在conda环境中查找
find $CONDA_PREFIX/share/trimmomatic/ -name "*.fa" 2>/dev/null

# 如果找不到，创建常用的Illumina接头文件
cat > data/TruSeq3-PE.fa << 'EOF'
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
```

#### 4.2 运行Trimmomatic进行数据清洗

**操作说明：**
使用Trimmomatic对双端测序数据进行清洗，包括接头去除、低质量碱基修剪和长度过滤。

**执行命令：**
```bash
# 运行Trimmomatic进行数据清洗
trimmomatic PE -threads 4 -phred33 \
    data/sample_R1.fastq.gz data/sample_R2.fastq.gz \
    results/trimmed/sample_R1_paired.fastq.gz results/trimmed/sample_R1_unpaired.fastq.gz \
    results/trimmed/sample_R2_paired.fastq.gz results/trimmed/sample_R2_unpaired.fastq.gz \
    ILLUMINACLIP:data/TruSeq3-PE.fa:2:30:10 \
    LEADING:3 TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:36 \
    2> logs/trimmomatic.log
```

**参数详细解释：**
- `PE`：双端测序模式
- `-threads 4`：使用4个线程
- `-phred33`：质量编码格式
- `ILLUMINACLIP:data/TruSeq3-PE.fa:2:30:10`：
  - 接头文件路径
  - 种子错配数：2
  - 回文阈值：30
  - 简单clip阈值：10
- `LEADING:3`：去除5'端质量低于3的碱基
- `TRAILING:3`：去除3'端质量低于3的碱基
- `SLIDINGWINDOW:4:15`：滑动窗口大小4，平均质量阈值15
- `MINLEN:36`：保留长度至少36bp的序列

**预期输出：**
```
results/trimmed/
├── sample_R1_paired.fastq.gz    # 清洗后的配对R1序列
├── sample_R1_unpaired.fastq.gz  # 清洗后的非配对R1序列
├── sample_R2_paired.fastq.gz    # 清洗后的配对R2序列
└── sample_R2_unpaired.fastq.gz  # 清洗后的非配对R2序列
```

**结果验证：**
```bash
# 检查清洗统计信息
tail -10 logs/trimmomatic.log

# 比较清洗前后的文件大小
ls -lh data/*.fastq.gz results/trimmed/*.fastq.gz
```

#### 4.3 清洗后质量评估

**操作说明：**
对清洗后的数据重新运行FastQC分析，评估清洗效果。

**执行命令：**
```bash
# 对清洗后的配对数据运行FastQC
fastqc -t 4 -o results/fastqc_clean/ results/trimmed/*_paired.fastq.gz

# 生成清洗后的MultiQC报告
multiqc -o results/multiqc_clean/ results/fastqc_clean/
```

**检查点：** 确认清洗后数据质量有明显改善，特别是质量分数分布和接头污染情况。

---

### 步骤5：质量控制效果对比分析

#### 5.1 统计清洗前后的数据变化

**操作说明：**
使用脚本统计清洗前后的序列数量、质量分布等关键指标的变化。

**执行命令：**
```bash
# 运行质量对比脚本
python scripts/compare_quality.py \
    --raw-dir results/fastqc_raw/ \
    --clean-dir results/fastqc_clean/ \
    --output results/quality_comparison.txt

# 查看对比结果
cat results/quality_comparison.txt
```

**预期输出示例：**
```
Quality Control Comparison Report
================================
Sample: sample_R1
Raw sequences: 1,000,000
Clean sequences: 950,000 (95.0% retained)
Raw Q30 percentage: 75.2%
Clean Q30 percentage: 92.8%
Raw average length: 150 bp
Clean average length: 142 bp

Sample: sample_R2
Raw sequences: 1,000,000
Clean sequences: 950,000 (95.0% retained)
Raw Q30 percentage: 72.8%
Clean Q30 percentage: 91.5%
Raw average length: 150 bp
Clean average length: 140 bp
```

#### 5.2 生成综合质量报告

**执行命令：**
```bash
# 创建包含清洗前后对比的综合报告
multiqc -o results/final_report/ \
    results/fastqc_raw/ results/fastqc_clean/ \
    --filename "comprehensive_qc_report"

# 查看最终报告
firefox results/final_report/comprehensive_qc_report.html &
```

**检查点：** 综合报告应显示清洗前后的质量改善情况，确认处理效果符合预期。

---

### 步骤6：批量处理脚本编写

#### 6.1 创建自动化处理脚本

**操作说明：**
编写一个完整的批量处理脚本，能够自动完成从原始数据到清洗完成的整个质量控制流程。

**创建脚本：**
```bash
# 创建批量处理脚本
cat > scripts/batch_qc.sh << 'EOF'
#!/bin/bash

# 批量质量控制脚本
# 作者：王运生
# 用法：bash batch_qc.sh <input_dir> <output_dir>

set -e  # 遇到错误立即退出

# 参数检查
if [ $# -ne 2 ]; then
    echo "用法: $0 <input_dir> <output_dir>"
    echo "示例: $0 raw_data qc_results"
    exit 1
fi

INPUT_DIR=$1
OUTPUT_DIR=$2
THREADS=4

echo "开始批量质量控制处理..."
echo "输入目录: $INPUT_DIR"
echo "输出目录: $OUTPUT_DIR"

# 创建输出目录结构
mkdir -p ${OUTPUT_DIR}/{fastqc_raw,fastqc_clean,multiqc,trimmed,logs}

# 步骤1: 原始数据FastQC分析
echo "步骤1: 运行FastQC分析原始数据..."
fastqc -t $THREADS -o ${OUTPUT_DIR}/fastqc_raw/ ${INPUT_DIR}/*.fastq.gz

# 步骤2: 生成原始数据MultiQC报告
echo "步骤2: 生成原始数据MultiQC报告..."
multiqc -o ${OUTPUT_DIR}/multiqc/ ${OUTPUT_DIR}/fastqc_raw/ \
    --filename "raw_data_report"

# 步骤3: 数据清洗
echo "步骤3: 执行数据清洗..."
for r1_file in ${INPUT_DIR}/*_R1.fastq.gz; do
    if [ -f "$r1_file" ]; then
        base=$(basename ${r1_file} _R1.fastq.gz)
        r2_file=${INPUT_DIR}/${base}_R2.fastq.gz
        
        if [ -f "$r2_file" ]; then
            echo "处理样本: $base"
            trimmomatic PE -threads $THREADS -phred33 \
                $r1_file $r2_file \
                ${OUTPUT_DIR}/trimmed/${base}_R1_paired.fastq.gz \
                ${OUTPUT_DIR}/trimmed/${base}_R1_unpaired.fastq.gz \
                ${OUTPUT_DIR}/trimmed/${base}_R2_paired.fastq.gz \
                ${OUTPUT_DIR}/trimmed/${base}_R2_unpaired.fastq.gz \
                ILLUMINACLIP:data/TruSeq3-PE.fa:2:30:10 \
                LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
                2>> ${OUTPUT_DIR}/logs/trimmomatic_${base}.log
        else
            echo "警告: 找不到配对文件 $r2_file"
        fi
    fi
done

# 步骤4: 清洗后数据FastQC分析
echo "步骤4: 分析清洗后数据..."
fastqc -t $THREADS -o ${OUTPUT_DIR}/fastqc_clean/ \
    ${OUTPUT_DIR}/trimmed/*_paired.fastq.gz

# 步骤5: 生成最终综合报告
echo "步骤5: 生成综合质量报告..."
multiqc -o ${OUTPUT_DIR}/multiqc/ \
    ${OUTPUT_DIR}/fastqc_raw/ ${OUTPUT_DIR}/fastqc_clean/ \
    --filename "comprehensive_qc_report"

echo "批量质量控制处理完成！"
echo "查看报告: ${OUTPUT_DIR}/multiqc/comprehensive_qc_report.html"
EOF

# 给脚本添加执行权限
chmod +x scripts/batch_qc.sh
```

#### 6.2 测试批量处理脚本

**执行命令：**
```bash
# 测试批量处理脚本
bash scripts/batch_qc.sh data results_batch

# 检查处理结果
ls -la results_batch/
```bash

**检查点：** 脚本应该能够成功处理所有数据文件，生成完整的质量控制报告。

---

## 预期结果

### 主要输出文件
1. **FastQC报告**：`results/fastqc_raw/` 和 `results/fastqc_clean/`
   - 内容：每个FASTQ文件的详细质量分析报告
   - 用途：识别数据质量问题和评估清洗效果

2. **MultiQC报告**：`results/multiqc/comprehensive_qc_report.html`
   - 内容：整合的多样本质量控制报告
   - 用途：批量样本质量比较和趋势分析

3. **清洗后数据**：`results/trimmed/*_paired.fastq.gz`
   - 内容：经过质量控制和清洗的高质量测序数据
   - 用途：后续生物信息学分析的输入数据

4. **处理日志**：`logs/trimmomatic.log`
   - 内容：数据清洗过程的详细日志
   - 用途：问题诊断和处理效果评估

### 关键结果指标
- **序列保留率**：应该在 85-95% 之间
- **Q30比例提升**：清洗后应提升至 90% 以上
- **平均质量分数**：应有明显提升
- **接头污染**：应完全去除或显著降低

### 成功标准
- [ ] 所有FastQC和MultiQC报告成功生成
- [ ] 数据清洗完成，保留率在合理范围内
- [ ] 清洗后质量指标有明显改善
- [ ] 批量处理脚本能够正常运行
- [ ] 生成的数据适合后续分析使用

## 故障排除

### 常见问题1：FastQC运行失败
**症状：** 
```
Exception in thread "main" java.lang.OutOfMemoryError: Java heap space
```
**原因：** Java内存不足
**解决方案：**
```bash
# 增加Java内存限制
export _JAVA_OPTIONS="-Xmx4g"
fastqc -t 4 -o results/fastqc_raw/ data/*.fastq.gz
```

### 常见问题2：Trimmomatic找不到接头文件
**症状：** 
```
Exception: Could not find adapter file: TruSeq3-PE.fa
```
**原因：** 接头文件路径不正确
**解决方案：**
```bash
# 使用绝对路径或确保文件存在
ls -la data/TruSeq3-PE.fa
# 或者使用系统安装的接头文件
find /usr -name "TruSeq3-PE.fa" 2>/dev/null
```

### 常见问题3：MultiQC报告为空
**症状：** MultiQC报告中没有数据
**原因：** FastQC结果目录路径错误或文件损坏
**解决方案：**
```bash
# 检查FastQC结果文件
ls -la results/fastqc_raw/*.zip
# 重新运行FastQC
rm -rf results/fastqc_raw/*
fastqc -t 4 -o results/fastqc_raw/ data/*.fastq.gz
```

### 常见问题4：清洗后序列数量过少
**症状：** 清洗后保留的序列数量不足50%
**原因：** 清洗参数过于严格或数据质量极差
**解决方案：**
```bash
# 放宽清洗参数
trimmomatic PE -threads 4 -phred33 \
    data/sample_R1.fastq.gz data/sample_R2.fastq.gz \
    results/trimmed/sample_R1_paired.fastq.gz results/trimmed/sample_R1_unpaired.fastq.gz \
    results/trimmed/sample_R2_paired.fastq.gz results/trimmed/sample_R2_unpaired.fastq.gz \
    ILLUMINACLIP:data/TruSeq3-PE.fa:2:30:10 \
    LEADING:2 TRAILING:2 \
    SLIDINGWINDOW:4:10 \
    MINLEN:30
```bash

### 获取帮助
如果遇到其他问题：
1. 检查错误日志：`cat logs/trimmomatic.log`
2. 查看软件帮助：`fastqc --help`, `multiqc --help`
3. 联系助教或老师：wangys@hunau.edu.cn

## 扩展练习

### 练习1：参数优化实验
**目标：** 比较不同Trimmomatic参数对清洗效果的影响
**任务：** 
1. 使用三组不同的参数设置清洗同一份数据
2. 比较清洗后的质量指标和保留率
3. 分析参数选择对结果的影响

**提示：** 
- 严格参数：SLIDINGWINDOW:4:20, MINLEN:50
- 中等参数：SLIDINGWINDOW:4:15, MINLEN:36
- 宽松参数：SLIDINGWINDOW:4:10, MINLEN:25

### 练习2：批量样本质量分析
**目标：** 分析多个样本的质量特征，识别异常样本
**任务：** 
1. 下载5-10个不同的测序样本
2. 进行批量质量控制分析
3. 识别质量异常的样本并分析原因

**提示：** 可以从SRA数据库下载不同项目的数据进行比较

### 练习3：自定义质量控制流程
**目标：** 根据特定需求设计质量控制流程
**任务：** 
1. 针对特定类型的测序数据（如RNA-seq, ChIP-seq）
2. 设计相应的质量控制策略
3. 编写自动化处理脚本

**提示：** 不同类型的测序数据可能需要不同的质量标准和处理策略

### 思考问题
1. 为什么测序数据的质量会在读长末端下降？如何在实验设计阶段减少这种现象？
2. 在什么情况下应该选择更严格的质量过滤标准？什么情况下可以放宽标准？
3. 如何平衡数据质量和数据量？过度清洗可能带来什么问题？

## 参考资料

### 相关文献
1. Andrews, S. (2010). FastQC: A quality control tool for high throughput sequence data. Babraham Bioinformatics.
2. Ewels, P., et al. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047-3048.
3. Bolger, A.M., et al. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics, 30(15), 2114-2120.
4. Chen, S., et al. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17), i884-i890.

### 在线资源
- FastQC官方文档：https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- MultiQC官方文档：https://multiqc.info/
- Trimmomatic用户手册：http://www.usadellab.org/cms/?page=trimmomatic
- FASTQ格式说明：https://en.wikipedia.org/wiki/FASTQ_format

### 软件文档
- FastQC Help文档：https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/
- MultiQC模块文档：https://multiqc.info/docs/#multiqc-modules
- Trimmomatic参数说明：http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf

## 附录

### 附录A：完整脚本文件
参见：`scripts/` 目录中的相关脚本文件
- `setup.sh`：环境设置脚本
- `download_data.sh`：数据下载脚本
- `run_fastqc.sh`：FastQC批量运行脚本
- `run_multiqc.py`：MultiQC分析脚本
- `quality_trim.sh`：数据清洗脚本
- `compare_quality.py`：质量对比脚本

### 附录B：Trimmomatic参数详解
```bash
# 接头去除参数
ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>

# 质量修剪参数
LEADING:<quality>     # 去除5'端低质量碱基
TRAILING:<quality>    # 去除3'端低质量碱基
SLIDINGWINDOW:<windowSize>:<requiredQuality>  # 滑动窗口修剪
MAXINFO:<targetLength>:<strictness>  # 最大信息修剪

# 长度过滤参数
MINLEN:<length>       # 最小长度过滤
MAXLEN:<length>       # 最大长度过滤
```

### 附录C：质量分数对照表
| Phred Score | 错误概率 | 准确率 | 质量等级 |
|-------------|----------|--------|----------|
| Q10 | 1 in 10 | 90% | 低 |
| Q20 | 1 in 100 | 99% | 中等 |
| Q30 | 1 in 1,000 | 99.9% | 高 |
| Q40 | 1 in 10,000 | 99.99% | 极高 |

---

**实验完成时间：** 预计 2 小时  
**难度等级：** 中级  
**最后更新：** 2025年1月