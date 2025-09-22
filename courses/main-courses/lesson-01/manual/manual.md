# 第1次课实践操作手册：测序技术原理与平台比较

## 课程信息
- **课程名称**：高通量测序数据分析
- **主讲教师**：王运生
- **联系邮箱**：wangys@hunau.edu.cn
- **办公室**：16教420室
- **上课地点**：105机房

## 实验目标
- 熟悉不同测序平台产生的数据格式和特征
- 学会分析测序数据的基本质量指标
- 掌握测序成本计算和平台选择方法
- 理解不同测序技术的优缺点和适用场景

## 环境要求

### 软件环境
| 软件名称 | 版本要求 | 安装方式 |
|---------|---------|---------|
| Python | 3.7+ | conda install python |
| pandas | 1.3+ | pip install pandas |
| matplotlib | 3.5+ | pip install matplotlib |
| seaborn | 0.11+ | pip install seaborn |
| biopython | 1.79+ | pip install biopython |

### 硬件要求
- **内存**：至少 4 GB RAM
- **存储空间**：至少 2 GB 可用空间
- **CPU**：双核以上处理器

### 数据准备
| 数据文件 | 大小 | 说明 |
|---------|------|------|
| illumina_sample.fastq | 50 MB | Illumina测序数据示例 |
| pacbio_sample.fastq | 30 MB | PacBio测序数据示例 |
| nanopore_sample.fastq | 40 MB | Nanopore测序数据示例 |
#
# 操作步骤

### 步骤1：环境设置和准备工作

#### 1.1 创建工作目录
```bash
# 创建本次实验的工作目录
mkdir -p ~/ngs-analysis/lesson-01
cd ~/ngs-analysis/lesson-01

# 创建子目录结构
mkdir -p {data,scripts,results,logs}
```

#### 1.2 检查软件环境
```bash
# 检查Python环境
python --version
pip list | grep -E "(pandas|matplotlib|seaborn|biopython)"

# 检查必要的Python包
python -c "import pandas, matplotlib, seaborn, Bio; print('所有包已安装')"
```

**预期输出：**
```
Python 3.8.x
所有包已安装
```

#### 1.3 下载和准备数据
```bash
# 下载实验数据（模拟数据）
cd data
wget -O illumina_sample.fastq "https://example.com/illumina_sample.fastq"
wget -O pacbio_sample.fastq "https://example.com/pacbio_sample.fastq"
wget -O nanopore_sample.fastq "https://example.com/nanopore_sample.fastq"

# 验证数据完整性
ls -lh *.fastq
```

**检查点：** 确认所有数据文件已正确下载并位于 `data/` 目录中。

---

### 步骤2：FASTQ格式数据分析

#### 2.1 理解FASTQ格式

**操作说明：**
FASTQ格式是测序数据的标准格式，包含序列信息和质量信息。每个序列记录包含4行：
1. 序列标识符（以@开头）
2. DNA序列
3. 分隔符（以+开头）
4. 质量值（ASCII编码）

**执行命令：**
```bash
# 查看FASTQ文件的前几行
head -20 data/illumina_sample.fastq
head -20 data/pacbio_sample.fastq
head -20 data/nanopore_sample.fastq
```

**结果验证：**
```bash
# 统计每个文件的序列数量
echo "Illumina序列数量:"
grep -c "^@" data/illumina_sample.fastq
echo "PacBio序列数量:"
grep -c "^@" data/pacbio_sample.fastq
echo "Nanopore序列数量:"
grep -c "^@" data/nanopore_sample.fastq
```
#### 2.2 序列长度分析

**执行命令：**
```bash
# 使用Python脚本分析序列长度分布
cd scripts
python analyze_read_length.py
```

**参数解释：**
- `analyze_read_length.py`：分析不同平台的读长分布
- 输出结果包括：平均长度、长度分布图、统计摘要

**预期输出：**
```
Illumina平均读长: 150 bp
PacBio平均读长: 12000 bp
Nanopore平均读长: 8500 bp
```

#### 2.3 质量值分析

**操作说明：**
质量值反映测序的准确性，Phred质量值越高表示错误率越低。

**执行命令：**
```bash
# 分析质量值分布
python analyze_quality.py
```

**结果验证：**
```bash
# 检查生成的质量分析图
ls -la results/quality_*.png
```

**检查点：** 完成序列长度和质量分析，生成相应的统计图表。

---

### 步骤3：平台特征比较分析

#### 3.1 读长分布比较

**操作说明：**
比较不同测序平台的读长分布特征，理解各平台的技术特点。

**执行命令：**
```bash
# 生成读长分布比较图
python compare_platforms.py --analysis read_length
```

**预期结果：**
- Illumina：读长集中在150bp左右，分布较窄
- PacBio：读长分布较宽，平均10-15kb
- Nanopore：读长变异大，可达数十kb

#### 3.2 错误率分析

**执行命令：**
```bash
# 分析不同平台的错误模式
python error_analysis.py
```

**关键指标：**
- **Illumina**：错误率<1%，主要是替换错误
- **PacBio**：错误率5-15%，随机分布
- **Nanopore**：错误率5-10%，插入/缺失为主

#### 3.3 数据产出量比较

**执行命令：**
```bash
# 计算数据产出统计
python throughput_analysis.py
```

**检查点：** 理解不同平台在读长、准确性、产出量方面的差异。---

###
 步骤4：测序成本分析

#### 4.1 成本计算模型

**操作说明：**
建立测序成本计算模型，考虑设备成本、试剂成本、人工成本等因素。

**执行命令：**
```bash
# 运行成本分析脚本
python cost_analysis.py --project_size 30G --coverage 30x
```

**参数解释：**
- `--project_size`：项目数据量需求
- `--coverage`：测序深度要求
- 输出：不同平台的总成本和每Gb成本

#### 4.2 成本效益分析

**执行命令：**
```bash
# 生成成本效益比较报告
python cost_benefit_analysis.py
```

**关键指标：**
- 每Gb成本
- 项目总成本
- 时间成本
- 质量权重成本

**预期输出：**
```
平台成本比较（30x人类基因组）:
Illumina: $1,500 (最经济)
PacBio: $8,000 (高质量)
Nanopore: $3,500 (平衡选择)
```

#### 4.3 投资回报分析

**执行命令：**
```bash
# 分析设备投资回报
python roi_analysis.py
```

**检查点：** 完成成本分析，理解不同平台的经济性考虑。

---

### 步骤5：平台选择决策分析

#### 5.1 需求评估矩阵

**操作说明：**
根据研究需求建立决策矩阵，评估不同平台的适用性。

**执行命令：**
```bash
# 生成平台选择决策矩阵
python platform_decision.py --application genome_assembly
```

**应用场景选项：**
- `genome_assembly`：基因组组装
- `rna_seq`：转录组分析
- `variant_calling`：变异检测
- `metagenomics`：宏基因组

#### 5.2 多因素决策分析

**执行命令：**
```bash
# 运行多因素决策分析
python multi_criteria_decision.py
```

**评估因素：**
- 技术指标（准确性、读长、通量）
- 经济指标（成本、效益）
- 实用指标（时间、难度、可用性）

**结果验证：**
```bash
# 查看决策分析结果
cat results/platform_recommendation.txt
```

**检查点：** 学会根据具体需求选择合适的测序平台。---

### 步
骤6：结果分析和解读

#### 6.1 结果文件说明
| 文件名 | 位置 | 内容说明 | 重要性 |
|-------|------|---------|--------|
| read_length_stats.csv | results/ | 读长统计数据 | 高 |
| quality_analysis.png | results/ | 质量分布图 | 高 |
| platform_comparison.png | results/ | 平台比较图 | 高 |
| cost_analysis.xlsx | results/ | 成本分析表 | 中 |
| decision_matrix.csv | results/ | 决策矩阵 | 中 |

#### 6.2 结果解读
**关键指标：**
- **读长分布**：正常范围因平台而异
- **质量分数**：Q20以上为可接受质量
- **成本效益比**：综合考虑质量和成本
- **适用性评分**：根据应用场景评估

**结果可视化：**
```bash
# 生成综合分析报告
python generate_report.py
```python

## 预期结果

### 主要输出文件
1. **读长分析结果**：`results/read_length_analysis.png`
   - 内容：不同平台读长分布对比图
   - 用途：理解平台技术特征

2. **质量分析结果**：`results/quality_distribution.png`
   - 内容：质量值分布和错误率分析
   - 用途：评估数据质量

3. **成本分析报告**：`results/cost_analysis_report.pdf`
   - 内容：详细的成本效益分析
   - 用途：项目预算和平台选择

### 关键结果指标
- **Illumina读长**：应该在 100-300 bp 之间
- **PacBio读长**：预期值约为 10-15 kb
- **Nanopore读长**：通常 5-50 kb，变异较大
- **质量分数**：Q30以上占比应>80%（Illumina）

### 成功标准
- [ ] 所有脚本执行无错误
- [ ] 生成了预期的分析图表
- [ ] 成本分析结果合理
- [ ] 平台选择建议明确

## 故障排除

### 常见问题1：Python包导入错误
**症状：** ImportError: No module named 'pandas'
**原因：** 缺少必要的Python包
**解决方案：**
```bash
# 安装缺少的包
pip install pandas matplotlib seaborn biopython
```bash

### 常见问题2：数据文件无法下载
**症状：** wget: command not found 或下载失败
**原因：** 网络问题或命令不可用
**解决方案：**
```bash
# 使用备用数据或手动下载
cp /shared/data/lesson01/* data/
```

### 常见问题3：脚本执行权限错误
**症状：** Permission denied
**原因：** 脚本文件没有执行权限
**解决方案：**
```bash
# 添加执行权限
chmod +x scripts/*.py
```python

### 获取帮助
如果遇到其他问题：
1. 检查错误日志：`cat logs/error.log`
2. 查看脚本帮助：`python script_name.py --help`
3. 联系助教或老师：wangys@hunau.edu.cn## 扩展练
习

### 练习1：自定义平台比较
**目标：** 根据特定研究需求定制平台比较分析
**任务：** 
1. 选择一个具体的研究场景（如植物基因组、微生物组等）
2. 调整分析参数，重新运行比较分析
3. 解释结果差异的原因
**提示：** 修改 `compare_platforms.py` 中的参数设置

### 练习2：成本敏感性分析
**目标：** 分析不同因素对测序成本的影响
**任务：** 
1. 改变项目规模（10G, 50G, 100G）
2. 改变测序深度（10x, 30x, 50x）
3. 分析成本变化趋势
**提示：** 使用 `cost_analysis.py` 的不同参数组合

### 练习3：新技术评估
**目标：** 评估新兴测序技术的潜在价值
**任务：** 
1. 研究一种新的测序技术（如蛋白质测序）
2. 预测其技术参数和成本
3. 与现有技术进行比较分析
**提示：** 查阅最新文献和技术报告

### 思考问题
1. 为什么不同测序平台会有不同的错误模式？
2. 在什么情况下会选择成本更高的测序平台？
3. 如何平衡测序成本和数据质量的关系？
4. 未来测序技术发展的主要方向是什么？

## 参考资料

### 相关文献
1. Goodwin S, et al. Coming of age: ten years of next-generation sequencing technologies. Nature Reviews Genetics, 2016.
2. Amarasinghe SL, et al. Opportunities and challenges in long-read sequencing data analysis. Genome Biology, 2020.
3. Logsdon GA, et al. Long-read human genome sequencing and its applications. Nature Reviews Genetics, 2020.

### 在线资源
- Illumina技术文档：https://www.illumina.com/techniques/sequencing.html
- PacBio技术介绍：https://www.pacb.com/smrt-science/
- Oxford Nanopore技术：https://nanoporetech.com/how-it-works
- NCBI测序技术指南：https://www.ncbi.nlm.nih.gov/guide/sequence-analysis/

### 软件文档
- Biopython文档：https://biopython.org/wiki/Documentation
- Pandas用户指南：https://pandas.pydata.org/docs/user_guide/
- Matplotlib教程：https://matplotlib.org/stable/tutorials/index.html

## 附录

### 附录A：完整脚本文件
参见：`scripts/` 目录中的相关脚本文件
- `analyze_read_length.py`：读长分析脚本
- `analyze_quality.py`：质量分析脚本
- `compare_platforms.py`：平台比较脚本
- `cost_analysis.py`：成本分析脚本
- `platform_decision.py`：平台选择脚本

### 附录B：配置文件模板
```python
# config.py - 分析配置文件
PLATFORMS = {
    'illumina': {'max_length': 300, 'avg_quality': 30},
    'pacbio': {'max_length': 100000, 'avg_quality': 20},
    'nanopore': {'max_length': 2000000, 'avg_quality': 15}
}

COST_PARAMS = {
    'illumina': {'per_gb': 25, 'setup': 1000},
    'pacbio': {'per_gb': 200, 'setup': 5000},
    'nanopore': {'per_gb': 100, 'setup': 2000}
}
```

### 附录C：数据格式说明
**FASTQ格式详解：**
```
@序列标识符 [可选描述]
ATCGATCGATCG...  (DNA序列)
+[可选重复标识符]
!"#$%&'()*+,...  (质量值，ASCII编码)
```

**质量值编码：**
- Phred+33编码（标准）
- ASCII字符对应质量分数
- 质量分数 = ASCII值 - 33

---

**实验完成时间：** 预计 2 小时  
**难度等级：** 初级  
**最后更新：** 2025年
