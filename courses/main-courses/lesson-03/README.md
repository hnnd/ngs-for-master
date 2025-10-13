# 第3次课：高通量测序序列比对算法与工具

## 课程概述

本次课程深入讲解高通量测序数据分析中最核心的步骤之一：序列比对。学生将学习BWA和Bowtie2两种主流比对工具的原理、使用方法和参数优化。

## 学习目标

1. **理解比对算法原理**
   - Burrows-Wheeler Transform (BWT)
   - FM-index索引结构
   - 种子延伸策略

2. **掌握比对工具使用**
   - BWA-MEM的使用和参数优化
   - Bowtie2的不同比对模式
   - SAM/BAM格式处理

3. **比对质量评估**
   - 比对统计分析
   - 质量控制标准
   - 结果可视化

## 课程材料

### 📊 课件 (slides/)
- \`slides.md\` - Marp格式的课程幻灯片
- 包含算法原理、工具介绍和实例演示

### 📖 实验手册 (manual/)
- \`manual.md\` - 详细的实验操作指南
- \`DATA_SOURCES.md\` - 数据来源和准备说明
- \`scripts/\` - 实验脚本目录
  - \`prepare_data.sh\` - 自动化数据准备脚本

### 📁 目录结构
\`\`\`
lesson-03/
├── README.md                    # 本文件
├── slides/
│   └── slides.md               # 课程幻灯片
├── manual/
│   ├── manual.md               # 实验手册
│   ├── DATA_SOURCES.md         # 数据说明文档
│   ├── data/                   # 测序数据目录
│   ├── reference/              # 参考基因组目录
│   ├── results/                # 结果输出目录
│   ├── logs/                   # 日志文件目录
│   └── scripts/                # 分析脚本
│       ├── prepare_data.sh     # 数据准备脚本
│       ├── alignment_stats.py  # 比对统计脚本
│       ├── visualize_results.py # 结果可视化脚本
│       └── generate_report.py  # 报告生成脚本
└── images/                     # 课程图片资源
\`\`\`

## 快速开始

### 1. 准备工作环境

\`\`\`bash
# 进入实验目录
cd lesson-03/manual

# 创建必要的子目录
mkdir -p data reference results logs scripts

# 检查软件安装
bwa 2>&1 | head -3
bowtie2 --version | head -1
samtools --version | head -1
\`\`\`

### 2. 准备实验数据

**推荐方式（使用自动化脚本）：**
\`\`\`bash
# 给脚本添加执行权限
chmod +x scripts/prepare_data.sh

# 运行数据准备脚本
bash scripts/prepare_data.sh

# 该脚本将自动：
# 1. 从Ensembl/UCSC下载人类22号染色体参考序列
# 2. 从NCBI SRA下载真实测序数据或生成模拟数据
# 3. 验证数据完整性
\`\`\`

**手动方式：**
详细步骤请参考 [manual/manual.md](manual/manual.md) 第1.4节

### 3. 开始实验

按照 [manual/manual.md](manual/manual.md) 中的步骤依次完成：

1. ✅ 环境设置和数据准备
2. 🔧 构建参考基因组索引
3. 🧬 使用BWA进行序列比对
4. 🎯 使用Bowtie2进行序列比对
5. 📊 比对结果统计和分析
6. 📈 结果可视化
7. ⚡ 性能比较和分析

## 软件要求

### 核心工具
- **BWA** (>= 0.7.17) - 序列比对工具
- **Bowtie2** (>= 2.4.0) - 序列比对工具
- **samtools** (>= 1.10) - SAM/BAM文件处理

### 辅助工具
- **Python 3.8+** - 脚本运行环境
  - pandas - 数据处理
  - matplotlib - 图表绘制
  - seaborn - 高级可视化
- **IGV** - 基因组浏览器（可选）

### 安装方法

**使用Conda（推荐）：**
\`\`\`bash
# 安装比对工具
conda install -c bioconda bwa bowtie2 samtools

# 安装数据下载工具（可选）
conda install -c bioconda sra-tools wgsim

# 安装Python依赖
pip install pandas matplotlib seaborn
\`\`\`

## 数据说明

### 参考基因组
- **来源**: Ensembl或UCSC公共数据库
- **版本**: GRCh38 (hg38)
- **内容**: 人类22号染色体前5Mb
- **用途**: 序列比对的参考序列

### 测序数据
提供三种获取方式：

1. **真实数据** (推荐用于完整体验)
   - 来源: NCBI SRA (SRR622461)
   - 类型: Illumina paired-end测序
   - 规模: 100万条reads

2. **模拟数据** (推荐用于快速实验)
   - 工具: wgsim
   - 优点: 快速生成，参数可控

3. **课程服务器** (如果可用)
   - 优点: 下载速度快
   - 需要: 内网访问权限

详细信息请参阅 [manual/DATA_SOURCES.md](manual/DATA_SOURCES.md)

## 预期结果

完成本次实验后，学生应该能够：

### 技能掌握
- ✅ 构建参考基因组索引
- ✅ 使用BWA和Bowtie2进行序列比对
- ✅ 处理和转换SAM/BAM文件
- ✅ 评估比对质量
- ✅ 优化比对参数

## 扩展资源

### 官方文档
- [BWA Manual](http://bio-bwa.sourceforge.net/bwa.shtml)
- [Bowtie2 Manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
- [SAMtools Documentation](http://www.htslib.org/doc/)

### 重要文献
1. Li & Durbin (2009). Fast and accurate short read alignment with Burrows-Wheeler transform. *Bioinformatics*.
2. Langmead & Salzberg (2012). Fast gapped-read alignment with Bowtie 2. *Nature Methods*.

## 故障排除

### 常见问题

**Q: 内存不足错误**
\`\`\`bash
# 减少线程数
bwa mem -t 2 ...
\`\`\`

**Q: 无法下载数据**
\`\`\`bash
# 使用模拟数据
bash scripts/prepare_data.sh --sequencing
\`\`\`

更多问题请参考 [manual/manual.md](manual/manual.md) 的故障排除章节。

## 联系方式

- **主讲教师**: 王运生
- **邮箱**: wangys@hunau.edu.cn
- **上课地点**: 105机房

---

**更新日期**: 2025-10-13
