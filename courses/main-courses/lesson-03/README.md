# 第3次课：高通量测序序列比对算法与工具

## 课程信息
- **课程名称**：高通量测序数据分析
- **主讲教师**：王运生教授
- **联系邮箱**：wangys@hunau.edu.cn
- **办公室**：16教420室
- **上课地点**：105机房
- **课程时长**：4学时（理论2学时 + 实践2学时）

## 学习目标

通过本次课程，学生将能够：

1. **理解比对算法原理**
   - 掌握序列比对的基本概念和算法分类
   - 理解BWA、Bowtie2等主流比对工具的算法原理
   - 了解不同算法的适用场景和性能特点

2. **掌握比对工具使用**
   - 学会使用BWA进行短读长序列比对
   - 掌握Bowtie2的参数设置和优化策略
   - 了解其他主流比对工具的特点和使用方法

3. **进行比对结果分析**
   - 理解SAM/BAM文件格式和内容
   - 掌握比对质量评估方法
   - 学会使用samtools进行结果处理和统计

4. **优化比对策略**
   - 了解比对参数对结果的影响
   - 掌握针对不同数据类型的比对策略选择
   - 学会比对结果的可视化和质量控制

## 课程内容

### 理论部分（2学时）
1. 序列比对算法基础
2. 主流比对工具原理与比较
3. 比对策略与参数优化
4. 比对质量评估方法

### 实践部分（2学时）
1. 参考基因组准备和索引构建
2. BWA和Bowtie2比对实战
3. SAM/BAM文件处理和统计
4. 比对结果可视化和质量评估

## 文件结构

```
lesson-03/
├── README.md                    # 本文件
├── slides/
│   └── slides.md               # 理论课幻灯片
├── manual/
│   ├── manual.md               # 实践操作手册
│   ├── scripts/                # 分析脚本
│   │   ├── setup.sh           # 环境设置
│   │   ├── build_index.sh     # 索引构建
│   │   ├── run_bwa.sh         # BWA比对
│   │   ├── run_bowtie2.sh     # Bowtie2比对
│   │   ├── process_sam.sh     # SAM处理
│   │   ├── alignment_stats.py  # 比对统计
│   │   └── visualize_results.py # 结果可视化
│   └── data/                   # 示例数据
│       ├── reference.fa        # 参考基因组片段
│       ├── sample_R1.fastq     # 测序数据1
│       └── sample_R2.fastq     # 测序数据2
└── images/                     # 课程专用图片
    ├── alignment_algorithms.svg # 比对算法原理图
    ├── tool_comparison.svg     # 工具比较图
    └── sam_format.svg          # SAM格式说明图
```

## 前置要求

### 软件环境
- BWA (>= 0.7.17)
- Bowtie2 (>= 2.4.0)
- samtools (>= 1.10)
- Python 3.8+ (pandas, matplotlib, seaborn)
- IGV (Integrative Genomics Viewer)

### 知识基础
- 熟悉Linux命令行操作
- 了解FASTQ文件格式
- 掌握基本的生物信息学概念

## 学习成果

完成本次课程后，学生应该能够：
- 独立进行高通量测序数据的比对分析
- 选择合适的比对工具和参数
- 评估和优化比对结果质量
- 处理和分析SAM/BAM格式文件

## 参考资料

1. Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics, 25, 1754-1760.
2. Langmead B. and Salzberg S. (2012) Fast gapped-read alignment with Bowtie 2. Nature Methods, 9, 357-359.
3. Li H. et al. (2009) The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25, 2078-2079.
4. BWA Manual: http://bio-bwa.sourceforge.net/bwa.shtml
5. Bowtie2 Manual: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml