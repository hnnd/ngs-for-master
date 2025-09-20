# R语言基础课程

## 课程信息
- **课程名称**：高通量测序数据分析 - R语言基础
- **主讲教师**：王运生教授
- **联系邮箱**：wangys@hunau.edu.cn
- **办公室**：16教420室
- **上课地点**：105机房

## 课程目标
通过本次课程，学生将能够：
1. 掌握R语言基础语法和数据结构
2. 熟练使用R进行数据处理和统计分析
3. 创建各种类型的数据可视化图表
4. 了解R在生物信息学中的应用
5. 使用R包进行生物数据分析

## 课程内容
- **理论课时**：1小时
- **实践课时**：3小时
- **幻灯片数量**：约25张

### 主要内容
1. R语言简介与环境配置
2. 基础语法与数据类型
3. 数据结构详解
4. 数据处理与统计分析
5. 数据可视化
6. 生物信息学R包介绍

## 文件结构
```
r-basics/
├── README.md                 # 课程说明
├── slides/                   # 幻灯片文件
│   ├── slides.md            # Marp源文件
│   └── slides.pdf           # PDF版本
├── manual/                   # 实践手册
│   ├── manual.md            # 手册主文件
│   └── scripts/             # 相关脚本
│       ├── setup.R          # 环境设置
│       ├── basic_exercises.R # 基础练习
│       ├── data_analysis.R   # 数据分析
│       ├── visualization.R   # 可视化
│       └── solutions.R      # 练习答案
└── data/                    # 示例数据
    ├── gene_expression.csv  # 基因表达数据
    ├── sample_info.txt      # 样本信息
    └── sequences.fasta      # 序列数据
```

## 环境要求
- R 4.0或更高版本
- RStudio（推荐）
- 必需R包：tidyverse, ggplot2, dplyr, readr, BiocManager

## 学习资源
- [R官方网站](https://www.r-project.org/)
- [RStudio文档](https://rstudio.com/resources/)
- [Bioconductor项目](https://bioconductor.org/)