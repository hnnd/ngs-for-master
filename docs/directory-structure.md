# 高通量测序数据分析课程教学材料目录结构说明

## 项目整体结构

```
ngs-course-materials/
├── README.md                    # 项目主说明文档
├── package.json                 # Node.js依赖配置
├── requirements.txt             # Python依赖配置
├── .gitignore                   # Git忽略文件配置
├── .gitattributes              # Git属性配置
│
├── templates/                   # 模板文件目录（本地开发用，不推送到远程）
│   ├── slide-template.md        # Marp幻灯片模板
│   ├── manual-template.md       # 实践手册模板
│   └── styles/                  # 样式文件
│       ├── course-theme.css     # 课程主题样式
│       ├── ngs-course.css       # NGS课程专用样式
│       └── code-highlight.css   # 代码高亮样式
│
├── assets/                      # 全局资源文件
│   ├── images/                  # 通用图片资源
│   │   ├── common/              # 通用图表和图标
│   │   ├── sequencing/          # 测序技术相关图片
│   │   ├── analysis/            # 分析流程图片
│   │   └── results/             # 结果展示图片
│   ├── logos/                   # 标识和Logo
│   │   ├── logo.png             # 课程Logo
│   │   └── README.md            # Logo使用说明
│   └── icons/                   # 图标资源
│       └── README.md            # 图标使用说明
│
├── courses/                     # 课程内容目录
│   ├── README.md                # 课程总体说明
│   ├── prep-courses/            # 预备课程
│   │   ├── linux-basics/        # Linux基础课程
│   │   ├── python-basics/       # Python基础课程
│   │   └── r-basics/            # R语言基础课程
│   └── main-courses/            # 正式课程
│       ├── lesson-01/           # 第1次课：测序技术原理
│       ├── lesson-02/           # 第2次课：质量控制
│       ├── lesson-03/           # 第3次课：序列比对
│       ├── lesson-04/           # 第4次课：变异检测
│       ├── lesson-05/           # 第5次课：转录组分析
│       ├── lesson-06/           # 第6次课：ChIP-seq分析
│       ├── lesson-07/           # 第7次课：单细胞分析
│       └── lesson-08/           # 第8次课：多组学整合
│
├── scripts/                     # 脚本工具目录
│   ├── README.md                # 脚本使用说明
│   ├── setup-repository.sh      # 仓库初始化脚本
│   ├── generators/              # 生成器脚本
│   │   └── .gitkeep
│   ├── validators/              # 验证脚本
│   │   └── .gitkeep
│   └── utilities/               # 工具脚本
│       └── .gitkeep
│
└── docs/                        # 文档目录
    ├── directory-structure.md   # 目录结构说明（本文档）
    ├── style-guide.md           # 样式指南
    ├── template-usage.md        # 模板使用说明
    └── maintenance.md           # 维护指南
```

## 单个课程目录结构

每个课程目录（如 `lesson-01/`）都遵循以下标准结构：

```
lesson-xx/
├── README.md                    # 课程概述和学习目标
│
├── slides/                      # 幻灯片文件
│   ├── slides.md                # Marp源文件
│   ├── slides.html              # 生成的HTML文件（可选）
│   └── slides.pdf               # 导出的PDF文件（可选）
│
├── manual/                      # 实践手册
│   ├── manual.md                # 手册主文件
│   ├── scripts/                 # 相关脚本
│   │   ├── setup.sh             # 环境设置脚本
│   │   ├── download_data.sh     # 数据下载脚本
│   │   ├── analysis_*.py        # 分析脚本（Python）
│   │   ├── analysis_*.R         # 分析脚本（R语言）
│   │   ├── visualization_*.py   # 可视化脚本
│   │   └── cleanup.sh           # 清理脚本
│   └── data/                    # 示例数据
│       ├── sample_data/         # 样本数据文件
│       ├── expected_results/    # 预期结果文件
│       └── README.md            # 数据说明
│
└── images/                      # 课程专用图片
    ├── diagrams/                # 原理图表
    ├── screenshots/             # 操作截图
    ├── results/                 # 结果图表
    └── workflows/               # 流程图
```

## 文件命名规范

### 1. 目录命名规范

- **课程目录**：使用 `lesson-XX` 格式，XX为两位数字（如 `lesson-01`）
- **预备课程**：使用描述性名称加连字符（如 `linux-basics`）
- **功能目录**：使用英文小写加连字符（如 `main-courses`）
- **资源目录**：使用复数形式（如 `images`, `scripts`）

### 2. 文件命名规范

#### 幻灯片文件
- **源文件**：`slides.md`
- **HTML文件**：`slides.html`
- **PDF文件**：`slides.pdf`

#### 手册文件
- **主文件**：`manual.md`
- **脚本文件**：使用功能描述加扩展名
  - `setup.sh` - 环境设置
  - `download_data.sh` - 数据下载
  - `run_analysis.py` - 运行分析
  - `generate_plots.R` - 生成图表
  - `cleanup.sh` - 清理文件

#### 图片文件
- **SVG图片**：使用描述性名称，小写加连字符
  - `sequencing-workflow.svg`
  - `quality-distribution.svg`
  - `alignment-statistics.svg`
- **PNG/JPG图片**：用于截图和照片
  - `fastqc-report-screenshot.png`
  - `igv-visualization.jpg`

#### 数据文件
- **示例数据**：使用描述性前缀
  - `sample_reads.fastq`
  - `reference_genome.fasta`
  - `gene_expression.csv`
- **结果文件**：使用 `expected_` 前缀
  - `expected_alignment.sam`
  - `expected_variants.vcf`

### 3. 版本控制规范

#### Git提交信息格式
```
类型(范围): 简短描述

详细描述（可选）

相关问题编号（可选）
```

**类型标识：**
- `feat`: 新功能
- `fix`: 修复问题
- `docs`: 文档更新
- `style`: 格式调整
- `refactor`: 代码重构
- `test`: 测试相关
- `chore`: 构建过程或辅助工具的变动

**示例：**
```
feat(lesson-01): 添加测序技术原理幻灯片

- 完成Illumina、PacBio、Nanopore平台介绍
- 添加技术对比表格和流程图
- 创建相关SVG图片资源

Closes #123
```

## 内容组织原则

### 1. 模块化设计
- 每个课程独立成模块，便于单独使用和维护
- 共享资源统一管理，避免重复
- 脚本和数据分离，提高复用性

### 2. 层次化结构
- 按照课程 → 类型 → 具体文件的层次组织
- 相同类型的文件放在同一目录下
- 使用README文件说明每个目录的用途

### 3. 标准化命名
- 统一的命名规范，便于查找和管理
- 描述性的文件名，见名知意
- 版本信息和时间戳的合理使用

### 4. 文档完整性
- 每个重要目录都有README说明
- 关键文件包含头部注释
- 使用说明和示例齐全

## 存储和备份策略

### 1. 版本控制
- 所有源文件纳入Git版本控制
- 生成的文件（HTML、PDF）可选择性提交
- 大型数据文件使用Git LFS管理

### 2. 文件分类存储
- **源文件**：Markdown、脚本、配置文件
- **生成文件**：HTML、PDF、图片
- **数据文件**：示例数据、结果文件
- **临时文件**：构建过程中的中间文件
- **模板文件**：仅本地使用，不推送到远程仓库

### 3. Templates目录特殊说明
`templates/` 目录包含开发用的模板文件和样式，具有以下特点：
- **本地开发**：仅在本地开发环境中使用
- **不推送远程**：通过.gitignore排除，不会推送到GitHub
- **个人定制**：开发者可以根据需要自定义模板
- **版本独立**：不受远程仓库版本控制影响

如需获取模板文件，请：
1. 从项目维护者处获取模板文件
2. 或参考现有课程内容创建自己的模板
3. 将模板文件放置在本地的 `templates/` 目录中

### 4. 备份机制
- 定期备份到远程仓库
- 重要节点创建标签（tag）
- 关键文件的本地备份

## 维护和更新

### 1. 定期检查
- 检查链接有效性
- 验证脚本可执行性
- 更新软件版本信息

### 2. 内容更新
- 跟踪技术发展，更新内容
- 根据教学反馈调整材料
- 修复发现的问题和错误

### 3. 结构优化
- 根据使用情况调整目录结构
- 优化文件组织和命名
- 改进文档和说明

## 使用建议

### 1. 新用户
- 先阅读主README文档
- 查看template-usage.md了解模板使用
- 从预备课程开始学习

### 2. 教师用户
- 根据课程安排选择相应材料
- 可以自定义调整内容和进度
- 使用提供的脚本工具提高效率

### 3. 开发者
- 遵循既定的目录结构和命名规范
- 使用模板创建新内容
- 提交前进行完整测试

这个目录结构设计确保了项目的可维护性、可扩展性和易用性，为高质量的教学材料提供了坚实的组织基础。