# 高通量测序数据分析课程教学材料目录结构规范

## 整体目录结构

```
ngs-course-materials/
├── templates/                    # 模板文件
│   ├── slide-template.md        # Marp幻灯片模板
│   ├── manual-template.md       # 实践手册模板
│   └── styles/                  # 样式文件
│       ├── course-theme.css     # 课程主题样式
│       └── code-highlight.css   # 代码高亮样式
├── assets/                      # 资源文件
│   ├── images/                  # SVG图片资源
│   │   ├── common/             # 通用图片
│   │   ├── sequencing/         # 测序相关图片
│   │   ├── analysis/           # 分析流程图片
│   │   └── results/            # 结果展示图片
│   ├── logos/                   # 学校/课程标识
│   └── icons/                   # 图标资源
├── courses/                     # 课程内容
│   ├── prep-courses/           # 预备课程
│   │   ├── linux-basics/       # Linux基础
│   │   ├── python-basics/      # Python基础
│   │   └── r-basics/           # R语言基础
│   └── main-courses/           # 正式课程
│       ├── lesson-01/          # 第1次课
│       ├── lesson-02/          # 第2次课
│       ├── lesson-03/          # 第3次课
│       ├── lesson-04/          # 第4次课
│       ├── lesson-05/          # 第5次课
│       ├── lesson-06/          # 第6次课
│       ├── lesson-07/          # 第7次课
│       └── lesson-08/          # 第8次课
├── scripts/                    # 脚本文件
│   ├── generators/             # 生成器脚本
│   ├── validators/             # 验证脚本
│   └── utilities/              # 工具脚本
└── docs/                       # 文档
    ├── style-guide.md          # 样式指南
    ├── template-usage.md       # 模板使用说明
    ├── directory-structure.md  # 本文档
    └── maintenance.md          # 维护指南
```

## 课程目录详细结构

### 预备课程目录结构
```
prep-courses/[course-name]/
├── slides/                     # 幻灯片文件
│   ├── slides.md              # Marp源文件
│   ├── slides.html            # 生成的HTML文件
│   └── slides.pdf             # 导出的PDF文件
├── manual/                    # 实践手册
│   ├── manual.md              # 手册主文件
│   ├── scripts/               # 相关脚本
│   │   ├── setup.sh          # 环境设置脚本
│   │   ├── exercises.py      # 练习脚本
│   │   └── solutions/        # 练习答案
│   └── data/                  # 示例数据
│       ├── sample-files/     # 样本文件
│       └── expected-outputs/ # 预期输出
├── images/                    # 课程专用图片
│   ├── concepts/             # 概念图解
│   ├── screenshots/          # 操作截图
│   └── examples/             # 示例图表
└── README.md                  # 课程说明文档
```

### 正式课程目录结构
```
main-courses/lesson-[XX]/
├── slides/                     # 幻灯片文件
│   ├── slides.md              # Marp源文件
│   ├── slides.html            # 生成的HTML文件
│   └── slides.pdf             # 导出的PDF文件
├── manual/                    # 实践手册
│   ├── manual.md              # 手册主文件
│   ├── scripts/               # 相关脚本
│   │   ├── setup.sh          # 环境设置脚本
│   │   ├── download_data.sh  # 数据下载脚本
│   │   ├── analysis.py       # 分析脚本
│   │   ├── analysis.R        # R分析脚本
│   │   ├── visualization.py  # 可视化脚本
│   │   └── cleanup.sh        # 清理脚本
│   └── data/                  # 示例数据
│       ├── raw-data/         # 原始数据
│       ├── processed-data/   # 处理后数据
│       └── expected-results/ # 预期结果
├── images/                    # 课程专用图片
│   ├── diagrams/             # 原理图表
│   ├── workflows/            # 流程图
│   ├── screenshots/          # 操作截图
│   └── results/              # 结果图表
└── README.md                  # 课程说明文档
```

## 文件命名规范

### 基本命名原则
1. **使用英文命名**：所有文件和目录名使用英文，避免中文字符
2. **小写字母**：使用小写字母，单词间用连字符(-)分隔
3. **描述性命名**：文件名应清楚描述文件内容和用途
4. **版本控制**：如需版本标识，使用v1.0格式

### 目录命名规范
- **课程目录**：`lesson-[XX]` (XX为两位数字，如lesson-01)
- **预备课程**：`[subject]-basics` (如linux-basics, python-basics)
- **功能目录**：使用复数形式 (如slides, scripts, images)

### 文件命名规范

#### 幻灯片文件
- **源文件**：`slides.md`
- **HTML文件**：`slides.html`
- **PDF文件**：`slides.pdf`
- **特殊主题**：`slides-[topic].md` (如slides-introduction.md)

#### 手册文件
- **主手册**：`manual.md`
- **分章节**：`manual-[section].md` (如manual-setup.md)

#### 脚本文件
- **环境设置**：`setup.sh`
- **数据下载**：`download_data.sh`
- **数据清理**：`cleanup.sh`
- **分析脚本**：`analysis_[type].[ext]` (如analysis_qc.py)
- **可视化脚本**：`visualization_[type].[ext]` (如visualization_plots.R)

#### 图片文件
- **SVG图片**：`[description].svg` (如sequencing-workflow.svg)
- **PNG图片**：`[description].png` (如quality-distribution.png)
- **截图文件**：`screenshot-[step].png` (如screenshot-fastqc.png)

#### 数据文件
- **原始数据**：保持原始文件名或使用描述性名称
- **处理数据**：`[type]_processed.[ext]` (如reads_processed.fastq)
- **结果文件**：`[analysis]_results.[ext]` (如alignment_results.sam)

### 特殊文件命名

#### 配置文件
- **Marp配置**：`marp.config.js`
- **样式文件**：`course-theme.css`, `code-highlight.css`
- **环境配置**：`environment.yml`, `requirements.txt`

#### 文档文件
- **说明文档**：`README.md`
- **使用指南**：`usage-guide.md`
- **故障排除**：`troubleshooting.md`

## 版本控制规范

### Git 提交规范
- **功能添加**：`feat: 添加第X次课幻灯片`
- **内容修复**：`fix: 修正Python脚本语法错误`
- **文档更新**：`docs: 更新使用说明`
- **样式调整**：`style: 优化幻灯片样式`

### 版本标签规范
- **主版本**：`v1.0.0` (完整课程材料第一版)
- **次版本**：`v1.1.0` (添加新课程或重大更新)
- **修订版本**：`v1.1.1` (错误修复和小幅改进)

## 文件大小限制

### 建议文件大小
- **Markdown文件**：< 1MB
- **SVG图片**：< 500KB
- **PNG图片**：< 2MB
- **PDF文件**：< 10MB
- **数据文件**：< 100MB (大文件使用外部链接)

### 大文件处理
- 使用Git LFS管理大型数据文件
- 提供数据下载脚本而非直接存储
- 压缩图片和PDF文件以减小体积

## 权限和访问控制

### 文件权限
- **可执行脚本**：755 (rwxr-xr-x)
- **数据文件**：644 (rw-r--r--)
- **配置文件**：644 (rw-r--r--)

### 访问控制
- **公开文件**：README.md, 使用说明等
- **受限文件**：答案脚本, 评分标准等
- **私有文件**：个人笔记, 草稿等

## 备份和同步

### 备份策略
- **每日备份**：自动备份到云存储
- **版本备份**：每个版本发布后创建完整备份
- **增量备份**：定期同步修改的文件

### 同步规范
- **主分支**：stable, 用于正式发布
- **开发分支**：develop, 用于日常开发
- **功能分支**：feature/[name], 用于新功能开发

这个目录结构规范确保了教学材料的组织性、可维护性和可扩展性，便于团队协作和长期维护。