# 教学材料模板使用指南

## 概述

本指南详细说明如何使用高通量测序数据分析课程的教学材料模板，包括幻灯片模板和实践手册模板的使用方法。

## 幻灯片模板使用

### 基本设置

#### 1. 文件头部配置
```yaml
---
marp: true
theme: ngs-course
paginate: true
header: '高通量测序数据分析 | 王运生教授'
footer: 'wangys@hunau.edu.cn | 16教420室'
---
```

**配置说明：**
- `marp: true`：启用Marp功能
- `theme: ngs-course`：使用课程专用主题
- `paginate: true`：显示页码
- `header`：页眉内容（课程名称和讲师）
- `footer`：页脚内容（联系方式和办公室）

#### 2. 课程信息注释
```html
<!-- 
课程幻灯片模板
课程名称：高通量测序数据分析
主讲教师：王运生教授
联系邮箱：wangys@hunau.edu.cn
办公室：16教420室
上课地点：105机房
-->
```

### 页面类型使用

#### 1. 标题页 (title)
```markdown
<!-- _class: title -->
# 测序技术原理与平台比较
## 高通量测序数据分析

**主讲教师：** 王运生教授  
**联系邮箱：** wangys@hunau.edu.cn  
**办公室：** 16教420室  
**上课地点：** 105机房  
```

**使用场景：** 每次课程的开始页面

#### 2. 目录页 (toc)
```markdown
<!-- _class: toc -->
# 本次课程内容

1. 测序技术发展历史
2. 主流测序平台原理
3. 平台特征比较分析
4. 实践操作演示

**学习目标：**
- 理解不同测序技术的基本原理
- 掌握平台选择的考虑因素
- 学会分析测序数据特征
```

**使用场景：** 课程内容概览和学习目标展示

#### 3. 内容页 (content)
```markdown
<!-- _class: content -->
# Illumina测序技术

## 基本原理

- 桥式PCR扩增
- 可逆终止子技术
- 荧光标记检测
- 边合成边测序

![Illumina原理图](../assets/images/illumina-principle.svg)
```

**使用场景：** 标准内容展示，适合文字和图片结合

#### 4. 多列页 (multi-column)
```markdown
<!-- _class: multi-column -->
# 测序平台比较

<div class="columns">
<div class="column">

## Illumina平台
- 读长：150-300bp
- 通量：高
- 准确率：>99%
- 成本：低

</div>
<div class="column">

## PacBio平台
- 读长：10-20kb
- 通量：中等
- 准确率：>99%
- 成本：高

</div>
</div>
```

**使用场景：** 对比内容、并列信息展示

#### 5. 图片页 (image)
```markdown
<!-- _class: image -->
# 测序数据质量分布

![质量分布图](../assets/images/quality-distribution.svg)
```

**使用场景：** 大图展示，最小化文字干扰

#### 6. 代码页 (code)
```markdown
<!-- _class: code -->
# FastQC质量控制

```bash
# 安装FastQC
conda install -c bioconda fastqc

# 运行质量控制
fastqc sample.fastq -o output_dir
```

```python
# Python处理FASTQ文件
from Bio import SeqIO
for record in SeqIO.parse("sample.fastq", "fastq"):
    print(f"ID: {record.id}, Length: {len(record)}")
```
```

**使用场景：** 代码展示和命令行操作

#### 7. 总结页 (summary)
```markdown
<!-- _class: summary -->
# 本次课程总结

## 主要内容回顾
- 测序技术发展历程
- 三大主流平台特点
- 平台选择考虑因素

## 下次课程预告
- 测序数据质量控制
- FastQC工具使用
- 数据预处理流程

**作业/练习：**
- 比较不同平台的技术参数
- 阅读相关文献资料
```

**使用场景：** 课程结束时的总结和预告

#### 8. 结束页 (end)
```markdown
<!-- _class: end -->
# 谢谢大家！

**有问题请联系：**
- 邮箱：wangys@hunau.edu.cn
- 办公室：16教420室
- 答疑时间：周三下午2-4点
```

**使用场景：** 课程结束页面

### 图片插入规范

#### SVG图片插入
```markdown
![图片描述](../assets/images/diagram-name.svg)
```

#### 图片尺寸控制
```html
<img src="../assets/images/large-diagram.svg" width="80%" alt="大型图表">
```

#### 图片居中
```html
<div style="text-align: center;">
<img src="../assets/images/centered-image.svg" alt="居中图片">
</div>
```

## 实践手册模板使用

### 基本结构

#### 1. 课程信息部分
```markdown
## 课程信息
- **课程名称**：高通量测序数据分析
- **主讲教师**：王运生
- **联系邮箱**：wangys@hunau.edu.cn
- **办公室**：16教420室
- **上课地点**：105机房
- **课程时间**：第X周 周X 第X-X节
```

**必填字段：** 所有字段都必须填写，确保信息完整

#### 2. 实验目标设定
```markdown
## 实验目标

### 主要目标
- 掌握FastQC质量控制工具的使用
- 理解测序数据质量评估指标
- 学会解读质量控制报告

### 预期成果
- 能够独立运行FastQC分析
- 能够识别常见的质量问题
- 能够制定数据清洗策略
```

**编写要点：** 目标要具体、可测量、可达成

#### 3. 环境要求表格
```markdown
### 软件环境
| 软件名称 | 版本要求 | 安装方式 | 说明 |
|---------|---------|---------|------|
| FastQC | ≥0.11.9 | `conda install -c bioconda fastqc` | 质量控制工具 |
| MultiQC | ≥1.9 | `pip install multiqc` | 报告整合工具 |
| Python | ≥3.7 | 系统自带或conda安装 | 脚本运行环境 |
```

**表格要求：** 信息完整，命令准确，说明清晰

### 操作步骤编写

#### 1. 步骤标题规范
```markdown
### 步骤1：环境设置和准备工作
#### 1.1 创建工作目录
#### 1.2 检查软件环境
#### 1.3 下载和准备数据
```

**编号规则：** 主步骤用数字，子步骤用小数点编号

#### 2. 命令块格式
```markdown
```bash
# 创建工作目录
mkdir -p ~/ngs-analysis/lesson-02
cd ~/ngs-analysis/lesson-02

# 创建子目录
mkdir -p {data,scripts,results,logs}
```
```

**格式要求：** 
- 使用正确的语言标识符
- 添加注释说明命令用途
- 保持缩进和格式一致

#### 3. 预期输出展示
```markdown
**预期输出：**
```
FastQC v0.11.9
Copyright (c) 2010-2020 Babraham Bioinformatics
Analysis complete for sample.fastq
```
```

**展示原则：** 显示关键输出信息，帮助学生验证结果

#### 4. 检查点设置
```markdown
**检查点：** 确认所有数据文件已正确下载并位于 `data/` 目录中。

验证方法：
```bash
ls -la data/
# 应该看到以下文件：
# sample_R1.fastq.gz
# sample_R2.fastq.gz
```
```

**设置原则：** 在关键步骤后设置检查点，确保学生跟上进度

### 故障排除编写

#### 1. 问题描述格式
```markdown
### 常见问题1：FastQC无法运行
**症状：** 运行fastqc命令时提示"command not found"
**原因：** FastQC未正确安装或未添加到PATH
**解决方案：**
```bash
# 重新安装FastQC
conda install -c bioconda fastqc

# 检查安装是否成功
which fastqc
fastqc --version
```
```

**编写要点：** 症状描述准确，原因分析合理，解决方案可操作

#### 2. 分级处理
- **常见问题**：学生经常遇到的问题
- **技术问题**：软件或系统相关问题
- **数据问题**：数据文件或格式问题

### 扩展练习设计

#### 1. 练习题格式
```markdown
### 练习1：批量质量控制
**目标：** 学会对多个文件进行批量质量控制
**任务：** 编写脚本对data目录下所有FASTQ文件运行FastQC
**提示：** 使用for循环和通配符

**参考代码：**
```bash
for file in data/*.fastq; do
    fastqc "$file" -o results/
done
```
```

**设计原则：** 循序渐进，从简单到复杂，提供适当提示

## 模板定制指南

### 样式定制

#### 1. 颜色主题修改
在 `course-theme.css` 中修改颜色变量：
```css
:root {
  --primary-color: #1e3a8a;      /* 主色调 */
  --secondary-color: #059669;     /* 辅助色 */
  --accent-color: #dc2626;        /* 强调色 */
}
```

#### 2. 字体设置
```css
:root {
  --main-font: 'Source Han Sans CN', 'Microsoft YaHei', sans-serif;
  --code-font: 'JetBrains Mono', 'Consolas', monospace;
}
```

### 内容定制

#### 1. 课程信息更新
在模板中全局替换以下信息：
- 讲师姓名：王运生教授
- 联系邮箱：wangys@hunau.edu.cn
- 办公室：16教420室
- 上课地点：105机房

#### 2. 学校标识添加
```markdown
![学校Logo](../assets/logos/university-logo.svg)
```

## 质量检查清单

### 幻灯片检查
- [ ] 页眉页脚信息正确
- [ ] 所有图片链接有效
- [ ] 代码语法高亮正常
- [ ] 页面布局美观统一
- [ ] 中文字符显示正常

### 手册检查
- [ ] 课程信息完整准确
- [ ] 操作步骤清晰可行
- [ ] 命令语法正确无误
- [ ] 预期输出真实有效
- [ ] 故障排除覆盖全面

### 文件检查
- [ ] 文件命名符合规范
- [ ] 目录结构完整正确
- [ ] 版本信息及时更新
- [ ] 备份文件已创建

## 常见问题解答

### Q1: 如何添加新的页面样式？
A: 在 `course-theme.css` 中添加新的class定义，然后在Markdown中使用 `<!-- _class: new-style -->` 应用。

### Q2: 图片显示不正常怎么办？
A: 检查图片路径是否正确，确保SVG文件中的中文字符已正确编码。

### Q3: 代码高亮不生效？
A: 确认代码块使用了正确的语言标识符，如 ```bash、```python、```r。

### Q4: 如何批量生成PDF文件？
A: 使用Marp CLI命令：`marp slides.md --pdf --theme course-theme.css`

### Q5: 如何获取最新版本的模板？
A: 从GitHub仓库获取最新版本：`git pull origin main` 或访问 https://github.com/hnnd/ngs-for-master

这个使用指南提供了完整的模板使用方法，确保教学材料的标准化和高质量制作。项目托管在 https://github.com/hnnd/ngs-for-master，欢迎贡献和反馈。