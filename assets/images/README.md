# 图片资源说明

## 目录结构

```
images/
├── common/              # 通用图表和图标
├── sequencing/          # 测序技术相关图片
├── analysis/            # 分析流程图片
└── results/             # 结果展示图片
```

## 图片规范

### 文件格式
- **SVG格式**：用于图表、流程图、原理图等矢量图形
- **PNG格式**：用于截图、界面图片等位图
- **JPG格式**：用于照片、复杂图像等

### 命名规范
- 使用小写字母和连字符
- 描述性命名，见名知意
- 包含图片类型信息

### 示例
- `sequencing-workflow.svg` - 测序工作流程图
- `quality-control-process.svg` - 质量控制流程
- `fastqc-report-screenshot.png` - FastQC报告截图
- `alignment-statistics-chart.svg` - 比对统计图表

## 使用说明

### 在幻灯片中使用
```markdown
![图片描述](../../../assets/images/category/image-name.svg)
```

### 在手册中使用
```markdown
![图片描述](../../../assets/images/category/image-name.png)
```

### 图片尺寸建议
- **幻灯片图片**：宽度不超过800px
- **手册图片**：宽度不超过600px
- **图标**：32x32px 或 64x64px

## 维护说明

- 定期检查图片链接有效性
- 更新过时的截图和界面图片
- 保持图片风格的一致性
- 优化图片文件大小