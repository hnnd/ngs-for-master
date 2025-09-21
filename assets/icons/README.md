# 图标资源说明

## 图标分类

### 功能图标
- `file-icon.svg` - 文件类型图标
- `folder-icon.svg` - 文件夹图标
- `download-icon.svg` - 下载图标
- `upload-icon.svg` - 上传图标
- `settings-icon.svg` - 设置图标

### 工具图标
- `terminal-icon.svg` - 终端/命令行图标
- `code-icon.svg` - 代码编辑图标
- `database-icon.svg` - 数据库图标
- `chart-icon.svg` - 图表分析图标
- `microscope-icon.svg` - 科研工具图标

### 状态图标
- `success-icon.svg` - 成功状态
- `warning-icon.svg` - 警告状态
- `error-icon.svg` - 错误状态
- `info-icon.svg` - 信息提示

## 设计规范

### 尺寸标准
- **小图标**：16x16px, 24x24px
- **中图标**：32x32px, 48x48px
- **大图标**：64x64px, 128x128px

### 风格要求
- 简洁明了，易于识别
- 统一的线条粗细（2px）
- 圆角半径统一（2px）
- 颜色使用课程主题色

### 颜色方案
- **主色**：#2E86AB（蓝色）
- **辅助色**：#A23B72（紫色）
- **中性色**：#F18F01（橙色）
- **灰色**：#C73E1D（深红）

## 使用方法

### 在Markdown中使用
```markdown
![图标描述](../../../assets/icons/icon-name.svg)
```

### 在HTML中使用
```html
<img src="../../../assets/icons/icon-name.svg" alt="图标描述" width="24" height="24">
```

### 在CSS中使用
```css
.icon {
    background-image: url('../../../assets/icons/icon-name.svg');
    width: 24px;
    height: 24px;
}
```

## 创建新图标

1. 使用SVG格式创建
2. 遵循设计规范
3. 优化文件大小
4. 添加适当的注释
5. 更新本README文档