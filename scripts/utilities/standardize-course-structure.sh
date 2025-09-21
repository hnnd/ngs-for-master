#!/bin/bash

# 标准化课程目录结构脚本
# 为所有课程创建统一的目录结构

echo "开始标准化课程目录结构..."

# 定义标准目录结构
STANDARD_DIRS=(
    "slides"
    "manual"
    "manual/scripts"
    "manual/data"
    "manual/data/sample_data"
    "manual/data/expected_results"
    "images"
    "images/diagrams"
    "images/screenshots"
    "images/results"
    "images/workflows"
)

# 处理主课程
echo "处理主课程目录..."
for lesson_dir in ngs-course-materials/courses/main-courses/lesson-*; do
    if [ -d "$lesson_dir" ]; then
        echo "处理 $(basename "$lesson_dir")..."
        
        # 创建标准目录
        for dir in "${STANDARD_DIRS[@]}"; do
            mkdir -p "$lesson_dir/$dir"
        done
        
        # 创建数据目录的README文件
        if [ ! -f "$lesson_dir/manual/data/README.md" ]; then
            cat > "$lesson_dir/manual/data/README.md" << EOF
# $(basename "$lesson_dir") 数据文件说明

## 目录结构

- \`sample_data/\` - 示例数据文件
- \`expected_results/\` - 预期结果文件

## 数据来源

请在此处说明数据的来源、格式和用途。

## 使用说明

1. 示例数据用于课程实践操作
2. 预期结果用于验证分析正确性
3. 数据文件较大时建议使用下载脚本获取

## 注意事项

- 请勿修改原始数据文件
- 分析结果保存在单独的输出目录
- 定期清理临时文件
EOF
        fi
        
        # 为images子目录创建.gitkeep文件
        for img_dir in diagrams screenshots results workflows; do
            if [ ! -f "$lesson_dir/images/$img_dir/.gitkeep" ]; then
                echo "# $(basename "$lesson_dir") $img_dir 图片目录" > "$lesson_dir/images/$img_dir/.gitkeep"
            fi
        done
    fi
done

# 处理预备课程
echo "处理预备课程目录..."
for prep_dir in ngs-course-materials/courses/prep-courses/*; do
    if [ -d "$prep_dir" ] && [ "$(basename "$prep_dir")" != "molecular-biology-basics" ]; then
        echo "处理 $(basename "$prep_dir")..."
        
        # 创建标准目录
        for dir in "${STANDARD_DIRS[@]}"; do
            mkdir -p "$prep_dir/$dir"
        done
        
        # 创建数据目录的README文件
        if [ ! -f "$prep_dir/manual/data/README.md" ]; then
            cat > "$prep_dir/manual/data/README.md" << EOF
# $(basename "$prep_dir") 数据文件说明

## 目录结构

- \`sample_data/\` - 示例数据文件
- \`expected_results/\` - 预期结果文件

## 数据来源

请在此处说明数据的来源、格式和用途。

## 使用说明

1. 示例数据用于课程实践操作
2. 预期结果用于验证分析正确性
3. 数据文件较大时建议使用下载脚本获取

## 注意事项

- 请勿修改原始数据文件
- 分析结果保存在单独的输出目录
- 定期清理临时文件
EOF
        fi
        
        # 为images子目录创建.gitkeep文件
        for img_dir in diagrams screenshots results workflows; do
            if [ ! -f "$prep_dir/images/$img_dir/.gitkeep" ]; then
                echo "# $(basename "$prep_dir") $img_dir 图片目录" > "$prep_dir/images/$img_dir/.gitkeep"
            fi
        done
    fi
done

echo "课程目录结构标准化完成！"

# 显示目录结构
echo ""
echo "标准化后的目录结构示例："
tree ngs-course-materials/courses/main-courses/lesson-01 2>/dev/null || find ngs-course-materials/courses/main-courses/lesson-01 -type d | sort