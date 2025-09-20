#!/bin/bash

# 高通量测序数据分析课程教学材料仓库设置脚本
# 作者：王运生教授
# 邮箱：wangys@hunau.edu.cn

set -e  # 遇到错误时退出

echo "🚀 开始设置高通量测序数据分析课程教学材料仓库..."

# 检查是否在正确的目录中
if [ ! -f "README.md" ] || [ ! -d "templates" ]; then
    echo "❌ 错误：请在 ngs-course-materials 目录中运行此脚本"
    exit 1
fi

# 检查Git是否已安装
if ! command -v git &> /dev/null; then
    echo "❌ 错误：Git 未安装，请先安装 Git"
    exit 1
fi

# 初始化Git仓库（如果尚未初始化）
if [ ! -d ".git" ]; then
    echo "📦 初始化 Git 仓库..."
    git init
    echo "✅ Git 仓库初始化完成"
else
    echo "📦 Git 仓库已存在"
fi

# 设置远程仓库
REMOTE_URL="https://github.com/hnnd/ngs-for-master.git"
if git remote get-url origin &> /dev/null; then
    CURRENT_URL=$(git remote get-url origin)
    if [ "$CURRENT_URL" != "$REMOTE_URL" ]; then
        echo "🔄 更新远程仓库地址..."
        git remote set-url origin "$REMOTE_URL"
        echo "✅ 远程仓库地址已更新为: $REMOTE_URL"
    else
        echo "📡 远程仓库地址已正确设置"
    fi
else
    echo "📡 添加远程仓库..."
    git remote add origin "$REMOTE_URL"
    echo "✅ 远程仓库已添加: $REMOTE_URL"
fi

# 检查Git LFS是否已安装
if command -v git-lfs &> /dev/null; then
    echo "📦 初始化 Git LFS..."
    git lfs install
    echo "✅ Git LFS 已初始化"
else
    echo "⚠️  警告：Git LFS 未安装，大文件管理可能受限"
    echo "   请安装 Git LFS: https://git-lfs.github.io/"
fi

# 设置Git配置（如果尚未设置）
if [ -z "$(git config user.name)" ]; then
    echo "👤 设置Git用户信息..."
    read -p "请输入您的姓名: " user_name
    read -p "请输入您的邮箱: " user_email
    git config user.name "$user_name"
    git config user.email "$user_email"
    echo "✅ Git用户信息已设置"
fi

# 创建必要的目录结构（如果不存在）
echo "📁 检查目录结构..."
REQUIRED_DIRS=(
    "assets/images/common"
    "assets/images/sequencing"
    "assets/images/analysis"
    "assets/images/results"
    "assets/logos"
    "assets/icons"
    "courses/prep-courses"
    "courses/main-courses"
    "scripts/generators"
    "scripts/validators"
    "scripts/utilities"
)

for dir in "${REQUIRED_DIRS[@]}"; do
    if [ ! -d "$dir" ]; then
        mkdir -p "$dir"
        echo "✅ 创建目录: $dir"
    fi
done

# 检查Node.js和npm（用于Marp）
if command -v npm &> /dev/null; then
    echo "📦 检查Marp CLI..."
    if ! npm list -g @marp-team/marp-cli &> /dev/null; then
        echo "🔧 安装Marp CLI..."
        npm install -g @marp-team/marp-cli
        echo "✅ Marp CLI 已安装"
    else
        echo "✅ Marp CLI 已安装"
    fi
else
    echo "⚠️  警告：Node.js/npm 未安装，无法使用Marp生成幻灯片"
    echo "   请安装 Node.js: https://nodejs.org/"
fi

# 检查Python环境
if command -v python3 &> /dev/null; then
    echo "🐍 Python 环境已就绪"
    
    # 创建requirements.txt（如果不存在）
    if [ ! -f "requirements.txt" ]; then
        echo "📝 创建 requirements.txt..."
        cat > requirements.txt << EOF
# 高通量测序数据分析课程教学材料依赖包

# 文档生成和处理
markdown>=3.4.0
pymdown-extensions>=9.0.0
mkdocs>=1.4.0
mkdocs-material>=8.0.0

# 数据处理和分析
pandas>=1.5.0
numpy>=1.21.0
matplotlib>=3.5.0
seaborn>=0.11.0

# 生物信息学相关
biopython>=1.79
pysam>=0.19.0

# 图像处理
pillow>=9.0.0
cairosvg>=2.5.0

# 实用工具
requests>=2.28.0
pyyaml>=6.0
click>=8.0.0

# 开发和测试工具
pytest>=7.0.0
black>=22.0.0
flake8>=4.0.0
EOF
        echo "✅ requirements.txt 已创建"
    fi
else
    echo "⚠️  警告：Python 未安装，脚本功能可能受限"
    echo "   请安装 Python 3.7+: https://www.python.org/"
fi

# 添加所有文件到Git
echo "📝 添加文件到Git..."
git add .

# 创建初始提交（如果是新仓库）
if [ -z "$(git log --oneline 2>/dev/null)" ]; then
    echo "💾 创建初始提交..."
    git commit -m "feat: 初始化高通量测序数据分析课程教学材料

- 创建完整的目录结构和模板系统
- 添加Marp幻灯片模板和CSS主题
- 建立实践手册标准模板
- 完善文档和使用指南
- 配置Git和Git LFS设置"
    echo "✅ 初始提交已创建"
fi

# 显示仓库状态
echo ""
echo "📊 仓库状态:"
echo "   远程仓库: $(git remote get-url origin)"
echo "   当前分支: $(git branch --show-current)"
echo "   提交数量: $(git rev-list --count HEAD 2>/dev/null || echo '0')"

echo ""
echo "🎉 仓库设置完成！"
echo ""
echo "📚 下一步操作:"
echo "   1. 推送到GitHub: git push -u origin main"
echo "   2. 开始创建课程内容"
echo "   3. 查看使用指南: docs/template-usage.md"
echo ""
echo "🔗 项目地址: https://github.com/hnnd/ngs-for-master"
echo "📧 联系方式: wangys@hunau.edu.cn"