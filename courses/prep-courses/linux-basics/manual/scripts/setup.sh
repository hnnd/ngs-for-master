#!/bin/bash
# Linux基础课程环境设置脚本
# 作者：王运生教授
# 日期：2025-12-20
# 用法：bash setup.sh

set -e  # 遇到错误立即退出
set -u  # 使用未定义变量时报错

echo "=== Linux基础课程环境设置 ==="
echo "开始时间: $(date)"
echo

# 检查是否为root用户
if [ "$EUID" -eq 0 ]; then
    echo "警告：不建议以root用户运行此脚本"
    read -p "是否继续？(y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

# 更新包列表
echo "1. 更新软件包列表..."
sudo apt update

# 安装必要的软件包
echo "2. 安装必要的软件包..."
PACKAGES=(
    "tree"
    "htop"
    "nano"
    "vim"
    "curl"
    "wget"
    "git"
    "bc"
    "zip"
    "unzip"
    "net-tools"
)

for package in "${PACKAGES[@]}"; do
    if ! dpkg -l | grep -q "^ii  $package "; then
        echo "   安装 $package..."
        sudo apt install -y "$package"
    else
        echo "   $package 已安装"
    fi
done

# 创建工作目录结构
echo "3. 创建工作目录结构..."
WORK_DIR="$HOME/linux-basics-lab"
mkdir -p "$WORK_DIR"/{data,scripts,results,logs,backup}

# 创建示例数据文件
echo "4. 创建示例数据文件..."
cd "$WORK_DIR/data"

# 创建FASTA文件
cat > sample_sequences.fasta << 'EOF'
>sequence1|length=64|GC=50%
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>sequence2|length=64|GC=50%
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
>sequence3|length=64|GC=25%
TTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAA
>sequence4|length=96|GC=66%
CCCGGGCCCGGGCCCGGGCCCGGGCCCGGGCCCGGGCCCGGGCCCGGGCCCGGGCCCGGGCCCGGGCCCGGGCCCGGGCCCGGGCCCGGGCCCGGG
>sequence5|length=48|GC=33%
AAATTTAAATTTAAATTTAAATTTAAATTTAAATTTAAATTTAAATTT
EOF

# 创建基因表达数据
cat > gene_expression.csv << 'EOF'
gene_id,sample1,sample2,sample3,sample4,description
GENE001,12.5,15.2,8.7,11.3,Transcription factor
GENE002,45.1,42.8,48.3,44.7,Metabolic enzyme
GENE003,2.1,1.8,2.5,2.0,Structural protein
GENE004,78.9,82.1,75.6,79.4,Heat shock protein
GENE005,0.5,0.3,0.8,0.6,Hypothetical protein
GENE006,23.4,25.1,21.8,24.6,Ribosomal protein
GENE007,67.2,69.8,65.4,68.1,DNA polymerase
GENE008,5.7,6.2,5.1,5.9,Cell wall protein
GENE009,91.3,89.7,93.2,90.8,Chaperone protein
GENE010,14.6,13.9,15.3,14.2,Transport protein
EOF

# 创建质量分数文件
cat > quality_scores.txt << 'EOF'
@Read1
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@Read2
HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
+
HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
@Read3
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
+
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
@Read4
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
@Read5
EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
+
EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
EOF

# 创建大文件示例（用于性能测试）
echo "5. 创建大文件示例..."
seq 1 10000 | awk '{print "line_" $1 "_data_" rand()}' > large_dataset.txt

# 设置文件权限
echo "6. 设置文件权限..."
chmod 644 *.{fasta,csv,txt}

# 创建符号链接示例
ln -sf sample_sequences.fasta sequences_link.fasta

# 验证安装
echo "7. 验证环境设置..."
cd "$WORK_DIR"

echo "   检查目录结构..."
tree -L 2

echo "   检查数据文件..."
ls -lh data/

echo "   检查软件版本..."
echo "   - Bash: $(bash --version | head -n1)"
echo "   - Tree: $(tree --version | head -n1)"
echo "   - Git: $(git --version)"

# 创建环境信息文件
echo "8. 创建环境信息文件..."
cat > logs/environment_info.txt << EOF
Linux基础课程环境信息
生成时间: $(date)
用户: $(whoami)
主目录: $HOME
工作目录: $WORK_DIR
操作系统: $(cat /etc/os-release | grep PRETTY_NAME | cut -d'"' -f2)
内核版本: $(uname -r)
Shell: $SHELL
已安装软件包:
$(dpkg -l | grep -E "tree|htop|nano|vim|curl|wget|git|bc")
EOF

echo
echo "=== 环境设置完成 ==="
echo "工作目录: $WORK_DIR"
echo "请运行以下命令开始实验:"
echo "cd $WORK_DIR"
echo "ls -la"
echo
echo "如有问题，请查看日志文件: $WORK_DIR/logs/environment_info.txt"
echo "完成时间: $(date)"