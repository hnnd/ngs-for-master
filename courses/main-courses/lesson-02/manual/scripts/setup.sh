#!/bin/bash

# 第2次课环境设置脚本
# 课程：高通量测序数据分析 - 测序数据质量控制与预处理
# 作者：王运生
# 日期：2025-01-20
# 用法：bash setup.sh

set -e  # 遇到错误立即退出
set -u  # 使用未定义变量时报错

echo "=========================================="
echo "第2次课：测序数据质量控制与预处理"
echo "环境设置脚本"
echo "=========================================="

# 检查操作系统
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    echo "检测到Linux系统"
    PACKAGE_MANAGER="apt-get"
elif [[ "$OSTYPE" == "darwin"* ]]; then
    echo "检测到macOS系统"
    PACKAGE_MANAGER="brew"
else
    echo "警告：未识别的操作系统，可能需要手动安装软件"
fi

# 创建工作目录结构
echo "创建工作目录结构..."
mkdir -p ~/ngs-analysis/lesson-02/{data,scripts,results,logs}
mkdir -p ~/ngs-analysis/lesson-02/results/{fastqc_raw,fastqc_clean,multiqc,trimmed}

echo "工作目录创建完成：~/ngs-analysis/lesson-02/"

# 检查并安装必要软件
echo "检查软件环境..."

# 检查Java环境
echo "检查Java环境..."
if command -v java &> /dev/null; then
    JAVA_VERSION=$(java -version 2>&1 | head -n 1 | cut -d'"' -f2)
    echo "Java已安装，版本：$JAVA_VERSION"
else
    echo "Java未安装，正在安装..."
    if [[ "$PACKAGE_MANAGER" == "apt-get" ]]; then
        sudo apt-get update
        sudo apt-get install -y openjdk-8-jdk
    elif [[ "$PACKAGE_MANAGER" == "brew" ]]; then
        brew install openjdk@8
    fi
fi

# 检查Python环境
echo "检查Python环境..."
if command -v python3 &> /dev/null; then
    PYTHON_VERSION=$(python3 --version)
    echo "Python已安装：$PYTHON_VERSION"
else
    echo "错误：Python3未安装，请先安装Python3"
    exit 1
fi

# 检查pip
if command -v pip3 &> /dev/null; then
    echo "pip3已安装"
else
    echo "安装pip3..."
    if [[ "$PACKAGE_MANAGER" == "apt-get" ]]; then
        sudo apt-get install -y python3-pip
    elif [[ "$PACKAGE_MANAGER" == "brew" ]]; then
        echo "pip应该随Python一起安装"
    fi
fi

# 检查conda（可选）
if command -v conda &> /dev/null; then
    echo "Conda已安装，可以使用conda安装生物信息学软件"
    CONDA_AVAILABLE=true
else
    echo "Conda未安装，将使用其他方式安装软件"
    CONDA_AVAILABLE=false
fi

# 安装FastQC
echo "检查FastQC..."
if command -v fastqc &> /dev/null; then
    FASTQC_VERSION=$(fastqc --version)
    echo "FastQC已安装：$FASTQC_VERSION"
else
    echo "安装FastQC..."
    if [[ "$CONDA_AVAILABLE" == true ]]; then
        conda install -c bioconda fastqc -y
    elif [[ "$PACKAGE_MANAGER" == "apt-get" ]]; then
        sudo apt-get install -y fastqc
    elif [[ "$PACKAGE_MANAGER" == "brew" ]]; then
        brew install fastqc
    else
        echo "请手动安装FastQC：https://www.bioinformatics.babraham.ac.uk/projects/fastqc/"
    fi
fi

# 安装MultiQC
echo "检查MultiQC..."
if command -v multiqc &> /dev/null; then
    MULTIQC_VERSION=$(multiqc --version)
    echo "MultiQC已安装：$MULTIQC_VERSION"
else
    echo "安装MultiQC..."
    pip3 install --user multiqc
    # 添加到PATH
    echo 'export PATH=$PATH:~/.local/bin' >> ~/.bashrc
    export PATH=$PATH:~/.local/bin
fi

# 安装Trimmomatic
echo "检查Trimmomatic..."
if [[ "$CONDA_AVAILABLE" == true ]]; then
    if conda list trimmomatic &> /dev/null; then
        echo "Trimmomatic已通过conda安装"
    else
        echo "通过conda安装Trimmomatic..."
        conda install -c bioconda trimmomatic -y
    fi
elif command -v trimmomatic &> /dev/null; then
    echo "Trimmomatic已安装"
else
    echo "下载Trimmomatic..."
    cd ~/ngs-analysis/lesson-02/
    if [ ! -f "Trimmomatic-0.39.jar" ]; then
        wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
        unzip Trimmomatic-0.39.zip
        mv Trimmomatic-0.39/trimmomatic-0.39.jar ./
        rm -rf Trimmomatic-0.39*
    fi
    
    # 创建trimmomatic命令别名
    echo "创建Trimmomatic启动脚本..."
    cat > ~/ngs-analysis/lesson-02/scripts/trimmomatic << 'EOF'
#!/bin/bash
java -jar ~/ngs-analysis/lesson-02/trimmomatic-0.39.jar "$@"
EOF
    chmod +x ~/ngs-analysis/lesson-02/scripts/trimmomatic
    echo 'export PATH=$PATH:~/ngs-analysis/lesson-02/scripts' >> ~/.bashrc
    export PATH=$PATH:~/ngs-analysis/lesson-02/scripts
fi

# 下载接头序列文件
echo "准备接头序列文件..."
cd ~/ngs-analysis/lesson-02/data/

if [ ! -f "TruSeq3-PE.fa" ]; then
    cat > TruSeq3-PE.fa << 'EOF'
>PrefixPE/1
TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PrefixPE/2
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
>PE1
TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PE1_rc
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
>PE2
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
>PE2_rc
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
EOF
    echo "TruSeq3-PE.fa接头文件已创建"
fi

# 创建测试数据（如果没有真实数据）
echo "准备测试数据..."
if [ ! -f "sample_R1.fastq.gz" ] && [ ! -f "sample_R2.fastq.gz" ]; then
    echo "创建模拟测试数据..."
    python3 << 'EOF'
import gzip
import random

def generate_fastq_read(read_id, length=150):
    """生成一条FASTQ格式的序列"""
    bases = ['A', 'T', 'G', 'C']
    sequence = ''.join(random.choices(bases, k=length))
    
    # 模拟质量分数，前面质量高，后面质量逐渐降低
    qualities = []
    for i in range(length):
        if i < length * 0.7:
            q = random.randint(30, 40)  # 高质量
        else:
            q = random.randint(15, 35)  # 质量逐渐降低
        qualities.append(chr(q + 33))  # Phred+33编码
    
    quality_string = ''.join(qualities)
    
    return f"@{read_id}\n{sequence}\n+\n{quality_string}\n"

# 生成R1文件
print("生成sample_R1.fastq.gz...")
with gzip.open('sample_R1.fastq.gz', 'wt') as f:
    for i in range(10000):  # 生成10000条序列
        read = generate_fastq_read(f"read_{i+1}/1")
        f.write(read)

# 生成R2文件
print("生成sample_R2.fastq.gz...")
with gzip.open('sample_R2.fastq.gz', 'wt') as f:
    for i in range(10000):  # 生成10000条序列
        read = generate_fastq_read(f"read_{i+1}/2")
        f.write(read)

print("测试数据生成完成")
EOF
fi

# 验证安装
echo "验证软件安装..."
echo "FastQC版本："
fastqc --version || echo "FastQC安装失败"

echo "MultiQC版本："
multiqc --version || echo "MultiQC安装失败"

echo "Trimmomatic版本："
if [[ "$CONDA_AVAILABLE" == true ]]; then
    trimmomatic -version || echo "Trimmomatic安装失败"
else
    java -jar ~/ngs-analysis/lesson-02/trimmomatic-0.39.jar -version || echo "Trimmomatic安装失败"
fi

echo "Java版本："
java -version

# 检查数据文件
echo "检查数据文件..."
ls -lh ~/ngs-analysis/lesson-02/data/

echo "=========================================="
echo "环境设置完成！"
echo "工作目录：~/ngs-analysis/lesson-02/"
echo "请运行以下命令激活环境变量："
echo "source ~/.bashrc"
echo "或者重新打开终端"
echo "=========================================="

# 创建环境检查脚本
cat > ~/ngs-analysis/lesson-02/scripts/check_environment.sh << 'EOF'
#!/bin/bash
echo "检查实验环境..."
echo "FastQC: $(which fastqc)"
echo "MultiQC: $(which multiqc)"
echo "Trimmomatic: $(which trimmomatic)"
echo "Java: $(which java)"
echo "Python3: $(which python3)"
echo "工作目录："
ls -la ~/ngs-analysis/lesson-02/
EOF

chmod +x ~/ngs-analysis/lesson-02/scripts/check_environment.sh

echo "环境检查脚本已创建：~/ngs-analysis/lesson-02/scripts/check_environment.sh"
echo "可以随时运行该脚本检查环境状态"