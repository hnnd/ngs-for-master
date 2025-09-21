#!/bin/bash
# 第8次课环境设置脚本
# 作者：王运生
# 日期：2025-01-21
# 用法：bash setup.sh

set -e  # 遇到错误立即退出
set -u  # 使用未定义变量时报错

echo "=== 第8次课：多组学数据整合与机器学习 - 环境设置 ==="

# 创建工作目录结构
echo "创建工作目录结构..."
mkdir -p ~/ngs-analysis/lesson-08/{data,scripts,results,plots,models}
cd ~/ngs-analysis/lesson-08

# 检查R环境
echo "检查R环境..."
if command -v R &> /dev/null; then
    echo "✓ R已安装: $(R --version | head -1)"
else
    echo "✗ R未安装，请先安装R"
    exit 1
fi

# 检查Python环境
echo "检查Python环境..."
if command -v python &> /dev/null; then
    echo "✓ Python已安装: $(python --version)"
else
    echo "✗ Python未安装，请先安装Python"
    exit 1
fi

# 检查pip
if command -v pip &> /dev/null; then
    echo "✓ pip已安装: $(pip --version)"
else
    echo "✗ pip未安装，请先安装pip"
    exit 1
fi

# 安装Python包
echo "安装Python包..."
pip install --quiet pandas numpy scikit-learn matplotlib seaborn
pip install --quiet tensorflow keras plotly umap-learn
pip install --quiet jupyter jupyterlab
pip install --quiet joblib

echo "✓ Python包安装完成"

# 创建R包安装脚本
echo "创建R包安装脚本..."
cat > install_r_packages.R << 'EOF'
# 安装必需的R包
cat("安装R包...\n")

# 设置CRAN镜像
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# 基础包
packages <- c("ggplot2", "pheatmap", "VennDiagram", "caret", "randomForest")
for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
        install.packages(pkg, dependencies = TRUE)
        cat("✓", pkg, "安装完成\n")
    } else {
        cat("✓", pkg, "已安装\n")
    }
}

# Bioconductor包
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

bioc_packages <- c("limma", "edgeR", "DESeq2", "MultiAssayExperiment")
for (pkg in bioc_packages) {
    if (!require(pkg, character.only = TRUE)) {
        BiocManager::install(pkg)
        cat("✓", pkg, "安装完成\n")
    } else {
        cat("✓", pkg, "已安装\n")
    }
}

# mixOmics包
if (!require("mixOmics", character.only = TRUE)) {
    install.packages("mixOmics")
    cat("✓ mixOmics 安装完成\n")
} else {
    cat("✓ mixOmics 已安装\n")
}

cat("所有R包安装完成！\n")
EOF

# 运行R包安装
echo "安装R包..."
Rscript install_r_packages.R

# 生成示例数据
echo "生成示例数据..."
python << 'EOF'
import pandas as pd
import numpy as np

# 设置随机种子
np.random.seed(42)

# 生成示例数据
n_samples = 50
n_genomics = 100
n_transcriptomics = 1000
n_proteomics = 500

# 基因组数据（SNP，0/1/2编码）
genomics_data = np.random.choice([0, 1, 2], size=(n_samples, n_genomics), p=[0.6, 0.3, 0.1])
genomics_df = pd.DataFrame(genomics_data, 
                          index=[f'Sample_{i+1}' for i in range(n_samples)],
                          columns=[f'SNP_{i+1}' for i in range(n_genomics)])

# 转录组数据（表达量，对数正态分布）
transcriptomics_data = np.random.lognormal(mean=5, sigma=2, size=(n_samples, n_transcriptomics))
transcriptomics_df = pd.DataFrame(transcriptomics_data,
                                 index=[f'Sample_{i+1}' for i in range(n_samples)],
                                 columns=[f'Gene_{i+1}' for i in range(n_transcriptomics)])

# 蛋白质组数据（丰度，正态分布）
proteomics_data = np.random.normal(loc=10, scale=3, size=(n_samples, n_proteomics))
proteomics_df = pd.DataFrame(proteomics_data,
                            index=[f'Sample_{i+1}' for i in range(n_samples)],
                            columns=[f'Protein_{i+1}' for i in range(n_proteomics)])

# 添加一些缺失值
transcriptomics_df.iloc[np.random.choice(n_samples, 5), np.random.choice(n_transcriptomics, 5)] = np.nan
proteomics_df.iloc[np.random.choice(n_samples, 3), np.random.choice(n_proteomics, 3)] = np.nan

# 样本信息
groups = ['group1'] * 25 + ['group2'] * 25
ages = np.random.randint(20, 80, n_samples)
genders = np.random.choice(['Male', 'Female'], n_samples)

metadata_df = pd.DataFrame({
    'group': groups,
    'age': ages,
    'gender': genders,
    'batch': np.random.choice(['batch1', 'batch2'], n_samples),
    'treatment': np.random.choice(['treated', 'control'], n_samples)
}, index=[f'Sample_{i+1}' for i in range(n_samples)])

# 保存数据
genomics_df.to_csv('data/genomics_data.csv')
transcriptomics_df.to_csv('data/transcriptomics_data.csv')
proteomics_df.to_csv('data/proteomics_data.csv')
metadata_df.to_csv('data/metadata.csv')

print("✓ 示例数据生成完成")
print(f"  - 基因组数据: {genomics_df.shape}")
print(f"  - 转录组数据: {transcriptomics_df.shape}")
print(f"  - 蛋白质组数据: {proteomics_df.shape}")
print(f"  - 样本信息: {metadata_df.shape}")
EOF

# 创建Jupyter配置
echo "配置Jupyter环境..."
jupyter --generate-config 2>/dev/null || true

# 创建快速启动脚本
cat > start_analysis.sh << 'EOF'
#!/bin/bash
# 快速启动分析环境
echo "启动多组学数据分析环境..."
cd ~/ngs-analysis/lesson-08

echo "可用的分析工具："
echo "1. 启动R: R"
echo "2. 启动Python: python"
echo "3. 启动Jupyter: jupyter lab"
echo "4. 查看数据: ls -la data/"

echo "开始分析！"
EOF

chmod +x start_analysis.sh

# 创建环境检查脚本
cat > check_environment.py << 'EOF'
#!/usr/bin/env python
"""环境检查脚本"""

import sys
import importlib

def check_package(package_name):
    try:
        importlib.import_module(package_name)
        return True
    except ImportError:
        return False

print("=== Python环境检查 ===")
print(f"Python版本: {sys.version}")

packages = [
    'pandas', 'numpy', 'sklearn', 'matplotlib', 'seaborn',
    'tensorflow', 'keras', 'plotly', 'umap', 'joblib'
]

for pkg in packages:
    status = "✓" if check_package(pkg) else "✗"
    print(f"{status} {pkg}")

print("\n=== 数据文件检查 ===")
import os
data_files = ['genomics_data.csv', 'transcriptomics_data.csv', 'proteomics_data.csv', 'metadata.csv']
for file in data_files:
    path = f'data/{file}'
    status = "✓" if os.path.exists(path) else "✗"
    if os.path.exists(path):
        size = os.path.getsize(path) / 1024  # KB
        print(f"{status} {file} ({size:.1f} KB)")
    else:
        print(f"{status} {file} (缺失)")

print("\n环境检查完成！")
EOF

chmod +x check_environment.py

# 运行环境检查
echo "运行环境检查..."
python check_environment.py

# 清理临时文件
rm -f install_r_packages.R

echo ""
echo "=== 环境设置完成 ==="
echo "工作目录: ~/ngs-analysis/lesson-08"
echo "启动分析: bash start_analysis.sh"
echo "环境检查: python check_environment.py"
echo ""
echo "准备开始第8次课的学习！"