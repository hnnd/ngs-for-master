#!/usr/bin/env python3
"""
Python环境设置和验证脚本
课程：高通量测序数据分析 - Python编程基础
作者：王运生教授
日期：2025年
"""

import sys
import subprocess
import importlib
import os

def check_python_version():
    """检查Python版本"""
    print("=== Python版本检查 ===")
    version = sys.version_info
    print(f"Python版本: {version.major}.{version.minor}.{version.micro}")
    
    if version.major >= 3 and version.minor >= 8:
        print("✓ Python版本符合要求 (3.8+)")
        return True
    else:
        print("✗ Python版本过低，需要3.8或更高版本")
        return False

def check_packages():
    """检查必需的Python包"""
    print("\n=== 必需包检查 ===")
    required_packages = [
        'numpy',
        'pandas', 
        'matplotlib',
        'seaborn',
        'jupyter'
    ]
    
    optional_packages = [
        'biopython'
    ]
    
    missing_packages = []
    
    # 检查必需包
    for package in required_packages:
        try:
            importlib.import_module(package)
            print(f"✓ {package} 已安装")
        except ImportError:
            print(f"✗ {package} 未安装")
            missing_packages.append(package)
    
    # 检查可选包
    for package in optional_packages:
        try:
            if package == 'biopython':
                importlib.import_module('Bio')
                print(f"✓ {package} 已安装")
            else:
                importlib.import_module(package)
                print(f"✓ {package} 已安装")
        except ImportError:
            print(f"⚠ {package} 未安装 (可选)")
    
    return missing_packages

def install_packages(packages):
    """安装缺失的包"""
    if not packages:
        return True
    
    print(f"\n=== 安装缺失的包 ===")
    print(f"需要安装: {', '.join(packages)}")
    
    try:
        # 升级pip
        print("升级pip...")
        subprocess.check_call([sys.executable, '-m', 'pip', 'install', '--upgrade', 'pip'])
        
        # 安装包
        for package in packages:
            print(f"安装 {package}...")
            if package == 'biopython':
                subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'biopython'])
            else:
                subprocess.check_call([sys.executable, '-m', 'pip', 'install', package])
        
        print("✓ 所有包安装完成")
        return True
    except subprocess.CalledProcessError as e:
        print(f"✗ 包安装失败: {e}")
        return False

def create_sample_data():
    """创建示例数据文件"""
    print("\n=== 创建示例数据文件 ===")
    
    # 创建data目录
    data_dir = "data"
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
        print(f"✓ 创建目录: {data_dir}")
    
    # 创建示例FASTA文件
    fasta_content = """>gene1|BRCA1|breast_cancer_gene
ATGAAATTTCCCGGGAAATAGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
CCCGGGAAATTTAAATTTCCCGGGAAATAGATCGATCGATCGATCGATCGAT
>gene2|TP53|tumor_suppressor
ATGCCCGGGAAATTTAAATTTCCCGGGAAATAGATCGATCGATCGATCGATC
GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
AAATTTCCCGGGAAATAGATCGATCGATCGATCGATCGATCGATCGATCGAT
>gene3|EGFR|growth_factor_receptor
ATGAAATTTCCCGGGAAATAGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
TTTCCCGGGAAATAGATCGATCGATCGATCGATCGATCGATCGATCGATCGA
>gene4|MYC|oncogene
ATGCCCGGGAAATTTAAATTTCCCGGGAAATAGATCGATCGATCGATCGATC
GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
CCCGGGAAATTTAAATTTCCCGGGAAATAGATCGATCGATCGATCGATCGAT
>gene5|KRAS|oncogene
ATGAAATTTCCCGGGAAATAGATCGATCGATCGATCGATCGATCGATCGATCG
CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT
AAATTTCCCGGGAAATAGATCGATCGATCGATCGATCGATCGATCGATCGAT"""
    
    fasta_file = os.path.join(data_dir, "sample_sequences.fasta")
    with open(fasta_file, "w") as f:
        f.write(fasta_content)
    print(f"✓ 创建文件: {fasta_file}")
    
    # 创建基因表达CSV文件
    import csv
    
    expression_data = [
        ["gene_name", "sample1", "sample2", "sample3", "sample4", "sample5"],
        ["BRCA1", "5.2", "4.8", "5.5", "4.9", "5.1"],
        ["BRCA2", "3.8", "4.1", "3.5", "3.9", "4.0"],
        ["TP53", "7.1", "6.8", "7.3", "7.0", "6.9"],
        ["EGFR", "4.5", "4.2", "4.8", "4.6", "4.4"],
        ["MYC", "6.3", "6.0", "6.5", "6.1", "6.2"],
        ["KRAS", "4.8", "5.1", "4.6", "4.9", "5.0"],
        ["PIK3CA", "5.5", "5.2", "5.8", "5.4", "5.6"]
    ]
    
    csv_file = os.path.join(data_dir, "gene_expression.csv")
    with open(csv_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(expression_data)
    print(f"✓ 创建文件: {csv_file}")
    
    # 创建质量分数文件
    quality_content = """# 测序质量分数数据
# 每行代表一个读段的质量分数
35 38 40 42 38 35 30 28 25 22 20 18 15 12 10
40 42 45 43 40 38 35 32 30 28 25 23 20 18 15
38 40 42 45 42 40 38 35 33 30 28 25 22 20 17
42 45 47 45 43 40 38 35 32 30 27 25 22 19 16
36 38 40 42 40 38 35 32 30 27 25 22 20 17 15
39 41 43 45 43 41 38 36 33 31 28 26 23 21 18
37 39 41 43 41 39 36 34 31 29 26 24 21 19 16
41 43 45 47 45 43 40 38 35 33 30 28 25 23 20"""
    
    quality_file = os.path.join(data_dir, "quality_scores.txt")
    with open(quality_file, "w") as f:
        f.write(quality_content)
    print(f"✓ 创建文件: {quality_file}")

def test_installation():
    """测试安装是否成功"""
    print("\n=== 安装测试 ===")
    
    try:
        # 测试NumPy
        import numpy as np
        arr = np.array([1, 2, 3, 4, 5])
        print(f"✓ NumPy测试通过: {arr.mean()}")
        
        # 测试Pandas
        import pandas as pd
        df = pd.DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]})
        print(f"✓ Pandas测试通过: DataFrame形状 {df.shape}")
        
        # 测试Matplotlib
        import matplotlib.pyplot as plt
        plt.figure()
        plt.close()
        print("✓ Matplotlib测试通过")
        
        # 测试BioPython（如果安装了）
        try:
            from Bio.Seq import Seq
            seq = Seq("ATCG")
            print(f"✓ BioPython测试通过: {seq}")
        except ImportError:
            print("⚠ BioPython未安装，跳过测试")
        
        return True
    except Exception as e:
        print(f"✗ 测试失败: {e}")
        return False

def main():
    """主函数"""
    print("Python环境设置和验证脚本")
    print("=" * 50)
    
    # 检查Python版本
    if not check_python_version():
        print("\n请升级Python版本后重新运行此脚本")
        return False
    
    # 检查包
    missing_packages = check_packages()
    
    # 安装缺失的包
    if missing_packages:
        install_choice = input(f"\n是否安装缺失的包? (y/n): ").lower()
        if install_choice == 'y':
            if not install_packages(missing_packages):
                print("包安装失败，请手动安装")
                return False
        else:
            print("跳过包安装")
    
    # 创建示例数据
    create_sample_data()
    
    # 测试安装
    if test_installation():
        print("\n" + "=" * 50)
        print("✓ 环境设置完成！可以开始Python编程学习")
        print("=" * 50)
        return True
    else:
        print("\n" + "=" * 50)
        print("✗ 环境设置存在问题，请检查错误信息")
        print("=" * 50)
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)