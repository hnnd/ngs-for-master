#!/usr/bin/env python3
"""
比对统计分析脚本
作者：王运生
日期：2025-01-20
用法：python alignment_stats.py
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import subprocess
import re
from collections import defaultdict

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

def run_command(cmd):
    """运行shell命令并返回输出"""
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        return result.stdout.strip()
    except Exception as e:
        print(f"命令执行错误: {cmd}")
        print(f"错误信息: {e}")
        return ""

def parse_flagstat(flagstat_file):
    """解析samtools flagstat输出"""
    stats = {}
    
    if not os.path.exists(flagstat_file):
        return stats
    
    with open(flagstat_file, 'r') as f:
        content = f.read()
    
    # 解析各种统计信息
    patterns = {
        'total_reads': r'(\d+) \+ \d+ in total',
        'mapped_reads': r'(\d+) \+ \d+ mapped \(',
        'paired_reads': r'(\d+) \+ \d+ paired in sequencing',
        'proper_pairs': r'(\d+) \+ \d+ properly paired',
        'singletons': r'(\d+) \+ \d+ singletons',
        'duplicates': r'(\d+) \+ \d+ duplicates'
    }
    
    for key, pattern in patterns.items():
        match = re.search(pattern, content)
        if match:
            stats[key] = int(match.group(1))
        else:
            stats[key] = 0
    
    # 计算比率
    if stats.get('total_reads', 0) > 0:
        stats['mapping_rate'] = stats.get('mapped_reads', 0) / stats['total_reads'] * 100
        stats['pairing_rate'] = stats.get('proper_pairs', 0) / stats['total_reads'] * 100
        stats['singleton_rate'] = stats.get('singletons', 0) / stats['total_reads'] * 100
        stats['duplicate_rate'] = stats.get('duplicates', 0) / stats['total_reads'] * 100
    
    return stats

def parse_mapq_distribution(mapq_file):
    """解析MAPQ分布文件"""
    mapq_dist = {}
    
    if not os.path.exists(mapq_file):
        return mapq_dist
    
    with open(mapq_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                count = int(parts[0])
                mapq = int(parts[1])
                mapq_dist[mapq] = count
    
    return mapq_dist

def parse_insert_sizes(insert_file):
    """解析插入片段大小文件"""
    insert_sizes = []
    
    if not os.path.exists(insert_file):
        return insert_sizes
    
    with open(insert_file, 'r') as f:
        for line in f:
            try:
                size = int(line.strip())
                if 0 < size < 2000:  # 过滤异常值
                    insert_sizes.append(size)
            except ValueError:
                continue
    
    return insert_sizes

def analyze_coverage(coverage_file):
    """分析覆盖度文件"""
    coverage_stats = {}
    
    if not os.path.exists(coverage_file):
        return coverage_stats
    
    coverages = []
    with open(coverage_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3:
                try:
                    cov = int(parts[2])
                    coverages.append(cov)
                except ValueError:
                    continue
    
    if coverages:
        coverage_stats = {
            'mean_coverage': np.mean(coverages),
            'median_coverage': np.median(coverages),
            'std_coverage': np.std(coverages),
            'min_coverage': np.min(coverages),
            'max_coverage': np.max(coverages),
            'covered_positions': len([c for c in coverages if c > 0]),
            'total_positions': len(coverages)
        }
        
        if coverage_stats['total_positions'] > 0:
            coverage_stats['coverage_rate'] = coverage_stats['covered_positions'] / coverage_stats['total_positions'] * 100
    
    return coverage_stats

def collect_alignment_statistics():
    """收集所有比对统计信息"""
    results_dir = Path("results")
    
    # 查找所有BAM文件
    bam_files = list(results_dir.glob("*_sorted.bam"))
    
    all_stats = {}
    
    for bam_file in bam_files:
        base_name = bam_file.stem.replace("_sorted", "")
        
        print(f"分析文件: {base_name}")
        
        # 解析flagstat文件
        flagstat_file = results_dir / f"{base_name}_stats.txt"
        stats = parse_flagstat(flagstat_file)
        
        # 解析MAPQ分布
        mapq_file = results_dir / f"{base_name}_mapq.txt"
        mapq_dist = parse_mapq_distribution(mapq_file)
        
        # 解析插入片段大小
        insert_file = results_dir / f"{base_name}_insert_sizes.txt"
        if not insert_file.exists():
            # 如果没有插入片段文件，尝试从BAM文件提取
            cmd = f"samtools view -f 2 {bam_file} | awk '$9 > 0 {{print $9}}' | head -10000 > {insert_file}"
            run_command(cmd)
        
        insert_sizes = parse_insert_sizes(insert_file)
        
        # 分析覆盖度
        coverage_file = results_dir / f"{base_name}_coverage.txt"
        coverage_stats = analyze_coverage(coverage_file)
        
        # 合并所有统计信息
        all_stats[base_name] = {
            'flagstat': stats,
            'mapq_distribution': mapq_dist,
            'insert_sizes': insert_sizes,
            'coverage': coverage_stats
        }
    
    return all_stats

def create_comparison_table(all_stats):
    """创建比对工具比较表"""
    comparison_data = []
    
    for sample_name, stats in all_stats.items():
        flagstat = stats['flagstat']
        coverage = stats['coverage']
        
        # 确定工具类型
        if 'bwa' in sample_name.lower():
            tool = 'BWA'
            if 'default' in sample_name:
                mode = 'Default'
            elif 'optimized' in sample_name:
                mode = 'Optimized'
            elif 'strict' in sample_name:
                mode = 'Strict'
            else:
                mode = 'Unknown'
        elif 'bowtie2' in sample_name.lower():
            tool = 'Bowtie2'
            if 'fast' in sample_name:
                mode = 'Fast'
            elif 'very_sensitive' in sample_name:
                mode = 'Very Sensitive'
            elif 'sensitive' in sample_name:
                mode = 'Sensitive'
            elif 'local' in sample_name:
                mode = 'Local'
            elif 'custom' in sample_name:
                mode = 'Custom'
            else:
                mode = 'Unknown'
        else:
            tool = 'Unknown'
            mode = 'Unknown'
        
        row = {
            'Sample': sample_name,
            'Tool': tool,
            'Mode': mode,
            'Total_Reads': flagstat.get('total_reads', 0),
            'Mapped_Reads': flagstat.get('mapped_reads', 0),
            'Mapping_Rate': flagstat.get('mapping_rate', 0),
            'Proper_Pairs': flagstat.get('proper_pairs', 0),
            'Pairing_Rate': flagstat.get('pairing_rate', 0),
            'Singletons': flagstat.get('singletons', 0),
            'Singleton_Rate': flagstat.get('singleton_rate', 0),
            'Mean_Coverage': coverage.get('mean_coverage', 0),
            'Coverage_Rate': coverage.get('coverage_rate', 0)
        }
        
        comparison_data.append(row)
    
    return pd.DataFrame(comparison_data)

def generate_statistics_report(all_stats):
    """生成统计报告"""
    report_file = "results/alignment_statistics.txt"
    
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write("序列比对统计分析报告\n")
        f.write("=" * 50 + "\n")
        f.write(f"生成时间: {pd.Timestamp.now()}\n\n")
        
        # 创建比较表
        df = create_comparison_table(all_stats)
        
        if not df.empty:
            f.write("1. 比对工具性能比较\n")
            f.write("-" * 30 + "\n")
            
            # 按工具分组统计
            for tool in df['Tool'].unique():
                tool_data = df[df['Tool'] == tool]
                f.write(f"\n{tool} 工具结果:\n")
                
                for _, row in tool_data.iterrows():
                    f.write(f"  {row['Mode']} 模式:\n")
                    f.write(f"    总reads数: {row['Total_Reads']:,}\n")
                    f.write(f"    比对率: {row['Mapping_Rate']:.2f}%\n")
                    f.write(f"    配对率: {row['Pairing_Rate']:.2f}%\n")
                    f.write(f"    平均覆盖度: {row['Mean_Coverage']:.2f}X\n")
                    f.write(f"    覆盖率: {row['Coverage_Rate']:.2f}%\n")
        
        f.write("\n2. 详细统计信息\n")
        f.write("-" * 30 + "\n")
        
        for sample_name, stats in all_stats.items():
            f.write(f"\n样本: {sample_name}\n")
            
            # Flagstat统计
            flagstat = stats['flagstat']
            f.write("  基本统计:\n")
            for key, value in flagstat.items():
                if isinstance(value, float):
                    f.write(f"    {key}: {value:.2f}\n")
                else:
                    f.write(f"    {key}: {value:,}\n")
            
            # MAPQ分布
            mapq_dist = stats['mapq_distribution']
            if mapq_dist:
                f.write("  MAPQ分布 (前10个):\n")
                sorted_mapq = sorted(mapq_dist.items(), key=lambda x: x[1], reverse=True)[:10]
                for mapq, count in sorted_mapq:
                    f.write(f"    MAPQ {mapq}: {count:,} reads\n")
            
            # 插入片段统计
            insert_sizes = stats['insert_sizes']
            if insert_sizes:
                f.write("  插入片段统计:\n")
                f.write(f"    平均值: {np.mean(insert_sizes):.2f}\n")
                f.write(f"    中位数: {np.median(insert_sizes):.2f}\n")
                f.write(f"    标准差: {np.std(insert_sizes):.2f}\n")
            
            # 覆盖度统计
            coverage = stats['coverage']
            if coverage:
                f.write("  覆盖度统计:\n")
                for key, value in coverage.items():
                    if isinstance(value, float):
                        f.write(f"    {key}: {value:.2f}\n")
                    else:
                        f.write(f"    {key}: {value:,}\n")
    
    print(f"统计报告已生成: {report_file}")
    return report_file

def create_visualization_plots(all_stats):
    """创建可视化图表"""
    # 设置图表样式
    plt.style.use('default')
    sns.set_palette("husl")
    
    # 创建比较数据框
    df = create_comparison_table(all_stats)
    
    if df.empty:
        print("没有找到比对统计数据")
        return
    
    # 1. 比对率比较
    plt.figure(figsize=(12, 8))
    
    plt.subplot(2, 2, 1)
    sns.barplot(data=df, x='Tool', y='Mapping_Rate', hue='Mode')
    plt.title('比对率比较')
    plt.ylabel('比对率 (%)')
    plt.xticks(rotation=45)
    
    # 2. 配对率比较
    plt.subplot(2, 2, 2)
    sns.barplot(data=df, x='Tool', y='Pairing_Rate', hue='Mode')
    plt.title('配对率比较')
    plt.ylabel('配对率 (%)')
    plt.xticks(rotation=45)
    
    # 3. 覆盖度比较
    plt.subplot(2, 2, 3)
    sns.barplot(data=df, x='Tool', y='Mean_Coverage', hue='Mode')
    plt.title('平均覆盖度比较')
    plt.ylabel('平均覆盖度 (X)')
    plt.xticks(rotation=45)
    
    # 4. 覆盖率比较
    plt.subplot(2, 2, 4)
    sns.barplot(data=df, x='Tool', y='Coverage_Rate', hue='Mode')
    plt.title('覆盖率比较')
    plt.ylabel('覆盖率 (%)')
    plt.xticks(rotation=45)
    
    plt.tight_layout()
    plt.savefig('results/alignment_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. MAPQ分布图
    plt.figure(figsize=(15, 10))
    
    # 收集所有MAPQ数据
    mapq_data = []
    for sample_name, stats in all_stats.items():
        mapq_dist = stats['mapq_distribution']
        for mapq, count in mapq_dist.items():
            mapq_data.extend([mapq] * count)
            if len(mapq_data) > 10000:  # 限制数据量
                break
    
    if mapq_data:
        plt.hist(mapq_data, bins=50, alpha=0.7, edgecolor='black')
        plt.title('MAPQ分布')
        plt.xlabel('MAPQ值')
        plt.ylabel('Reads数量')
        plt.grid(True, alpha=0.3)
        
        plt.savefig('results/mapq_distribution.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    # 3. 插入片段大小分布
    plt.figure(figsize=(12, 8))
    
    insert_data = {}
    for sample_name, stats in all_stats.items():
        insert_sizes = stats['insert_sizes']
        if insert_sizes and len(insert_sizes) > 100:  # 只显示有足够数据的样本
            insert_data[sample_name] = insert_sizes[:1000]  # 限制数据量
    
    if insert_data:
        for i, (sample_name, sizes) in enumerate(insert_data.items()):
            plt.hist(sizes, bins=50, alpha=0.6, label=sample_name[:20])
        
        plt.title('插入片段大小分布')
        plt.xlabel('插入片段大小 (bp)')
        plt.ylabel('频率')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plt.savefig('results/insert_size_distribution.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    print("可视化图表已生成:")
    print("  - results/alignment_comparison.png")
    print("  - results/mapq_distribution.png")
    print("  - results/insert_size_distribution.png")

def main():
    """主函数"""
    print("开始比对统计分析...")
    
    # 检查结果目录
    if not os.path.exists("results"):
        print("错误：results目录不存在")
        sys.exit(1)
    
    # 收集统计信息
    print("收集比对统计信息...")
    all_stats = collect_alignment_statistics()
    
    if not all_stats:
        print("未找到比对结果文件")
        sys.exit(1)
    
    print(f"找到 {len(all_stats)} 个比对结果")
    
    # 生成统计报告
    print("生成统计报告...")
    report_file = generate_statistics_report(all_stats)
    
    # 创建可视化图表
    print("创建可视化图表...")
    try:
        create_visualization_plots(all_stats)
    except Exception as e:
        print(f"创建图表时出错: {e}")
        print("继续执行其他分析...")
    
    # 创建CSV格式的比较表
    print("生成CSV比较表...")
    df = create_comparison_table(all_stats)
    if not df.empty:
        csv_file = "results/alignment_comparison.csv"
        df.to_csv(csv_file, index=False, encoding='utf-8')
        print(f"CSV比较表已生成: {csv_file}")
    
    print("\n比对统计分析完成！")
    print("生成的文件:")
    print(f"  - {report_file}")
    print("  - results/alignment_comparison.csv")
    print("  - results/alignment_comparison.png")
    print("  - results/mapq_distribution.png")
    print("  - results/insert_size_distribution.png")

if __name__ == "__main__":
    main()