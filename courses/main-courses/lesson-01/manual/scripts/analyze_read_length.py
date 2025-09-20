#!/usr/bin/env python3
"""
读长分析脚本
分析不同测序平台的读长分布特征
作者：王运生
日期：2025年
"""

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
import numpy as np

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

def parse_fastq(filename):
    """解析FASTQ文件，返回序列长度列表"""
    lengths = []
    try:
        with open(filename, 'r') as handle:
            for record in SeqIO.parse(handle, "fastq"):
                lengths.append(len(record.seq))
    except FileNotFoundError:
        print(f"警告：文件 {filename} 不存在，使用模拟数据")
        # 生成模拟数据
        if 'illumina' in filename:
            lengths = np.random.normal(150, 10, 10000).astype(int)
            lengths = lengths[lengths > 0]
        elif 'pacbio' in filename:
            lengths = np.random.lognormal(9, 0.5, 5000).astype(int)
        elif 'nanopore' in filename:
            lengths = np.random.lognormal(8.5, 0.8, 3000).astype(int)
    
    return lengths

def analyze_platform_data():
    """分析各平台数据"""
    platforms = {
        'Illumina': '../data/illumina_sample.fastq',
        'PacBio': '../data/pacbio_sample.fastq', 
        'Nanopore': '../data/nanopore_sample.fastq'
    }
    
    results = {}
    
    for platform, filename in platforms.items():
        print(f"分析 {platform} 数据...")
        lengths = parse_fastq(filename)
        
        if lengths:
            results[platform] = {
                'lengths': lengths,
                'count': len(lengths),
                'mean': np.mean(lengths),
                'median': np.median(lengths),
                'std': np.std(lengths),
                'min': np.min(lengths),
                'max': np.max(lengths)
            }
            
            print(f"{platform} 统计结果:")
            print(f"  序列数量: {results[platform]['count']:,}")
            print(f"  平均长度: {results[platform]['mean']:.1f} bp")
            print(f"  中位数长度: {results[platform]['median']:.1f} bp")
            print(f"  标准差: {results[platform]['std']:.1f}")
            print(f"  长度范围: {results[platform]['min']}-{results[platform]['max']} bp")
            print()
    
    return results

def create_length_distribution_plot(results):
    """创建读长分布图"""
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('测序平台读长分布分析', fontsize=16, fontweight='bold')
    
    # 子图1：读长分布直方图
    ax1 = axes[0, 0]
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
    
    for i, (platform, data) in enumerate(results.items()):
        ax1.hist(data['lengths'], bins=50, alpha=0.7, 
                label=platform, color=colors[i], density=True)
    
    ax1.set_xlabel('读长 (bp)')
    ax1.set_ylabel('密度')
    ax1.set_title('读长分布直方图')
    ax1.legend()
    ax1.set_xscale('log')
    
    # 子图2：箱线图
    ax2 = axes[0, 1]
    box_data = [data['lengths'] for data in results.values()]
    box_labels = list(results.keys())
    
    bp = ax2.boxplot(box_data, labels=box_labels, patch_artist=True)
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    ax2.set_ylabel('读长 (bp)')
    ax2.set_title('读长分布箱线图')
    ax2.set_yscale('log')
    
    # 子图3：累积分布
    ax3 = axes[1, 0]
    for i, (platform, data) in enumerate(results.items()):
        sorted_lengths = np.sort(data['lengths'])
        cumulative = np.arange(1, len(sorted_lengths) + 1) / len(sorted_lengths)
        ax3.plot(sorted_lengths, cumulative, label=platform, 
                color=colors[i], linewidth=2)
    
    ax3.set_xlabel('读长 (bp)')
    ax3.set_ylabel('累积概率')
    ax3.set_title('读长累积分布')
    ax3.legend()
    ax3.set_xscale('log')
    ax3.grid(True, alpha=0.3)
    
    # 子图4：统计摘要
    ax4 = axes[1, 1]
    ax4.axis('off')
    
    summary_text = "平台统计摘要\n\n"
    for platform, data in results.items():
        summary_text += f"{platform}:\n"
        summary_text += f"  平均长度: {data['mean']:.0f} bp\n"
        summary_text += f"  中位数: {data['median']:.0f} bp\n"
        summary_text += f"  序列数: {data['count']:,}\n\n"
    
    ax4.text(0.1, 0.9, summary_text, transform=ax4.transAxes, 
            fontsize=12, verticalalignment='top', fontfamily='monospace')
    
    plt.tight_layout()
    
    # 保存图片
    os.makedirs('../results', exist_ok=True)
    plt.savefig('../results/read_length_analysis.png', dpi=300, bbox_inches='tight')
    print("读长分析图已保存到: ../results/read_length_analysis.png")
    
    return fig

def save_statistics(results):
    """保存统计结果到CSV文件"""
    stats_data = []
    
    for platform, data in results.items():
        stats_data.append({
            '平台': platform,
            '序列数量': data['count'],
            '平均长度(bp)': round(data['mean'], 1),
            '中位数长度(bp)': round(data['median'], 1),
            '标准差': round(data['std'], 1),
            '最小长度(bp)': data['min'],
            '最大长度(bp)': data['max']
        })
    
    df = pd.DataFrame(stats_data)
    
    os.makedirs('../results', exist_ok=True)
    df.to_csv('../results/read_length_stats.csv', index=False, encoding='utf-8-sig')
    print("统计结果已保存到: ../results/read_length_stats.csv")
    
    return df

def main():
    """主函数"""
    print("开始分析测序平台读长分布...")
    print("=" * 50)
    
    # 分析数据
    results = analyze_platform_data()
    
    if not results:
        print("错误：没有找到有效的数据文件")
        return
    
    # 创建可视化图表
    create_length_distribution_plot(results)
    
    # 保存统计结果
    save_statistics(results)
    
    print("=" * 50)
    print("读长分析完成！")
    print("结果文件:")
    print("  - 分析图表: ../results/read_length_analysis.png")
    print("  - 统计数据: ../results/read_length_stats.csv")

if __name__ == "__main__":
    main()