#!/usr/bin/env python3
"""
测序平台比较分析脚本
比较不同测序平台的技术特征
作者：王运生
日期：2025年
"""

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import argparse
from Bio import SeqIO

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

def analyze_platform_characteristics():
    """分析各平台的技术特征"""
    platforms = {
        'Illumina': '../data/illumina_sample.fastq',
        'PacBio': '../data/pacbio_sample.fastq',
        'Nanopore': '../data/nanopore_sample.fastq'
    }
    
    results = {}
    
    for platform, filename in platforms.items():
        print(f"分析 {platform} 平台特征...")
        
        try:
            # 读取FASTQ文件
            sequences = []
            qualities = []
            
            with open(filename, 'r') as handle:
                for record in SeqIO.parse(handle, "fastq"):
                    sequences.append(len(record.seq))
                    qualities.extend(record.letter_annotations["phred_quality"])
            
            # 计算统计指标
            results[platform] = {
                'read_count': len(sequences),
                'total_bases': sum(sequences),
                'mean_length': np.mean(sequences),
                'median_length': np.median(sequences),
                'min_length': np.min(sequences),
                'max_length': np.max(sequences),
                'length_std': np.std(sequences),
                'mean_quality': np.mean(qualities),
                'median_quality': np.median(qualities),
                'q20_percent': np.sum(np.array(qualities) >= 20) / len(qualities) * 100,
                'q30_percent': np.sum(np.array(qualities) >= 30) / len(qualities) * 100,
                'sequences': sequences[:1000],  # 保存前1000个用于可视化
                'qualities': qualities[:10000]   # 保存前10000个质量值
            }
            
        except FileNotFoundError:
            print(f"警告：文件 {filename} 不存在，使用模拟数据")
            # 生成模拟数据
            if platform == 'Illumina':
                sequences = list(np.random.normal(150, 10, 10000).astype(int))
                sequences = [max(50, min(300, x)) for x in sequences]
                qualities = list(np.random.normal(35, 5, 100000).astype(int))
                qualities = [max(2, min(40, x)) for x in qualities]
            elif platform == 'PacBio':
                sequences = list(np.random.lognormal(9, 0.5, 1000).astype(int))
                qualities = list(np.random.normal(20, 8, 50000).astype(int))
                qualities = [max(2, min(35, x)) for x in qualities]
            else:  # Nanopore
                sequences = list(np.random.lognormal(8.5, 0.8, 500).astype(int))
                qualities = list(np.random.normal(15, 6, 30000).astype(int))
                qualities = [max(2, min(30, x)) for x in qualities]
            
            results[platform] = {
                'read_count': len(sequences),
                'total_bases': sum(sequences),
                'mean_length': np.mean(sequences),
                'median_length': np.median(sequences),
                'min_length': np.min(sequences),
                'max_length': np.max(sequences),
                'length_std': np.std(sequences),
                'mean_quality': np.mean(qualities),
                'median_quality': np.median(qualities),
                'q20_percent': np.sum(np.array(qualities) >= 20) / len(qualities) * 100,
                'q30_percent': np.sum(np.array(qualities) >= 30) / len(qualities) * 100,
                'sequences': sequences,
                'qualities': qualities
            }
    
    return results

def create_comprehensive_comparison(results, analysis_type='read_length'):
    """创建综合比较图表"""
    if analysis_type == 'read_length':
        create_read_length_comparison(results)
    elif analysis_type == 'quality':
        create_quality_comparison(results)
    elif analysis_type == 'throughput':
        create_throughput_comparison(results)
    else:
        create_all_comparisons(results)

def create_read_length_comparison(results):
    """创建读长比较图表"""
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('测序平台读长特征比较', fontsize=16, fontweight='bold')
    
    platforms = list(results.keys())
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
    
    # 子图1：读长分布对比
    ax1 = axes[0, 0]
    for i, (platform, data) in enumerate(results.items()):
        ax1.hist(data['sequences'], bins=50, alpha=0.7, 
                label=platform, color=colors[i], density=True)
    ax1.set_xlabel('读长 (bp)')
    ax1.set_ylabel('密度')
    ax1.set_title('读长分布对比')
    ax1.legend()
    ax1.set_xscale('log')
    
    # 子图2：读长统计比较
    ax2 = axes[0, 1]
    metrics = ['mean_length', 'median_length', 'max_length']
    metric_names = ['平均长度', '中位数长度', '最大长度']
    
    x = np.arange(len(platforms))
    width = 0.25
    
    for i, metric in enumerate(metrics):
        values = [results[p][metric] for p in platforms]
        ax2.bar(x + i * width, values, width, label=metric_names[i], alpha=0.8)
    
    ax2.set_xlabel('测序平台')
    ax2.set_ylabel('读长 (bp)')
    ax2.set_title('读长统计指标比较')
    ax2.set_xticks(x + width)
    ax2.set_xticklabels(platforms)
    ax2.legend()
    ax2.set_yscale('log')
    
    # 子图3：读长变异系数
    ax3 = axes[1, 0]
    cv_values = [results[p]['length_std'] / results[p]['mean_length'] for p in platforms]
    bars = ax3.bar(platforms, cv_values, color=colors, alpha=0.8)
    ax3.set_ylabel('变异系数')
    ax3.set_title('读长变异性比较')
    
    for bar, value in zip(bars, cv_values):
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                f'{value:.2f}', ha='center', va='bottom')
    
    # 子图4：数据产出比较
    ax4 = axes[1, 1]
    throughput_gb = [results[p]['total_bases'] / 1e9 for p in platforms]
    bars = ax4.bar(platforms, throughput_gb, color=colors, alpha=0.8)
    ax4.set_ylabel('数据产出 (GB)')
    ax4.set_title('数据产出量比较')
    
    for bar, value in zip(bars, throughput_gb):
        ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                f'{value:.2f}', ha='center', va='bottom')
    
    plt.tight_layout()
    
    os.makedirs('../results', exist_ok=True)
    plt.savefig('../results/platform_read_length_comparison.png', dpi=300, bbox_inches='tight')
    print("读长比较图已保存到: ../results/platform_read_length_comparison.png")

def create_quality_comparison(results):
    """创建质量比较图表"""
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('测序平台质量特征比较', fontsize=16, fontweight='bold')
    
    platforms = list(results.keys())
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
    
    # 子图1：质量分布对比
    ax1 = axes[0, 0]
    for i, (platform, data) in enumerate(results.items()):
        sample_qualities = np.random.choice(data['qualities'], 
                                          min(5000, len(data['qualities'])), 
                                          replace=False)
        ax1.hist(sample_qualities, bins=range(0, 45), alpha=0.7,
                label=platform, color=colors[i], density=True)
    ax1.set_xlabel('质量分数 (Phred)')
    ax1.set_ylabel('密度')
    ax1.set_title('质量分数分布对比')
    ax1.legend()
    
    # 子图2：高质量碱基比例
    ax2 = axes[0, 1]
    q20_values = [results[p]['q20_percent'] for p in platforms]
    q30_values = [results[p]['q30_percent'] for p in platforms]
    
    x = np.arange(len(platforms))
    width = 0.35
    
    ax2.bar(x - width/2, q20_values, width, label='Q20+', alpha=0.8, color='skyblue')
    ax2.bar(x + width/2, q30_values, width, label='Q30+', alpha=0.8, color='lightcoral')
    
    ax2.set_xlabel('测序平台')
    ax2.set_ylabel('百分比 (%)')
    ax2.set_title('高质量碱基比例比较')
    ax2.set_xticks(x)
    ax2.set_xticklabels(platforms)
    ax2.legend()
    
    # 子图3：平均质量比较
    ax3 = axes[1, 0]
    mean_qualities = [results[p]['mean_quality'] for p in platforms]
    bars = ax3.bar(platforms, mean_qualities, color=colors, alpha=0.8)
    ax3.set_ylabel('平均质量分数')
    ax3.set_title('平均质量比较')
    
    for bar, value in zip(bars, mean_qualities):
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                f'{value:.1f}', ha='center', va='bottom')
    
    # 子图4：质量-长度关系（散点图）
    ax4 = axes[1, 1]
    for i, (platform, data) in enumerate(results.items()):
        # 采样数据点
        sample_size = min(1000, len(data['sequences']))
        sample_indices = np.random.choice(len(data['sequences']), sample_size, replace=False)
        
        lengths = [data['sequences'][i] for i in sample_indices]
        # 为每个序列计算平均质量（简化处理）
        avg_qualities = [results[platform]['mean_quality'] + np.random.normal(0, 2) 
                        for _ in range(sample_size)]
        
        ax4.scatter(lengths, avg_qualities, alpha=0.6, label=platform, 
                   color=colors[i], s=20)
    
    ax4.set_xlabel('读长 (bp)')
    ax4.set_ylabel('平均质量分数')
    ax4.set_title('读长-质量关系')
    ax4.set_xscale('log')
    ax4.legend()
    
    plt.tight_layout()
    plt.savefig('../results/platform_quality_comparison.png', dpi=300, bbox_inches='tight')
    print("质量比较图已保存到: ../results/platform_quality_comparison.png")

def save_comparison_summary(results):
    """保存比较摘要"""
    summary_data = []
    
    for platform, data in results.items():
        summary_data.append({
            '平台': platform,
            '序列数量': data['read_count'],
            '总碱基数': data['total_bases'],
            '平均读长(bp)': round(data['mean_length'], 1),
            '中位数读长(bp)': round(data['median_length'], 1),
            '最大读长(bp)': data['max_length'],
            '读长标准差': round(data['length_std'], 1),
            '平均质量': round(data['mean_quality'], 1),
            'Q20以上(%)': round(data['q20_percent'], 1),
            'Q30以上(%)': round(data['q30_percent'], 1)
        })
    
    df = pd.DataFrame(summary_data)
    
    os.makedirs('../results', exist_ok=True)
    df.to_csv('../results/platform_comparison_summary.csv', index=False, encoding='utf-8-sig')
    print("比较摘要已保存到: ../results/platform_comparison_summary.csv")
    
    return df

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='测序平台比较分析')
    parser.add_argument('--analysis', type=str, default='all',
                       choices=['read_length', 'quality', 'throughput', 'all'],
                       help='分析类型')
    
    args = parser.parse_args()
    
    print("测序平台比较分析")
    print("=" * 50)
    
    # 分析平台特征
    results = analyze_platform_characteristics()
    
    # 显示基本统计
    print("\n平台基本统计:")
    for platform, data in results.items():
        print(f"\n{platform}:")
        print(f"  序列数量: {data['read_count']:,}")
        print(f"  平均读长: {data['mean_length']:.1f} bp")
        print(f"  平均质量: {data['mean_quality']:.1f}")
        print(f"  Q30以上: {data['q30_percent']:.1f}%")
    
    # 创建比较图表
    if args.analysis == 'all':
        create_read_length_comparison(results)
        create_quality_comparison(results)
    else:
        create_comprehensive_comparison(results, args.analysis)
    
    # 保存摘要
    save_comparison_summary(results)
    
    print("\n" + "=" * 50)
    print("平台比较分析完成！")
    print("结果文件:")
    print("  - 读长比较: ../results/platform_read_length_comparison.png")
    print("  - 质量比较: ../results/platform_quality_comparison.png")
    print("  - 比较摘要: ../results/platform_comparison_summary.csv")

if __name__ == "__main__":
    main()