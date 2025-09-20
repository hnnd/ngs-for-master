#!/usr/bin/env python3
"""
质量分析脚本
分析不同测序平台的质量分布特征
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

def parse_quality_scores(filename):
    """解析FASTQ文件，返回质量分数列表"""
    quality_scores = []
    
    try:
        with open(filename, 'r') as handle:
            for record in SeqIO.parse(handle, "fastq"):
                # 获取质量分数（Phred+33编码）
                scores = record.letter_annotations["phred_quality"]
                quality_scores.extend(scores)
    except FileNotFoundError:
        print(f"警告：文件 {filename} 不存在，使用模拟数据")
        # 生成模拟质量数据
        if 'illumina' in filename:
            # Illumina通常质量较高
            quality_scores = np.random.normal(35, 5, 100000).astype(int)
            quality_scores = np.clip(quality_scores, 2, 40)
        elif 'pacbio' in filename:
            # PacBio质量中等
            quality_scores = np.random.normal(20, 8, 50000).astype(int)
            quality_scores = np.clip(quality_scores, 2, 35)
        elif 'nanopore' in filename:
            # Nanopore质量较低但在改善
            quality_scores = np.random.normal(15, 6, 30000).astype(int)
            quality_scores = np.clip(quality_scores, 2, 30)
    
    return quality_scores

def calculate_quality_metrics(quality_scores):
    """计算质量指标"""
    if not quality_scores:
        return {}
    
    scores_array = np.array(quality_scores)
    
    metrics = {
        'mean_quality': np.mean(scores_array),
        'median_quality': np.median(scores_array),
        'q20_percent': np.sum(scores_array >= 20) / len(scores_array) * 100,
        'q30_percent': np.sum(scores_array >= 30) / len(scores_array) * 100,
        'q40_percent': np.sum(scores_array >= 40) / len(scores_array) * 100,
        'low_quality_percent': np.sum(scores_array < 10) / len(scores_array) * 100,
        'total_bases': len(scores_array)
    }
    
    return metrics

def analyze_platform_quality():
    """分析各平台质量数据"""
    platforms = {
        'Illumina': '../data/illumina_sample.fastq',
        'PacBio': '../data/pacbio_sample.fastq',
        'Nanopore': '../data/nanopore_sample.fastq'
    }
    
    results = {}
    
    for platform, filename in platforms.items():
        print(f"分析 {platform} 质量数据...")
        quality_scores = parse_quality_scores(filename)
        
        if quality_scores:
            metrics = calculate_quality_metrics(quality_scores)
            results[platform] = {
                'quality_scores': quality_scores,
                'metrics': metrics
            }
            
            print(f"{platform} 质量统计:")
            print(f"  总碱基数: {metrics['total_bases']:,}")
            print(f"  平均质量: {metrics['mean_quality']:.1f}")
            print(f"  中位数质量: {metrics['median_quality']:.1f}")
            print(f"  Q20以上: {metrics['q20_percent']:.1f}%")
            print(f"  Q30以上: {metrics['q30_percent']:.1f}%")
            print(f"  低质量(<Q10): {metrics['low_quality_percent']:.1f}%")
            print()
    
    return results

def create_quality_plots(results):
    """创建质量分析图表"""
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('测序平台质量分析', fontsize=16, fontweight='bold')
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
    
    # 子图1：质量分布直方图
    ax1 = axes[0, 0]
    for i, (platform, data) in enumerate(results.items()):
        # 采样以减少计算量
        sample_scores = np.random.choice(data['quality_scores'], 
                                       min(10000, len(data['quality_scores'])), 
                                       replace=False)
        ax1.hist(sample_scores, bins=range(0, 45), alpha=0.7, 
                label=platform, color=colors[i], density=True)
    
    ax1.set_xlabel('质量分数 (Phred)')
    ax1.set_ylabel('密度')
    ax1.set_title('质量分数分布')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 子图2：质量指标比较
    ax2 = axes[0, 1]
    platforms = list(results.keys())
    q20_values = [results[p]['metrics']['q20_percent'] for p in platforms]
    q30_values = [results[p]['metrics']['q30_percent'] for p in platforms]
    
    x = np.arange(len(platforms))
    width = 0.35
    
    ax2.bar(x - width/2, q20_values, width, label='Q20+', alpha=0.8, color='skyblue')
    ax2.bar(x + width/2, q30_values, width, label='Q30+', alpha=0.8, color='lightcoral')
    
    ax2.set_xlabel('测序平台')
    ax2.set_ylabel('百分比 (%)')
    ax2.set_title('高质量碱基比例')
    ax2.set_xticks(x)
    ax2.set_xticklabels(platforms)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 子图3：平均质量比较
    ax3 = axes[1, 0]
    mean_qualities = [results[p]['metrics']['mean_quality'] for p in platforms]
    
    bars = ax3.bar(platforms, mean_qualities, color=colors, alpha=0.8)
    ax3.set_ylabel('平均质量分数')
    ax3.set_title('平台平均质量比较')
    ax3.grid(True, alpha=0.3)
    
    # 在柱状图上添加数值标签
    for bar, value in zip(bars, mean_qualities):
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                f'{value:.1f}', ha='center', va='bottom')
    
    # 子图4：质量统计摘要
    ax4 = axes[1, 1]
    ax4.axis('off')
    
    summary_text = "质量统计摘要\n\n"
    for platform, data in results.items():
        metrics = data['metrics']
        summary_text += f"{platform}:\n"
        summary_text += f"  平均质量: {metrics['mean_quality']:.1f}\n"
        summary_text += f"  Q30+: {metrics['q30_percent']:.1f}%\n"
        summary_text += f"  总碱基: {metrics['total_bases']:,}\n\n"
    
    ax4.text(0.1, 0.9, summary_text, transform=ax4.transAxes,
            fontsize=12, verticalalignment='top', fontfamily='monospace')
    
    plt.tight_layout()
    
    # 保存图片
    os.makedirs('../results', exist_ok=True)
    plt.savefig('../results/quality_analysis.png', dpi=300, bbox_inches='tight')
    print("质量分析图已保存到: ../results/quality_analysis.png")
    
    return fig

def save_quality_metrics(results):
    """保存质量指标到CSV文件"""
    metrics_data = []
    
    for platform, data in results.items():
        metrics = data['metrics']
        metrics_data.append({
            '平台': platform,
            '总碱基数': metrics['total_bases'],
            '平均质量': round(metrics['mean_quality'], 2),
            '中位数质量': round(metrics['median_quality'], 2),
            'Q20以上(%)': round(metrics['q20_percent'], 2),
            'Q30以上(%)': round(metrics['q30_percent'], 2),
            'Q40以上(%)': round(metrics['q40_percent'], 2),
            '低质量<Q10(%)': round(metrics['low_quality_percent'], 2)
        })
    
    df = pd.DataFrame(metrics_data)
    
    os.makedirs('../results', exist_ok=True)
    df.to_csv('../results/quality_metrics.csv', index=False, encoding='utf-8-sig')
    print("质量指标已保存到: ../results/quality_metrics.csv")
    
    return df

def main():
    """主函数"""
    print("开始分析测序平台质量分布...")
    print("=" * 50)
    
    # 分析质量数据
    results = analyze_platform_quality()
    
    if not results:
        print("错误：没有找到有效的数据文件")
        return
    
    # 创建质量分析图表
    create_quality_plots(results)
    
    # 保存质量指标
    save_quality_metrics(results)
    
    print("=" * 50)
    print("质量分析完成！")
    print("结果文件:")
    print("  - 分析图表: ../results/quality_analysis.png")
    print("  - 质量指标: ../results/quality_metrics.csv")

if __name__ == "__main__":
    main()