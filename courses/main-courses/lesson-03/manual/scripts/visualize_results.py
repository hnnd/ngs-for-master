#!/usr/bin/env python3
"""
比对结果可视化脚本
作者：王运生
日期：2025-01-20
用法：python visualize_results.py
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
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体和样式
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

# 设置颜色主题
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f']
sns.set_palette(colors)

def load_alignment_data():
    """加载比对数据"""
    results_dir = Path("results")
    
    # 查找比较表文件
    comparison_file = results_dir / "alignment_comparison.csv"
    
    if comparison_file.exists():
        df = pd.read_csv(comparison_file)
        return df
    else:
        print("未找到比较数据文件，尝试从统计文件重新生成...")
        return generate_comparison_data()

def generate_comparison_data():
    """从统计文件生成比较数据"""
    results_dir = Path("results")
    
    # 查找所有统计文件
    stats_files = list(results_dir.glob("*_stats.txt"))
    
    comparison_data = []
    
    for stats_file in stats_files:
        base_name = stats_file.stem.replace("_stats", "")
        
        # 解析统计文件
        stats = parse_flagstat_file(stats_file)
        
        if stats:
            # 确定工具和模式
            tool, mode = identify_tool_and_mode(base_name)
            
            row = {
                'Sample': base_name,
                'Tool': tool,
                'Mode': mode,
                'Total_Reads': stats.get('total_reads', 0),
                'Mapped_Reads': stats.get('mapped_reads', 0),
                'Mapping_Rate': stats.get('mapping_rate', 0),
                'Proper_Pairs': stats.get('proper_pairs', 0),
                'Pairing_Rate': stats.get('pairing_rate', 0),
                'Singletons': stats.get('singletons', 0),
                'Singleton_Rate': stats.get('singleton_rate', 0)
            }
            
            comparison_data.append(row)
    
    return pd.DataFrame(comparison_data)

def parse_flagstat_file(stats_file):
    """解析flagstat文件"""
    stats = {}
    
    try:
        with open(stats_file, 'r') as f:
            content = f.read()
        
        # 解析各种统计信息
        patterns = {
            'total_reads': r'(\d+) \+ \d+ in total',
            'mapped_reads': r'(\d+) \+ \d+ mapped \(',
            'paired_reads': r'(\d+) \+ \d+ paired in sequencing',
            'proper_pairs': r'(\d+) \+ \d+ properly paired',
            'singletons': r'(\d+) \+ \d+ singletons'
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
        
    except Exception as e:
        print(f"解析文件 {stats_file} 时出错: {e}")
    
    return stats

def identify_tool_and_mode(sample_name):
    """识别工具和模式"""
    sample_lower = sample_name.lower()
    
    if 'bwa' in sample_lower:
        tool = 'BWA'
        if 'default' in sample_lower:
            mode = 'Default'
        elif 'optimized' in sample_lower:
            mode = 'Optimized'
        elif 'strict' in sample_lower:
            mode = 'Strict'
        else:
            mode = 'Unknown'
    elif 'bowtie2' in sample_lower:
        tool = 'Bowtie2'
        if 'fast' in sample_lower:
            mode = 'Fast'
        elif 'very_sensitive' in sample_lower:
            mode = 'Very Sensitive'
        elif 'sensitive' in sample_lower and 'very' not in sample_lower:
            mode = 'Sensitive'
        elif 'local' in sample_lower:
            mode = 'Local'
        elif 'custom' in sample_lower:
            mode = 'Custom'
        else:
            mode = 'Unknown'
    else:
        tool = 'Unknown'
        mode = 'Unknown'
    
    return tool, mode

def create_performance_comparison(df):
    """创建性能比较图"""
    if df.empty:
        print("没有数据可供可视化")
        return
    
    # 创建子图
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('序列比对工具性能比较', fontsize=16, fontweight='bold')
    
    # 1. 比对率比较
    ax1 = axes[0, 0]
    if 'Tool' in df.columns and 'Mode' in df.columns:
        sns.barplot(data=df, x='Tool', y='Mapping_Rate', hue='Mode', ax=ax1)
    else:
        sns.barplot(data=df, x='Sample', y='Mapping_Rate', ax=ax1)
        ax1.tick_params(axis='x', rotation=45)
    ax1.set_title('比对率比较')
    ax1.set_ylabel('比对率 (%)')
    ax1.grid(True, alpha=0.3)
    
    # 2. 配对率比较
    ax2 = axes[0, 1]
    if 'Tool' in df.columns and 'Mode' in df.columns:
        sns.barplot(data=df, x='Tool', y='Pairing_Rate', hue='Mode', ax=ax2)
    else:
        sns.barplot(data=df, x='Sample', y='Pairing_Rate', ax=ax2)
        ax2.tick_params(axis='x', rotation=45)
    ax2.set_title('配对率比较')
    ax2.set_ylabel('配对率 (%)')
    ax2.grid(True, alpha=0.3)
    
    # 3. 单端比对率比较
    ax3 = axes[0, 2]
    if 'Tool' in df.columns and 'Mode' in df.columns:
        sns.barplot(data=df, x='Tool', y='Singleton_Rate', hue='Mode', ax=ax3)
    else:
        sns.barplot(data=df, x='Sample', y='Singleton_Rate', ax=ax3)
        ax3.tick_params(axis='x', rotation=45)
    ax3.set_title('单端比对率比较')
    ax3.set_ylabel('单端比对率 (%)')
    ax3.grid(True, alpha=0.3)
    
    # 4. 总reads数比较
    ax4 = axes[1, 0]
    if 'Tool' in df.columns and 'Mode' in df.columns:
        sns.barplot(data=df, x='Tool', y='Total_Reads', hue='Mode', ax=ax4)
    else:
        sns.barplot(data=df, x='Sample', y='Total_Reads', ax=ax4)
        ax4.tick_params(axis='x', rotation=45)
    ax4.set_title('总reads数比较')
    ax4.set_ylabel('总reads数')
    ax4.grid(True, alpha=0.3)
    
    # 5. 比对成功reads数比较
    ax5 = axes[1, 1]
    if 'Tool' in df.columns and 'Mode' in df.columns:
        sns.barplot(data=df, x='Tool', y='Mapped_Reads', hue='Mode', ax=ax5)
    else:
        sns.barplot(data=df, x='Sample', y='Mapped_Reads', ax=ax5)
        ax5.tick_params(axis='x', rotation=45)
    ax5.set_title('比对成功reads数比较')
    ax5.set_ylabel('比对成功reads数')
    ax5.grid(True, alpha=0.3)
    
    # 6. 正确配对reads数比较
    ax6 = axes[1, 2]
    if 'Tool' in df.columns and 'Mode' in df.columns:
        sns.barplot(data=df, x='Tool', y='Proper_Pairs', hue='Mode', ax=ax6)
    else:
        sns.barplot(data=df, x='Sample', y='Proper_Pairs', ax=ax6)
        ax6.tick_params(axis='x', rotation=45)
    ax6.set_title('正确配对reads数比较')
    ax6.set_ylabel('正确配对reads数')
    ax6.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('results/performance_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("性能比较图已生成: results/performance_comparison.png")

def create_mapq_analysis():
    """创建MAPQ分析图"""
    results_dir = Path("results")
    mapq_files = list(results_dir.glob("*_mapq.txt"))
    
    if not mapq_files:
        print("未找到MAPQ分布文件")
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('MAPQ质量分析', fontsize=16, fontweight='bold')
    
    # 收集所有MAPQ数据
    all_mapq_data = {}
    
    for mapq_file in mapq_files[:4]:  # 最多显示4个样本
        sample_name = mapq_file.stem.replace("_mapq", "")
        
        mapq_values = []
        mapq_counts = []
        
        try:
            with open(mapq_file, 'r') as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) == 2:
                        count = int(parts[0])
                        mapq = int(parts[1])
                        mapq_values.append(mapq)
                        mapq_counts.append(count)
            
            all_mapq_data[sample_name] = (mapq_values, mapq_counts)
        except Exception as e:
            print(f"读取MAPQ文件 {mapq_file} 时出错: {e}")
    
    # 绘制MAPQ分布
    for i, (sample_name, (mapq_values, mapq_counts)) in enumerate(all_mapq_data.items()):
        if i >= 4:
            break
        
        ax = axes[i // 2, i % 2]
        
        # 创建直方图数据
        mapq_expanded = []
        for mapq, count in zip(mapq_values, mapq_counts):
            mapq_expanded.extend([mapq] * min(count, 1000))  # 限制数据量
        
        if mapq_expanded:
            ax.hist(mapq_expanded, bins=50, alpha=0.7, edgecolor='black')
            ax.set_title(f'{sample_name} MAPQ分布')
            ax.set_xlabel('MAPQ值')
            ax.set_ylabel('Reads数量')
            ax.grid(True, alpha=0.3)
            
            # 添加统计信息
            mean_mapq = np.mean(mapq_expanded)
            ax.axvline(mean_mapq, color='red', linestyle='--', 
                      label=f'平均值: {mean_mapq:.1f}')
            ax.legend()
    
    plt.tight_layout()
    plt.savefig('results/mapq_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("MAPQ分析图已生成: results/mapq_analysis.png")

def create_coverage_analysis():
    """创建覆盖度分析图"""
    results_dir = Path("results")
    coverage_files = list(results_dir.glob("*_coverage.txt"))
    
    if not coverage_files:
        print("未找到覆盖度文件")
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('基因组覆盖度分析', fontsize=16, fontweight='bold')
    
    coverage_data = {}
    
    # 读取覆盖度数据（采样以减少内存使用）
    for cov_file in coverage_files[:4]:
        sample_name = cov_file.stem.replace("_coverage", "")
        
        coverages = []
        positions = []
        
        try:
            with open(cov_file, 'r') as f:
                for i, line in enumerate(f):
                    if i % 100 == 0:  # 每100行采样一次
                        parts = line.strip().split()
                        if len(parts) >= 3:
                            try:
                                pos = int(parts[1])
                                cov = int(parts[2])
                                positions.append(pos)
                                coverages.append(cov)
                            except ValueError:
                                continue
                    
                    if len(coverages) >= 10000:  # 限制数据量
                        break
            
            coverage_data[sample_name] = (positions, coverages)
        except Exception as e:
            print(f"读取覆盖度文件 {cov_file} 时出错: {e}")
    
    # 绘制覆盖度图
    for i, (sample_name, (positions, coverages)) in enumerate(coverage_data.items()):
        if i >= 4:
            break
        
        ax = axes[i // 2, i % 2]
        
        if positions and coverages:
            # 覆盖度分布直方图
            ax.hist(coverages, bins=50, alpha=0.7, edgecolor='black')
            ax.set_title(f'{sample_name} 覆盖度分布')
            ax.set_xlabel('覆盖度')
            ax.set_ylabel('位点数量')
            ax.grid(True, alpha=0.3)
            
            # 添加统计信息
            mean_cov = np.mean(coverages)
            median_cov = np.median(coverages)
            ax.axvline(mean_cov, color='red', linestyle='--', 
                      label=f'平均值: {mean_cov:.1f}')
            ax.axvline(median_cov, color='orange', linestyle='--', 
                      label=f'中位数: {median_cov:.1f}')
            ax.legend()
    
    plt.tight_layout()
    plt.savefig('results/coverage_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("覆盖度分析图已生成: results/coverage_analysis.png")

def create_insert_size_analysis():
    """创建插入片段大小分析图"""
    results_dir = Path("results")
    insert_files = list(results_dir.glob("*_insert_sizes.txt"))
    
    if not insert_files:
        print("未找到插入片段大小文件")
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('插入片段大小分析', fontsize=16, fontweight='bold')
    
    insert_data = {}
    
    # 读取插入片段数据
    for insert_file in insert_files[:4]:
        sample_name = insert_file.stem.replace("_insert_sizes", "")
        
        insert_sizes = []
        
        try:
            with open(insert_file, 'r') as f:
                for line in f:
                    try:
                        size = int(line.strip())
                        if 0 < size < 2000:  # 过滤异常值
                            insert_sizes.append(size)
                        
                        if len(insert_sizes) >= 10000:  # 限制数据量
                            break
                    except ValueError:
                        continue
            
            insert_data[sample_name] = insert_sizes
        except Exception as e:
            print(f"读取插入片段文件 {insert_file} 时出错: {e}")
    
    # 绘制插入片段分布图
    for i, (sample_name, insert_sizes) in enumerate(insert_data.items()):
        if i >= 4:
            break
        
        ax = axes[i // 2, i % 2]
        
        if insert_sizes:
            ax.hist(insert_sizes, bins=50, alpha=0.7, edgecolor='black')
            ax.set_title(f'{sample_name} 插入片段分布')
            ax.set_xlabel('插入片段大小 (bp)')
            ax.set_ylabel('频率')
            ax.grid(True, alpha=0.3)
            
            # 添加统计信息
            mean_size = np.mean(insert_sizes)
            median_size = np.median(insert_sizes)
            std_size = np.std(insert_sizes)
            
            ax.axvline(mean_size, color='red', linestyle='--', 
                      label=f'平均值: {mean_size:.0f}')
            ax.axvline(median_size, color='orange', linestyle='--', 
                      label=f'中位数: {median_size:.0f}')
            
            # 添加文本信息
            ax.text(0.7, 0.8, f'标准差: {std_size:.0f}', 
                   transform=ax.transAxes, fontsize=10,
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            
            ax.legend()
    
    plt.tight_layout()
    plt.savefig('results/insert_size_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("插入片段分析图已生成: results/insert_size_analysis.png")

def create_summary_dashboard():
    """创建总结仪表板"""
    # 加载比较数据
    df = load_alignment_data()
    
    if df.empty:
        print("无法创建总结仪表板：没有比较数据")
        return
    
    fig, axes = plt.subplots(2, 3, figsize=(20, 12))
    fig.suptitle('序列比对分析总结仪表板', fontsize=18, fontweight='bold')
    
    # 1. 工具性能雷达图（如果有多个工具）
    ax1 = axes[0, 0]
    if 'Tool' in df.columns and len(df['Tool'].unique()) > 1:
        # 计算每个工具的平均性能
        tool_performance = df.groupby('Tool').agg({
            'Mapping_Rate': 'mean',
            'Pairing_Rate': 'mean',
            'Singleton_Rate': 'mean'
        }).reset_index()
        
        x = range(len(tool_performance))
        width = 0.25
        
        ax1.bar([i - width for i in x], tool_performance['Mapping_Rate'], 
               width, label='比对率', alpha=0.8)
        ax1.bar(x, tool_performance['Pairing_Rate'], 
               width, label='配对率', alpha=0.8)
        ax1.bar([i + width for i in x], tool_performance['Singleton_Rate'], 
               width, label='单端比对率', alpha=0.8)
        
        ax1.set_xlabel('比对工具')
        ax1.set_ylabel('百分比 (%)')
        ax1.set_title('工具性能对比')
        ax1.set_xticks(x)
        ax1.set_xticklabels(tool_performance['Tool'])
        ax1.legend()
        ax1.grid(True, alpha=0.3)
    else:
        ax1.text(0.5, 0.5, '需要多个工具数据\n才能显示对比', 
                ha='center', va='center', transform=ax1.transAxes)
        ax1.set_title('工具性能对比')
    
    # 2. 比对质量分布
    ax2 = axes[0, 1]
    mapping_rates = df['Mapping_Rate'].dropna()
    if len(mapping_rates) > 0:
        ax2.hist(mapping_rates, bins=10, alpha=0.7, edgecolor='black')
        ax2.axvline(mapping_rates.mean(), color='red', linestyle='--', 
                   label=f'平均值: {mapping_rates.mean():.1f}%')
        ax2.set_xlabel('比对率 (%)')
        ax2.set_ylabel('样本数量')
        ax2.set_title('比对率分布')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
    
    # 3. 配对质量分布
    ax3 = axes[0, 2]
    pairing_rates = df['Pairing_Rate'].dropna()
    if len(pairing_rates) > 0:
        ax3.hist(pairing_rates, bins=10, alpha=0.7, edgecolor='black')
        ax3.axvline(pairing_rates.mean(), color='red', linestyle='--', 
                   label=f'平均值: {pairing_rates.mean():.1f}%')
        ax3.set_xlabel('配对率 (%)')
        ax3.set_ylabel('样本数量')
        ax3.set_title('配对率分布')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
    
    # 4. 样本比较热图
    ax4 = axes[1, 0]
    if len(df) > 1:
        # 选择数值列创建热图
        numeric_cols = ['Mapping_Rate', 'Pairing_Rate', 'Singleton_Rate']
        heatmap_data = df[numeric_cols].T
        
        im = ax4.imshow(heatmap_data, cmap='RdYlBu_r', aspect='auto')
        ax4.set_xticks(range(len(df)))
        ax4.set_xticklabels([s[:15] + '...' if len(s) > 15 else s 
                            for s in df['Sample']], rotation=45)
        ax4.set_yticks(range(len(numeric_cols)))
        ax4.set_yticklabels(['比对率', '配对率', '单端比对率'])
        ax4.set_title('样本性能热图')
        
        # 添加颜色条
        cbar = plt.colorbar(im, ax=ax4)
        cbar.set_label('百分比 (%)')
    else:
        ax4.text(0.5, 0.5, '需要多个样本\n才能显示热图', 
                ha='center', va='center', transform=ax4.transAxes)
        ax4.set_title('样本性能热图')
    
    # 5. 数据质量评估
    ax5 = axes[1, 1]
    quality_scores = []
    quality_labels = []
    
    for _, row in df.iterrows():
        # 计算综合质量分数
        score = (row['Mapping_Rate'] * 0.4 + 
                row['Pairing_Rate'] * 0.4 + 
                (100 - row['Singleton_Rate']) * 0.2)
        quality_scores.append(score)
        quality_labels.append(row['Sample'][:10])
    
    if quality_scores:
        bars = ax5.bar(range(len(quality_scores)), quality_scores, alpha=0.7)
        ax5.set_xlabel('样本')
        ax5.set_ylabel('质量分数')
        ax5.set_title('综合质量评估')
        ax5.set_xticks(range(len(quality_labels)))
        ax5.set_xticklabels(quality_labels, rotation=45)
        ax5.grid(True, alpha=0.3)
        
        # 颜色编码
        for i, bar in enumerate(bars):
            if quality_scores[i] >= 90:
                bar.set_color('green')
            elif quality_scores[i] >= 80:
                bar.set_color('orange')
            else:
                bar.set_color('red')
    
    # 6. 统计摘要
    ax6 = axes[1, 2]
    ax6.axis('off')
    
    # 创建统计摘要文本
    summary_text = f"""
统计摘要

样本数量: {len(df)}
工具数量: {len(df['Tool'].unique()) if 'Tool' in df.columns else 'N/A'}

平均比对率: {df['Mapping_Rate'].mean():.1f}%
最高比对率: {df['Mapping_Rate'].max():.1f}%
最低比对率: {df['Mapping_Rate'].min():.1f}%

平均配对率: {df['Pairing_Rate'].mean():.1f}%
最高配对率: {df['Pairing_Rate'].max():.1f}%
最低配对率: {df['Pairing_Rate'].min():.1f}%

推荐设置:
最佳工具: {df.loc[df['Mapping_Rate'].idxmax(), 'Tool'] if 'Tool' in df.columns else 'N/A'}
最佳模式: {df.loc[df['Mapping_Rate'].idxmax(), 'Mode'] if 'Mode' in df.columns else 'N/A'}
"""
    
    ax6.text(0.1, 0.9, summary_text, transform=ax6.transAxes, 
            fontsize=11, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig('results/summary_dashboard.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("总结仪表板已生成: results/summary_dashboard.png")

def main():
    """主函数"""
    print("开始生成比对结果可视化...")
    
    # 检查结果目录
    if not os.path.exists("results"):
        print("错误：results目录不存在")
        sys.exit(1)
    
    # 加载数据
    print("加载比对数据...")
    df = load_alignment_data()
    
    if df.empty:
        print("未找到比对数据，无法生成可视化")
        sys.exit(1)
    
    print(f"找到 {len(df)} 个样本的数据")
    
    # 创建各种可视化图表
    print("创建性能比较图...")
    create_performance_comparison(df)
    
    print("创建MAPQ分析图...")
    create_mapq_analysis()
    
    print("创建覆盖度分析图...")
    create_coverage_analysis()
    
    print("创建插入片段分析图...")
    create_insert_size_analysis()
    
    print("创建总结仪表板...")
    create_summary_dashboard()
    
    print("\n可视化完成！生成的图表文件:")
    print("  - results/performance_comparison.png")
    print("  - results/mapq_analysis.png")
    print("  - results/coverage_analysis.png")
    print("  - results/insert_size_analysis.png")
    print("  - results/summary_dashboard.png")
    
    print("\n建议使用图片查看器或在实验报告中查看这些图表。")

if __name__ == "__main__":
    main()