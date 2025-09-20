#!/usr/bin/env python3
"""
绘制Ti/Tv比值分析图
作者：王运生
日期：2025-01-20
用法：python3 plot_titv_ratio.py input.vcf output_dir/
"""

import sys
import argparse
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from collections import defaultdict
import os

# 设置中文字体
matplotlib.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus'] = False

def classify_variant(ref, alt):
    """分类变异类型"""
    # 只处理单核苷酸变异
    if len(ref) != 1 or len(alt) != 1:
        return 'other'
    
    # 转换 (Transition): A<->G, C<->T
    transitions = {
        ('A', 'G'), ('G', 'A'),
        ('C', 'T'), ('T', 'C')
    }
    
    # 颠换 (Transversion): A/G<->C/T
    transversions = {
        ('A', 'C'), ('C', 'A'),
        ('A', 'T'), ('T', 'A'),
        ('G', 'C'), ('C', 'G'),
        ('G', 'T'), ('T', 'G')
    }
    
    variant_pair = (ref.upper(), alt.upper())
    
    if variant_pair in transitions:
        return 'transition'
    elif variant_pair in transversions:
        return 'transversion'
    else:
        return 'other'

def parse_vcf_titv(vcf_file):
    """解析VCF文件计算Ti/Tv信息"""
    
    # 全局统计
    global_stats = {
        'transitions': 0,
        'transversions': 0,
        'other': 0
    }
    
    # 详细变异类型统计
    transition_types = defaultdict(int)
    transversion_types = defaultdict(int)
    
    # 染色体统计
    chromosome_stats = defaultdict(lambda: {'ti': 0, 'tv': 0})
    
    # 质量分层统计
    quality_bins = [(0, 30), (30, 50), (50, 100), (100, float('inf'))]
    quality_stats = {f"{low}-{high if high != float('inf') else '∞'}": {'ti': 0, 'tv': 0} 
                    for low, high in quality_bins}
    
    try:
        with open(vcf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 6:
                    continue
                
                chrom = fields[0]
                ref = fields[3]
                alt = fields[4]
                qual = fields[5]
                
                # 解析质量分数
                try:
                    qual_score = float(qual) if qual != '.' else 0
                except ValueError:
                    qual_score = 0
                
                # 分类变异
                variant_type = classify_variant(ref, alt)
                
                if variant_type == 'transition':
                    global_stats['transitions'] += 1
                    transition_types[f"{ref}>{alt}"] += 1
                    chromosome_stats[chrom]['ti'] += 1
                    
                    # 质量分层
                    for bin_name, (low, high) in zip(quality_stats.keys(), quality_bins):
                        if low <= qual_score < high:
                            quality_stats[bin_name]['ti'] += 1
                            break
                
                elif variant_type == 'transversion':
                    global_stats['transversions'] += 1
                    transversion_types[f"{ref}>{alt}"] += 1
                    chromosome_stats[chrom]['tv'] += 1
                    
                    # 质量分层
                    for bin_name, (low, high) in zip(quality_stats.keys(), quality_bins):
                        if low <= qual_score < high:
                            quality_stats[bin_name]['tv'] += 1
                            break
                
                else:
                    global_stats['other'] += 1
    
    except FileNotFoundError:
        print(f"错误: 找不到文件 {vcf_file}")
        return None
    except Exception as e:
        print(f"错误: 处理文件时出现问题 - {e}")
        return None
    
    return {
        'global_stats': global_stats,
        'transition_types': dict(transition_types),
        'transversion_types': dict(transversion_types),
        'chromosome_stats': dict(chromosome_stats),
        'quality_stats': quality_stats
    }

def plot_titv_analysis(data, output_dir):
    """绘制Ti/Tv分析图"""
    
    global_stats = data['global_stats']
    transition_types = data['transition_types']
    transversion_types = data['transversion_types']
    chromosome_stats = data['chromosome_stats']
    quality_stats = data['quality_stats']
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 计算全局Ti/Tv比值
    if global_stats['transversions'] > 0:
        global_titv = global_stats['transitions'] / global_stats['transversions']
    else:
        global_titv = float('inf')
    
    # 1. Ti/Tv总体统计图
    plt.figure(figsize=(15, 10))
    
    # 1.1 Ti vs Tv柱状图
    plt.subplot(2, 3, 1)
    categories = ['转换 (Ti)', '颠换 (Tv)']
    counts = [global_stats['transitions'], global_stats['transversions']]
    colors = ['skyblue', 'lightcoral']
    
    bars = plt.bar(categories, counts, color=colors, alpha=0.7, edgecolor='black')
    plt.ylabel('变异数量')
    plt.title(f'Ti/Tv统计\n比值: {global_titv:.3f}')
    
    # 添加数值标签
    for bar, count in zip(bars, counts):
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(counts)*0.01,
                f'{count:,}', ha='center', va='bottom')
    
    # 添加质量评估
    if 2.0 <= global_titv <= 2.1:
        quality_text = "优秀"
        quality_color = "green"
    elif 1.8 <= global_titv <= 2.3:
        quality_text = "良好"
        quality_color = "orange"
    else:
        quality_text = "需检查"
        quality_color = "red"
    
    plt.text(0.5, 0.95, f'质量评估: {quality_text}', transform=plt.gca().transAxes,
             ha='center', va='top', bbox=dict(boxstyle='round', facecolor=quality_color, alpha=0.3))
    
    # 1.2 详细变异类型分布
    plt.subplot(2, 3, 2)
    
    # 合并转换和颠换类型
    all_types = {}
    all_types.update(transition_types)
    all_types.update(transversion_types)
    
    # 按数量排序
    sorted_types = sorted(all_types.items(), key=lambda x: x[1], reverse=True)
    
    if len(sorted_types) > 0:
        types = [item[0] for item in sorted_types]
        type_counts = [item[1] for item in sorted_types]
        
        # 区分转换和颠换的颜色
        type_colors = []
        for variant_type in types:
            ref, alt = variant_type.split('>')
            if classify_variant(ref, alt) == 'transition':
                type_colors.append('skyblue')
            else:
                type_colors.append('lightcoral')
        
        plt.bar(range(len(types)), type_counts, color=type_colors, alpha=0.7, edgecolor='black')
        plt.xlabel('变异类型')
        plt.ylabel('数量')
        plt.title('详细变异类型分布')
        plt.xticks(range(len(types)), types, rotation=45)
    
    # 1.3 染色体Ti/Tv比值
    plt.subplot(2, 3, 3)
    
    # 计算各染色体的Ti/Tv比值
    chrom_titv = {}
    for chrom, counts in chromosome_stats.items():
        if counts['tv'] > 0:
            chrom_titv[chrom] = counts['ti'] / counts['tv']
    
    # 按染色体编号排序
    def sort_chromosome(chrom):
        if chrom.isdigit():
            return int(chrom)
        elif chrom == 'X':
            return 23
        elif chrom == 'Y':
            return 24
        elif chrom == 'MT' or chrom == 'M':
            return 25
        else:
            return 26
    
    sorted_chroms = sorted(chrom_titv.items(), key=lambda x: sort_chromosome(x[0]))
    
    if len(sorted_chroms) > 0:
        chroms = [item[0] for item in sorted_chroms[:22]]  # 只显示前22个
        ratios = [item[1] for item in sorted_chroms[:22]]
        
        # 根据比值设置颜色
        colors = ['green' if 2.0 <= r <= 2.1 else 'orange' if 1.8 <= r <= 2.3 else 'red' 
                 for r in ratios]
        
        plt.bar(range(len(chroms)), ratios, color=colors, alpha=0.7, edgecolor='black')
        plt.axhline(y=2.0, color='green', linestyle='--', alpha=0.5, label='正常范围')
        plt.axhline(y=2.1, color='green', linestyle='--', alpha=0.5)
        plt.xlabel('染色体')
        plt.ylabel('Ti/Tv比值')
        plt.title('各染色体Ti/Tv比值')
        plt.xticks(range(len(chroms)), chroms, rotation=45)
        plt.legend()
    
    # 1.4 质量分层Ti/Tv比值
    plt.subplot(2, 3, 4)
    
    quality_bins = list(quality_stats.keys())
    quality_ratios = []
    quality_counts = []
    
    for bin_name in quality_bins:
        ti_count = quality_stats[bin_name]['ti']
        tv_count = quality_stats[bin_name]['tv']
        total_count = ti_count + tv_count
        
        if tv_count > 0:
            ratio = ti_count / tv_count
        else:
            ratio = 0
        
        quality_ratios.append(ratio)
        quality_counts.append(total_count)
    
    # 根据变异数量设置柱子宽度
    max_count = max(quality_counts) if quality_counts else 1
    widths = [0.8 * (count / max_count) for count in quality_counts]
    
    bars = plt.bar(range(len(quality_bins)), quality_ratios, 
                   width=widths, alpha=0.7, edgecolor='black')
    
    # 根据比值设置颜色
    for i, (bar, ratio) in enumerate(zip(bars, quality_ratios)):
        if 2.0 <= ratio <= 2.1:
            bar.set_color('green')
        elif 1.8 <= ratio <= 2.3:
            bar.set_color('orange')
        else:
            bar.set_color('red')
    
    plt.axhline(y=2.0, color='green', linestyle='--', alpha=0.5)
    plt.axhline(y=2.1, color='green', linestyle='--', alpha=0.5)
    plt.xlabel('质量分数区间')
    plt.ylabel('Ti/Tv比值')
    plt.title('不同质量区间Ti/Tv比值')
    plt.xticks(range(len(quality_bins)), quality_bins, rotation=45)
    
    # 添加变异数量标签
    for i, (bar, count) in enumerate(zip(bars, quality_counts)):
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05,
                f'n={count}', ha='center', va='bottom', fontsize=8)
    
    # 1.5 Ti/Tv比值分布直方图
    plt.subplot(2, 3, 5)
    
    # 收集所有染色体的Ti/Tv比值
    all_ratios = [ratio for ratio in chrom_titv.values() if 0 < ratio < 10]  # 过滤异常值
    
    if all_ratios:
        plt.hist(all_ratios, bins=20, alpha=0.7, color='lightblue', edgecolor='black')
        plt.axvline(x=global_titv, color='red', linestyle='-', linewidth=2, label=f'全局比值: {global_titv:.3f}')
        plt.axvline(x=2.0, color='green', linestyle='--', alpha=0.7, label='正常范围')
        plt.axvline(x=2.1, color='green', linestyle='--', alpha=0.7)
        plt.xlabel('Ti/Tv比值')
        plt.ylabel('染色体数量')
        plt.title('Ti/Tv比值分布')
        plt.legend()
    
    # 1.6 累积变异数量
    plt.subplot(2, 3, 6)
    
    # 按染色体大小排序显示累积变异数
    chrom_totals = {chrom: counts['ti'] + counts['tv'] 
                   for chrom, counts in chromosome_stats.items()}
    sorted_chrom_totals = sorted(chrom_totals.items(), key=lambda x: sort_chromosome(x[0]))
    
    if sorted_chrom_totals:
        chroms = [item[0] for item in sorted_chrom_totals[:22]]
        totals = [item[1] for item in sorted_chrom_totals[:22]]
        
        # 分别显示Ti和Tv
        ti_counts = [chromosome_stats[chrom]['ti'] for chrom in chroms]
        tv_counts = [chromosome_stats[chrom]['tv'] for chrom in chroms]
        
        plt.bar(range(len(chroms)), ti_counts, label='转换 (Ti)', 
               color='skyblue', alpha=0.7, edgecolor='black')
        plt.bar(range(len(chroms)), tv_counts, bottom=ti_counts, label='颠换 (Tv)', 
               color='lightcoral', alpha=0.7, edgecolor='black')
        
        plt.xlabel('染色体')
        plt.ylabel('变异数量')
        plt.title('各染色体变异数量分布')
        plt.xticks(range(len(chroms)), chroms, rotation=45)
        plt.legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'titv_analysis.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'titv_analysis.pdf'), bbox_inches='tight')
    print(f"Ti/Tv分析图已保存到: {output_dir}/titv_analysis.png")

def generate_titv_report(data, output_dir):
    """生成Ti/Tv分析报告"""
    
    global_stats = data['global_stats']
    chromosome_stats = data['chromosome_stats']
    quality_stats = data['quality_stats']
    
    # 计算全局Ti/Tv比值
    if global_stats['transversions'] > 0:
        global_titv = global_stats['transitions'] / global_stats['transversions']
    else:
        global_titv = float('inf')
    
    report_file = os.path.join(output_dir, 'titv_report.txt')
    
    with open(report_file, 'w') as f:
        f.write("=== Ti/Tv比值分析报告 ===\n\n")
        
        # 全局统计
        f.write("全局统计:\n")
        f.write(f"  转换 (Transitions): {global_stats['transitions']:,}\n")
        f.write(f"  颠换 (Transversions): {global_stats['transversions']:,}\n")
        f.write(f"  Ti/Tv比值: {global_titv:.3f}\n")
        f.write(f"  其他变异: {global_stats['other']:,}\n\n")
        
        # 质量评估
        f.write("质量评估:\n")
        if 2.0 <= global_titv <= 2.1:
            f.write("  ✓ Ti/Tv比值优秀 (2.0-2.1)\n")
            f.write("  变异检测质量很好\n")
        elif 1.8 <= global_titv <= 2.3:
            f.write("  ! Ti/Tv比值良好 (1.8-2.3)\n")
            f.write("  变异检测质量可接受\n")
        else:
            f.write("  ✗ Ti/Tv比值异常\n")
            if global_titv < 1.8:
                f.write("  比值过低，可能存在大量假阳性颠换变异\n")
            else:
                f.write("  比值过高，可能过度过滤了颠换变异\n")
        f.write("\n")
        
        # 染色体统计
        f.write("染色体Ti/Tv比值 (前10个):\n")
        chrom_titv = {}
        for chrom, counts in chromosome_stats.items():
            if counts['tv'] > 0:
                chrom_titv[chrom] = counts['ti'] / counts['tv']
        
        # 按变异总数排序
        sorted_chroms = sorted(chrom_titv.items(), 
                             key=lambda x: chromosome_stats[x[0]]['ti'] + chromosome_stats[x[0]]['tv'], 
                             reverse=True)[:10]
        
        for chrom, ratio in sorted_chroms:
            ti_count = chromosome_stats[chrom]['ti']
            tv_count = chromosome_stats[chrom]['tv']
            total = ti_count + tv_count
            f.write(f"  {chrom}: {ratio:.3f} (Ti:{ti_count}, Tv:{tv_count}, 总计:{total})\n")
        f.write("\n")
        
        # 质量分层统计
        f.write("不同质量区间Ti/Tv比值:\n")
        for bin_name, counts in quality_stats.items():
            ti_count = counts['ti']
            tv_count = counts['tv']
            total = ti_count + tv_count
            
            if tv_count > 0:
                ratio = ti_count / tv_count
                f.write(f"  QUAL {bin_name}: {ratio:.3f} (Ti:{ti_count}, Tv:{tv_count}, 总计:{total})\n")
            else:
                f.write(f"  QUAL {bin_name}: N/A (Ti:{ti_count}, Tv:{tv_count}, 总计:{total})\n")
        f.write("\n")
        
        # 建议
        f.write("建议:\n")
        if global_titv < 1.8:
            f.write("  - 检查变异过滤参数，可能过于宽松\n")
            f.write("  - 增加质量过滤阈值\n")
            f.write("  - 检查测序数据质量\n")
        elif global_titv > 2.3:
            f.write("  - 检查变异过滤参数，可能过于严格\n")
            f.write("  - 适当放宽过滤标准\n")
            f.write("  - 检查是否存在系统性偏倚\n")
        else:
            f.write("  - Ti/Tv比值正常，变异检测质量良好\n")
            f.write("  - 可以进行后续分析\n")
    
    print(f"Ti/Tv分析报告已保存到: {report_file}")

def main():
    parser = argparse.ArgumentParser(description='绘制Ti/Tv比值分析图')
    parser.add_argument('vcf_file', help='输入的VCF文件')
    parser.add_argument('output_dir', help='输出目录')
    
    args = parser.parse_args()
    
    # 解析VCF文件
    print("解析VCF文件...")
    data = parse_vcf_titv(args.vcf_file)
    
    if data is None:
        sys.exit(1)
    
    # 绘制图表
    print("生成Ti/Tv分析图...")
    plot_titv_analysis(data, args.output_dir)
    
    # 生成报告
    print("生成Ti/Tv分析报告...")
    generate_titv_report(data, args.output_dir)
    
    print("完成!")

if __name__ == "__main__":
    main()