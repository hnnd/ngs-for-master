#!/usr/bin/env python3
"""
绘制VCF文件中变异质量分布图
作者：王运生
日期：2025-01-20
用法：python3 plot_quality_distribution.py input.vcf output_dir/
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

def parse_vcf_quality(vcf_file):
    """解析VCF文件中的质量信息"""
    
    qualities = []
    depths = []
    variant_types = {'SNP': [], 'InDel': []}
    chromosome_counts = defaultdict(int)
    
    try:
        with open(vcf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 8:
                    continue
                
                chrom = fields[0]
                pos = int(fields[1])
                ref = fields[3]
                alt = fields[4]
                qual = fields[5]
                info = fields[7]
                
                # 解析质量分数
                try:
                    qual_score = float(qual) if qual != '.' else 0
                    qualities.append(qual_score)
                except ValueError:
                    continue
                
                # 解析深度信息
                depth = None
                for info_item in info.split(';'):
                    if info_item.startswith('DP='):
                        try:
                            depth = int(info_item.split('=')[1])
                            depths.append(depth)
                        except ValueError:
                            pass
                        break
                
                # 分类变异类型
                if len(ref) == 1 and len(alt) == 1:
                    variant_types['SNP'].append(qual_score)
                else:
                    variant_types['InDel'].append(qual_score)
                
                # 统计染色体
                chromosome_counts[chrom] += 1
    
    except FileNotFoundError:
        print(f"错误: 找不到文件 {vcf_file}")
        return None
    except Exception as e:
        print(f"错误: 处理文件时出现问题 - {e}")
        return None
    
    return {
        'qualities': qualities,
        'depths': depths,
        'variant_types': variant_types,
        'chromosome_counts': chromosome_counts
    }

def plot_quality_distribution(data, output_dir):
    """绘制质量分布图"""
    
    qualities = data['qualities']
    depths = data['depths']
    variant_types = data['variant_types']
    chromosome_counts = data['chromosome_counts']
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. 质量分数分布直方图
    plt.figure(figsize=(12, 8))
    
    plt.subplot(2, 2, 1)
    plt.hist(qualities, bins=50, alpha=0.7, color='skyblue', edgecolor='black')
    plt.axvline(x=30, color='red', linestyle='--', label='QUAL=30阈值')
    plt.xlabel('变异质量分数 (QUAL)')
    plt.ylabel('变异数量')
    plt.title('变异质量分数分布')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 添加统计信息
    mean_qual = np.mean(qualities)
    median_qual = np.median(qualities)
    high_qual_count = sum(1 for q in qualities if q >= 30)
    high_qual_percent = (high_qual_count / len(qualities)) * 100
    
    plt.text(0.02, 0.98, f'平均值: {mean_qual:.1f}\n中位数: {median_qual:.1f}\nQUAL≥30: {high_qual_percent:.1f}%', 
             transform=plt.gca().transAxes, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # 2. 覆盖深度分布
    if depths:
        plt.subplot(2, 2, 2)
        # 限制深度范围以便更好显示
        filtered_depths = [d for d in depths if d <= 200]
        plt.hist(filtered_depths, bins=50, alpha=0.7, color='lightgreen', edgecolor='black')
        plt.axvline(x=20, color='red', linestyle='--', label='DP=20阈值')
        plt.xlabel('覆盖深度 (DP)')
        plt.ylabel('变异数量')
        plt.title('覆盖深度分布')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # 添加统计信息
        mean_depth = np.mean(depths)
        median_depth = np.median(depths)
        high_depth_count = sum(1 for d in depths if d >= 20)
        high_depth_percent = (high_depth_count / len(depths)) * 100
        
        plt.text(0.02, 0.98, f'平均值: {mean_depth:.1f}\n中位数: {median_depth:.1f}\nDP≥20: {high_depth_percent:.1f}%', 
                 transform=plt.gca().transAxes, verticalalignment='top',
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # 3. SNP vs InDel质量比较
    plt.subplot(2, 2, 3)
    snp_quals = variant_types['SNP']
    indel_quals = variant_types['InDel']
    
    plt.hist([snp_quals, indel_quals], bins=30, alpha=0.7, 
             label=['SNP', 'InDel'], color=['blue', 'orange'])
    plt.xlabel('变异质量分数')
    plt.ylabel('变异数量')
    plt.title('SNP vs InDel质量分布')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 4. 染色体变异数量分布
    plt.subplot(2, 2, 4)
    
    # 只显示前20个染色体
    sorted_chroms = sorted(chromosome_counts.items(), 
                          key=lambda x: int(x[0]) if x[0].isdigit() else 999)[:20]
    
    chroms = [item[0] for item in sorted_chroms]
    counts = [item[1] for item in sorted_chroms]
    
    plt.bar(chroms, counts, alpha=0.7, color='coral')
    plt.xlabel('染色体')
    plt.ylabel('变异数量')
    plt.title('各染色体变异数量分布')
    plt.xticks(rotation=45)
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'quality_distribution.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'quality_distribution.pdf'), bbox_inches='tight')
    print(f"质量分布图已保存到: {output_dir}/quality_distribution.png")
    
    # 5. 质量分数累积分布图
    plt.figure(figsize=(10, 6))
    
    sorted_quals = sorted(qualities)
    cumulative_percent = np.arange(1, len(sorted_quals) + 1) / len(sorted_quals) * 100
    
    plt.plot(sorted_quals, cumulative_percent, linewidth=2, color='blue')
    plt.axvline(x=30, color='red', linestyle='--', label='QUAL=30阈值')
    plt.axhline(y=90, color='green', linestyle='--', label='90%变异')
    
    plt.xlabel('变异质量分数 (QUAL)')
    plt.ylabel('累积百分比 (%)')
    plt.title('变异质量分数累积分布')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 添加关键点标注
    qual_90 = sorted_quals[int(len(sorted_quals) * 0.9)]
    plt.annotate(f'90%变异质量: {qual_90:.1f}', 
                xy=(qual_90, 90), xytext=(qual_90 + 20, 80),
                arrowprops=dict(arrowstyle='->', color='green'))
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'quality_cumulative.png'), dpi=300, bbox_inches='tight')
    print(f"累积分布图已保存到: {output_dir}/quality_cumulative.png")
    
    # 6. 质量vs深度散点图
    if depths and len(depths) == len(qualities):
        plt.figure(figsize=(10, 6))
        
        # 随机采样以避免过度绘制
        if len(qualities) > 10000:
            indices = np.random.choice(len(qualities), 10000, replace=False)
            sample_quals = [qualities[i] for i in indices]
            sample_depths = [depths[i] for i in indices]
        else:
            sample_quals = qualities
            sample_depths = depths
        
        plt.scatter(sample_depths, sample_quals, alpha=0.5, s=1)
        plt.axhline(y=30, color='red', linestyle='--', label='QUAL=30')
        plt.axvline(x=20, color='green', linestyle='--', label='DP=20')
        
        plt.xlabel('覆盖深度 (DP)')
        plt.ylabel('变异质量分数 (QUAL)')
        plt.title('变异质量 vs 覆盖深度')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # 限制显示范围
        plt.xlim(0, 200)
        plt.ylim(0, 200)
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'quality_vs_depth.png'), dpi=300, bbox_inches='tight')
        print(f"质量vs深度图已保存到: {output_dir}/quality_vs_depth.png")

def generate_summary_report(data, output_dir):
    """生成质量控制总结报告"""
    
    qualities = data['qualities']
    depths = data['depths']
    variant_types = data['variant_types']
    
    report_file = os.path.join(output_dir, 'quality_summary.txt')
    
    with open(report_file, 'w') as f:
        f.write("=== 变异质量控制总结报告 ===\n\n")
        
        # 基本统计
        f.write("基本统计:\n")
        f.write(f"  总变异数量: {len(qualities):,}\n")
        f.write(f"  SNP数量: {len(variant_types['SNP']):,}\n")
        f.write(f"  InDel数量: {len(variant_types['InDel']):,}\n")
        f.write(f"  SNP/InDel比值: {len(variant_types['SNP'])/len(variant_types['InDel']):.2f}\n\n")
        
        # 质量统计
        f.write("质量分数统计:\n")
        f.write(f"  平均质量: {np.mean(qualities):.2f}\n")
        f.write(f"  中位数质量: {np.median(qualities):.2f}\n")
        f.write(f"  最小质量: {np.min(qualities):.2f}\n")
        f.write(f"  最大质量: {np.max(qualities):.2f}\n")
        
        high_qual_count = sum(1 for q in qualities if q >= 30)
        high_qual_percent = (high_qual_count / len(qualities)) * 100
        f.write(f"  QUAL≥30变异: {high_qual_count:,} ({high_qual_percent:.1f}%)\n\n")
        
        # 深度统计
        if depths:
            f.write("覆盖深度统计:\n")
            f.write(f"  平均深度: {np.mean(depths):.2f}\n")
            f.write(f"  中位数深度: {np.median(depths):.2f}\n")
            f.write(f"  最小深度: {np.min(depths)}\n")
            f.write(f"  最大深度: {np.max(depths)}\n")
            
            high_depth_count = sum(1 for d in depths if d >= 20)
            high_depth_percent = (high_depth_count / len(depths)) * 100
            f.write(f"  DP≥20变异: {high_depth_count:,} ({high_depth_percent:.1f}%)\n\n")
        
        # 质量评估
        f.write("质量评估:\n")
        
        if high_qual_percent >= 90:
            f.write("  ✓ 高质量变异比例优秀 (≥90%)\n")
        elif high_qual_percent >= 80:
            f.write("  ! 高质量变异比例良好 (80-90%)\n")
        else:
            f.write("  ✗ 高质量变异比例需要改进 (<80%)\n")
        
        if depths:
            if high_depth_percent >= 90:
                f.write("  ✓ 高深度变异比例优秀 (≥90%)\n")
            elif high_depth_percent >= 80:
                f.write("  ! 高深度变异比例良好 (80-90%)\n")
            else:
                f.write("  ✗ 高深度变异比例需要改进 (<80%)\n")
        
        # 建议
        f.write("\n建议:\n")
        if high_qual_percent < 80:
            f.write("  - 考虑调整变异过滤参数\n")
            f.write("  - 检查原始数据质量\n")
        if depths and high_depth_percent < 80:
            f.write("  - 检查测序深度是否足够\n")
            f.write("  - 考虑增加测序深度\n")
        if high_qual_percent >= 90 and (not depths or high_depth_percent >= 90):
            f.write("  - 数据质量优秀，可以进行后续分析\n")
    
    print(f"质量控制报告已保存到: {report_file}")

def main():
    parser = argparse.ArgumentParser(description='绘制VCF文件质量分布图')
    parser.add_argument('vcf_file', help='输入的VCF文件')
    parser.add_argument('output_dir', help='输出目录')
    
    args = parser.parse_args()
    
    # 解析VCF文件
    print("解析VCF文件...")
    data = parse_vcf_quality(args.vcf_file)
    
    if data is None:
        sys.exit(1)
    
    # 绘制图表
    print("生成质量分布图...")
    plot_quality_distribution(data, args.output_dir)
    
    # 生成报告
    print("生成质量控制报告...")
    generate_summary_report(data, args.output_dir)
    
    print("完成!")

if __name__ == "__main__":
    main()