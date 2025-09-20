#!/usr/bin/env python3
"""
绘制变异密度分布图
作者：王运生
日期：2025-01-20
用法：python3 plot_variant_density.py input.vcf output_dir/
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

def parse_vcf_positions(vcf_file):
    """解析VCF文件中的变异位置信息"""
    
    variants = []
    chromosome_variants = defaultdict(list)
    variant_types = defaultdict(list)
    
    try:
        with open(vcf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 5:
                    continue
                
                chrom = fields[0]
                pos = int(fields[1])
                ref = fields[3]
                alt = fields[4]
                
                # 分类变异类型
                if len(ref) == 1 and len(alt) == 1:
                    variant_type = 'SNP'
                elif len(ref) > len(alt):
                    variant_type = 'Deletion'
                elif len(ref) < len(alt):
                    variant_type = 'Insertion'
                else:
                    variant_type = 'Complex'
                
                variant_info = {
                    'chrom': chrom,
                    'pos': pos,
                    'ref': ref,
                    'alt': alt,
                    'type': variant_type
                }
                
                variants.append(variant_info)
                chromosome_variants[chrom].append(pos)
                variant_types[variant_type].append((chrom, pos))
    
    except FileNotFoundError:
        print(f"错误: 找不到文件 {vcf_file}")
        return None
    except Exception as e:
        print(f"错误: 处理文件时出现问题 - {e}")
        return None
    
    return {
        'variants': variants,
        'chromosome_variants': dict(chromosome_variants),
        'variant_types': dict(variant_types)
    }

def calculate_variant_density(positions, window_size=1000000):
    """计算变异密度"""
    if not positions:
        return [], []
    
    min_pos = min(positions)
    max_pos = max(positions)
    
    # 创建窗口
    windows = []
    densities = []
    
    current_pos = min_pos
    while current_pos <= max_pos:
        window_end = current_pos + window_size
        
        # 计算窗口内的变异数量
        count = sum(1 for pos in positions if current_pos <= pos < window_end)
        density = count / (window_size / 1000000)  # 每Mb的变异数
        
        windows.append(current_pos + window_size // 2)  # 窗口中心
        densities.append(density)
        
        current_pos += window_size
    
    return windows, densities

def plot_variant_density(data, output_dir, window_size=1000000):
    """绘制变异密度图"""
    
    chromosome_variants = data['chromosome_variants']
    variant_types = data['variant_types']
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. 全基因组变异密度图
    plt.figure(figsize=(16, 10))
    
    # 1.1 各染色体变异密度
    plt.subplot(3, 2, 1)
    
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
    
    sorted_chroms = sorted(chromosome_variants.items(), key=lambda x: sort_chromosome(x[0]))
    
    # 计算每个染色体的变异密度
    chrom_names = []
    chrom_densities = []
    chrom_counts = []
    
    for chrom, positions in sorted_chroms[:24]:  # 只显示前24个染色体
        if positions:
            # 估算染色体长度（使用最大位置）
            chrom_length = max(positions)
            density = len(positions) / (chrom_length / 1000000)  # 每Mb变异数
            
            chrom_names.append(chrom)
            chrom_densities.append(density)
            chrom_counts.append(len(positions))
    
    bars = plt.bar(range(len(chrom_names)), chrom_densities, 
                   alpha=0.7, color='skyblue', edgecolor='black')
    plt.xlabel('染色体')
    plt.ylabel('变异密度 (变异数/Mb)')
    plt.title('各染色体变异密度')
    plt.xticks(range(len(chrom_names)), chrom_names, rotation=45)
    
    # 添加变异数量标签
    for i, (bar, count) in enumerate(zip(bars, chrom_counts)):
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(chrom_densities)*0.01,
                f'{count}', ha='center', va='bottom', fontsize=8)
    
    # 1.2 变异类型密度比较
    plt.subplot(3, 2, 2)
    
    type_densities = {}
    type_counts = {}
    
    for variant_type, positions in variant_types.items():
        if positions:
            # 计算总的基因组覆盖范围
            all_positions = [pos for chrom, pos in positions]
            if all_positions:
                # 使用人类基因组大小估算 (~3Gb)
                genome_size = 3000  # Mb
                density = len(all_positions) / genome_size
                type_densities[variant_type] = density
                type_counts[variant_type] = len(all_positions)
    
    if type_densities:
        types = list(type_densities.keys())
        densities = list(type_densities.values())
        counts = [type_counts[t] for t in types]
        
        colors = ['lightblue', 'lightgreen', 'lightcoral', 'lightyellow'][:len(types)]
        bars = plt.bar(types, densities, color=colors, alpha=0.7, edgecolor='black')
        plt.ylabel('变异密度 (变异数/Mb)')
        plt.title('不同类型变异密度')
        plt.xticks(rotation=45)
        
        # 添加数量标签
        for bar, count in zip(bars, counts):
            plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(densities)*0.01,
                    f'{count:,}', ha='center', va='bottom', fontsize=9)
    
    # 1.3 染色体1的详细密度分布（示例）
    plt.subplot(3, 2, 3)
    
    if '1' in chromosome_variants and chromosome_variants['1']:
        positions = chromosome_variants['1']
        windows, densities = calculate_variant_density(positions, window_size)
        
        plt.plot(np.array(windows) / 1000000, densities, 'b-', linewidth=1, alpha=0.7)
        plt.fill_between(np.array(windows) / 1000000, densities, alpha=0.3)
        plt.xlabel('位置 (Mb)')
        plt.ylabel('变异密度 (变异数/Mb)')
        plt.title(f'染色体1变异密度分布 (窗口: {window_size//1000}kb)')
        plt.grid(True, alpha=0.3)
    
    # 1.4 变异间距分布
    plt.subplot(3, 2, 4)
    
    all_intervals = []
    for chrom, positions in chromosome_variants.items():
        if len(positions) > 1:
            sorted_pos = sorted(positions)
            intervals = [sorted_pos[i+1] - sorted_pos[i] for i in range(len(sorted_pos)-1)]
            all_intervals.extend(intervals)
    
    if all_intervals:
        # 过滤极端值
        filtered_intervals = [interval for interval in all_intervals if interval < 100000]
        
        plt.hist(filtered_intervals, bins=50, alpha=0.7, color='lightgreen', edgecolor='black')
        plt.xlabel('变异间距 (bp)')
        plt.ylabel('频次')
        plt.title('变异间距分布')
        plt.yscale('log')
        
        # 添加统计信息
        mean_interval = np.mean(all_intervals)
        median_interval = np.median(all_intervals)
        plt.axvline(x=mean_interval, color='red', linestyle='--', label=f'平均: {mean_interval:.0f}bp')
        plt.axvline(x=median_interval, color='green', linestyle='--', label=f'中位数: {median_interval:.0f}bp')
        plt.legend()
    
    # 1.5 变异密度热图（前10个染色体）
    plt.subplot(3, 2, 5)
    
    # 创建密度矩阵
    density_matrix = []
    chrom_labels = []
    
    for chrom, positions in sorted_chroms[:10]:
        if positions:
            windows, densities = calculate_variant_density(positions, window_size * 5)  # 使用更大窗口
            # 标准化到固定长度
            if len(densities) > 50:
                densities = densities[:50]
            elif len(densities) < 50:
                densities.extend([0] * (50 - len(densities)))
            
            density_matrix.append(densities)
            chrom_labels.append(chrom)
    
    if density_matrix:
        plt.imshow(density_matrix, cmap='YlOrRd', aspect='auto', interpolation='nearest')
        plt.colorbar(label='变异密度')
        plt.xlabel('基因组位置 (窗口)')
        plt.ylabel('染色体')
        plt.title('变异密度热图')
        plt.yticks(range(len(chrom_labels)), chrom_labels)
    
    # 1.6 累积变异分布
    plt.subplot(3, 2, 6)
    
    # 选择几个主要染色体显示累积分布
    major_chroms = ['1', '2', '3', 'X']
    colors = ['blue', 'red', 'green', 'orange']
    
    for i, chrom in enumerate(major_chroms):
        if chrom in chromosome_variants and chromosome_variants[chrom]:
            positions = sorted(chromosome_variants[chrom])
            cumulative = np.arange(1, len(positions) + 1)
            
            plt.plot(np.array(positions) / 1000000, cumulative, 
                    color=colors[i], label=f'Chr{chrom}', linewidth=2)
    
    plt.xlabel('染色体位置 (Mb)')
    plt.ylabel('累积变异数')
    plt.title('主要染色体累积变异分布')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'variant_density.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'variant_density.pdf'), bbox_inches='tight')
    print(f"变异密度图已保存到: {output_dir}/variant_density.png")
    
    # 2. 单独的全基因组概览图
    plt.figure(figsize=(20, 6))
    
    # 创建全基因组变异分布图
    x_offset = 0
    chrom_boundaries = []
    all_x_positions = []
    all_y_positions = []
    
    # 染色体长度估算（Mb）
    chrom_lengths = {
        '1': 249, '2': 243, '3': 198, '4': 191, '5': 181, '6': 171,
        '7': 159, '8': 146, '9': 138, '10': 134, '11': 135, '12': 133,
        '13': 115, '14': 107, '15': 102, '16': 90, '17': 83, '18': 80,
        '19': 59, '20': 63, '21': 48, '22': 51, 'X': 155, 'Y': 59
    }
    
    for chrom, positions in sorted_chroms:
        if chrom in chrom_lengths and positions:
            # 转换位置到全基因组坐标
            chrom_x = [x_offset + pos/1000000 for pos in positions]
            chrom_y = [1] * len(positions)  # 所有变异在同一水平线上
            
            all_x_positions.extend(chrom_x)
            all_y_positions.extend(chrom_y)
            
            # 记录染色体边界
            chrom_boundaries.append((x_offset, x_offset + chrom_lengths[chrom], chrom))
            x_offset += chrom_lengths[chrom] + 10  # 染色体间隔
    
    # 绘制变异点
    plt.scatter(all_x_positions, all_y_positions, s=0.1, alpha=0.6, c='blue')
    
    # 添加染色体分隔线和标签
    for start, end, chrom in chrom_boundaries:
        plt.axvline(x=start, color='red', linestyle='-', alpha=0.3)
        plt.axvline(x=end, color='red', linestyle='-', alpha=0.3)
        plt.text((start + end) / 2, 1.1, chrom, ha='center', va='bottom', fontsize=10)
    
    plt.xlabel('基因组位置 (Mb)')
    plt.ylabel('')
    plt.title('全基因组变异分布概览')
    plt.ylim(0.5, 1.5)
    plt.yticks([])
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'genome_overview.png'), dpi=300, bbox_inches='tight')
    print(f"全基因组概览图已保存到: {output_dir}/genome_overview.png")

def generate_density_report(data, output_dir):
    """生成变异密度分析报告"""
    
    chromosome_variants = data['chromosome_variants']
    variant_types = data['variant_types']
    
    report_file = os.path.join(output_dir, 'density_report.txt')
    
    with open(report_file, 'w') as f:
        f.write("=== 变异密度分析报告 ===\n\n")
        
        # 总体统计
        total_variants = sum(len(positions) for positions in chromosome_variants.values())
        f.write(f"总变异数量: {total_variants:,}\n")
        f.write(f"涉及染色体数量: {len(chromosome_variants)}\n\n")
        
        # 变异类型统计
        f.write("变异类型分布:\n")
        for variant_type, positions in variant_types.items():
            count = len(positions)
            percentage = (count / total_variants) * 100 if total_variants > 0 else 0
            f.write(f"  {variant_type}: {count:,} ({percentage:.1f}%)\n")
        f.write("\n")
        
        # 染色体密度统计
        f.write("各染色体变异密度 (前15个):\n")
        
        # 染色体长度估算（Mb）
        chrom_lengths = {
            '1': 249, '2': 243, '3': 198, '4': 191, '5': 181, '6': 171,
            '7': 159, '8': 146, '9': 138, '10': 134, '11': 135, '12': 133,
            '13': 115, '14': 107, '15': 102, '16': 90, '17': 83, '18': 80,
            '19': 59, '20': 63, '21': 48, '22': 51, 'X': 155, 'Y': 59
        }
        
        chrom_densities = []
        for chrom, positions in chromosome_variants.items():
            if chrom in chrom_lengths and positions:
                density = len(positions) / chrom_lengths[chrom]
                chrom_densities.append((chrom, len(positions), density))
        
        # 按密度排序
        chrom_densities.sort(key=lambda x: x[2], reverse=True)
        
        for chrom, count, density in chrom_densities[:15]:
            f.write(f"  Chr{chrom}: {density:.1f} 变异/Mb ({count:,} 变异)\n")
        f.write("\n")
        
        # 变异间距统计
        f.write("变异间距统计:\n")
        all_intervals = []
        for chrom, positions in chromosome_variants.items():
            if len(positions) > 1:
                sorted_pos = sorted(positions)
                intervals = [sorted_pos[i+1] - sorted_pos[i] for i in range(len(sorted_pos)-1)]
                all_intervals.extend(intervals)
        
        if all_intervals:
            f.write(f"  平均间距: {np.mean(all_intervals):.0f} bp\n")
            f.write(f"  中位数间距: {np.median(all_intervals):.0f} bp\n")
            f.write(f"  最小间距: {np.min(all_intervals)} bp\n")
            f.write(f"  最大间距: {np.max(all_intervals):,} bp\n")
        f.write("\n")
        
        # 密度评估
        f.write("密度评估:\n")
        
        # 计算全基因组平均密度
        total_genome_size = sum(chrom_lengths.values())  # 约3000 Mb
        global_density = total_variants / total_genome_size
        
        f.write(f"  全基因组平均密度: {global_density:.1f} 变异/Mb\n")
        
        if global_density > 20:
            f.write("  ✓ 变异密度正常，符合人类基因组预期\n")
        elif global_density > 10:
            f.write("  ! 变异密度偏低，可能过度过滤\n")
        else:
            f.write("  ✗ 变异密度过低，需要检查分析流程\n")
        
        # 染色体密度均匀性
        if chrom_densities:
            densities = [x[2] for x in chrom_densities]
            cv = np.std(densities) / np.mean(densities)  # 变异系数
            
            f.write(f"  染色体间密度变异系数: {cv:.3f}\n")
            if cv < 0.3:
                f.write("  ✓ 染色体间变异密度分布均匀\n")
            elif cv < 0.5:
                f.write("  ! 染色体间变异密度存在一定差异\n")
            else:
                f.write("  ✗ 染色体间变异密度差异较大\n")
        
        f.write("\n")
        
        # 建议
        f.write("建议:\n")
        if global_density < 10:
            f.write("  - 检查变异过滤参数是否过于严格\n")
            f.write("  - 验证输入数据的完整性\n")
        elif global_density > 50:
            f.write("  - 检查是否存在假阳性变异\n")
            f.write("  - 考虑增加质量过滤\n")
        else:
            f.write("  - 变异密度分布正常\n")
            f.write("  - 可以进行后续分析\n")
    
    print(f"变异密度报告已保存到: {report_file}")

def main():
    parser = argparse.ArgumentParser(description='绘制变异密度分布图')
    parser.add_argument('vcf_file', help='输入的VCF文件')
    parser.add_argument('output_dir', help='输出目录')
    parser.add_argument('--window-size', type=int, default=1000000, 
                       help='密度计算窗口大小 (默认1Mb)')
    
    args = parser.parse_args()
    
    # 解析VCF文件
    print("解析VCF文件...")
    data = parse_vcf_positions(args.vcf_file)
    
    if data is None:
        sys.exit(1)
    
    # 绘制图表
    print("生成变异密度图...")
    plot_variant_density(data, args.output_dir, args.window_size)
    
    # 生成报告
    print("生成变异密度报告...")
    generate_density_report(data, args.output_dir)
    
    print("完成!")

if __name__ == "__main__":
    main()