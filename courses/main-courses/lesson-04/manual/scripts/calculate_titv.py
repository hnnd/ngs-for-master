#!/usr/bin/env python3
"""
计算VCF文件中的Ti/Tv比值
作者：王运生
日期：2025-01-20
用法：python3 calculate_titv.py input.vcf
"""

import sys
import argparse
from collections import defaultdict

def parse_vcf_line(line):
    """解析VCF文件行"""
    if line.startswith('#'):
        return None
    
    fields = line.strip().split('\t')
    if len(fields) < 5:
        return None
    
    chrom = fields[0]
    pos = int(fields[1])
    ref = fields[3]
    alt = fields[4]
    
    return {
        'chrom': chrom,
        'pos': pos,
        'ref': ref,
        'alt': alt
    }

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

def calculate_titv_ratio(vcf_file):
    """计算Ti/Tv比值"""
    
    stats = {
        'total_variants': 0,
        'snps': 0,
        'indels': 0,
        'transitions': 0,
        'transversions': 0,
        'other': 0
    }
    
    # 详细统计
    transition_counts = defaultdict(int)
    transversion_counts = defaultdict(int)
    chromosome_stats = defaultdict(lambda: {'ti': 0, 'tv': 0})
    
    try:
        with open(vcf_file, 'r') as f:
            for line in f:
                variant = parse_vcf_line(line)
                if variant is None:
                    continue
                
                stats['total_variants'] += 1
                
                ref = variant['ref']
                alt = variant['alt']
                chrom = variant['chrom']
                
                # 区分SNP和InDel
                if len(ref) == 1 and len(alt) == 1:
                    stats['snps'] += 1
                    
                    # 分类变异类型
                    variant_type = classify_variant(ref, alt)
                    
                    if variant_type == 'transition':
                        stats['transitions'] += 1
                        transition_counts[f"{ref}>{alt}"] += 1
                        chromosome_stats[chrom]['ti'] += 1
                    elif variant_type == 'transversion':
                        stats['transversions'] += 1
                        transversion_counts[f"{ref}>{alt}"] += 1
                        chromosome_stats[chrom]['tv'] += 1
                    else:
                        stats['other'] += 1
                else:
                    stats['indels'] += 1
    
    except FileNotFoundError:
        print(f"错误: 找不到文件 {vcf_file}")
        return None
    except Exception as e:
        print(f"错误: 处理文件时出现问题 - {e}")
        return None
    
    # 计算比值
    if stats['transversions'] > 0:
        titv_ratio = stats['transitions'] / stats['transversions']
    else:
        titv_ratio = float('inf')
    
    return {
        'stats': stats,
        'titv_ratio': titv_ratio,
        'transition_counts': dict(transition_counts),
        'transversion_counts': dict(transversion_counts),
        'chromosome_stats': dict(chromosome_stats)
    }

def print_results(results):
    """打印结果"""
    if results is None:
        return
    
    stats = results['stats']
    titv_ratio = results['titv_ratio']
    
    print("=== Ti/Tv比值分析结果 ===")
    print()
    
    # 基本统计
    print("基本统计:")
    print(f"  总变异数量: {stats['total_variants']:,}")
    print(f"  SNP数量: {stats['snps']:,}")
    print(f"  InDel数量: {stats['indels']:,}")
    print(f"  其他变异: {stats['other']:,}")
    print()
    
    # Ti/Tv统计
    print("Ti/Tv统计:")
    print(f"  转换 (Transitions): {stats['transitions']:,}")
    print(f"  颠换 (Transversions): {stats['transversions']:,}")
    print(f"  Ti/Tv比值: {titv_ratio:.3f}")
    print()
    
    # 质量评估
    print("质量评估:")
    if 2.0 <= titv_ratio <= 2.1:
        quality = "优秀"
        color = "✓"
    elif 1.8 <= titv_ratio <= 2.3:
        quality = "良好"
        color = "!"
    else:
        quality = "需要检查"
        color = "✗"
    
    print(f"  {color} 全基因组Ti/Tv比值: {quality}")
    print(f"    正常范围: 2.0-2.1")
    print(f"    当前值: {titv_ratio:.3f}")
    print()
    
    # 详细转换统计
    if results['transition_counts']:
        print("转换类型统计:")
        for variant_type, count in sorted(results['transition_counts'].items()):
            percentage = (count / stats['transitions']) * 100
            print(f"  {variant_type}: {count:,} ({percentage:.1f}%)")
        print()
    
    if results['transversion_counts']:
        print("颠换类型统计:")
        for variant_type, count in sorted(results['transversion_counts'].items()):
            percentage = (count / stats['transversions']) * 100
            print(f"  {variant_type}: {count:,} ({percentage:.1f}%)")
        print()
    
    # 染色体统计（显示前10个）
    chrom_stats = results['chromosome_stats']
    if chrom_stats:
        print("各染色体Ti/Tv比值 (前10个):")
        sorted_chroms = sorted(chrom_stats.items(), 
                             key=lambda x: x[1]['ti'] + x[1]['tv'], 
                             reverse=True)[:10]
        
        for chrom, counts in sorted_chroms:
            if counts['tv'] > 0:
                chrom_titv = counts['ti'] / counts['tv']
                total_vars = counts['ti'] + counts['tv']
                print(f"  {chrom}: {chrom_titv:.3f} (Ti:{counts['ti']}, Tv:{counts['tv']}, 总计:{total_vars})")
        print()
    
    # 建议
    print("建议:")
    if titv_ratio < 1.8:
        print("  - Ti/Tv比值过低，可能存在大量假阳性变异")
        print("  - 建议检查变异过滤参数")
        print("  - 检查测序数据质量")
    elif titv_ratio > 2.3:
        print("  - Ti/Tv比值过高，可能过度过滤了真实变异")
        print("  - 建议放宽过滤标准")
        print("  - 检查是否存在系统性偏倚")
    else:
        print("  - Ti/Tv比值在正常范围内")
        print("  - 变异检测质量良好")

def main():
    parser = argparse.ArgumentParser(description='计算VCF文件的Ti/Tv比值')
    parser.add_argument('vcf_file', help='输入的VCF文件')
    parser.add_argument('--output', '-o', help='输出结果到文件')
    
    args = parser.parse_args()
    
    # 计算Ti/Tv比值
    results = calculate_titv_ratio(args.vcf_file)
    
    if results is None:
        sys.exit(1)
    
    # 输出结果
    if args.output:
        # 重定向输出到文件
        import contextlib
        with open(args.output, 'w') as f:
            with contextlib.redirect_stdout(f):
                print_results(results)
        print(f"结果已保存到: {args.output}")
    else:
        print_results(results)

if __name__ == "__main__":
    main()