#!/usr/bin/env python3
"""
从VEP注释结果中筛选有害变异
作者：王运生
日期：2025-01-20
用法：python3 filter_damaging_variants.py input.txt > output.txt
"""

import sys
import argparse
import re

def parse_vep_line(line):
    """解析VEP注释行"""
    if line.startswith('#'):
        return None
    
    fields = line.strip().split('\t')
    if len(fields) < 14:  # VEP输出至少有14列
        return None
    
    # VEP标准输出格式
    variant_info = {
        'uploaded_variation': fields[0],
        'location': fields[1],
        'allele': fields[2],
        'gene': fields[3],
        'feature': fields[4],
        'feature_type': fields[5],
        'consequence': fields[6],
        'cdna_position': fields[7],
        'cds_position': fields[8],
        'protein_position': fields[9],
        'amino_acids': fields[10],
        'codons': fields[11],
        'existing_variation': fields[12],
        'extra': fields[13] if len(fields) > 13 else ''
    }
    
    return variant_info

def parse_extra_field(extra_field):
    """解析VEP的Extra字段"""
    extra_info = {}
    
    if not extra_field:
        return extra_info
    
    # Extra字段格式: KEY=VALUE;KEY=VALUE
    pairs = extra_field.split(';')
    for pair in pairs:
        if '=' in pair:
            key, value = pair.split('=', 1)
            extra_info[key] = value
    
    return extra_info

def is_damaging_consequence(consequence):
    """判断变异后果是否有害"""
    # 高影响后果
    high_impact = {
        'transcript_ablation',
        'splice_acceptor_variant',
        'splice_donor_variant',
        'stop_gained',
        'frameshift_variant',
        'stop_lost',
        'start_lost',
        'transcript_amplification'
    }
    
    # 中等影响后果
    moderate_impact = {
        'inframe_insertion',
        'inframe_deletion',
        'missense_variant',
        'protein_altering_variant'
    }
    
    # 检查是否包含有害后果
    consequences = consequence.split(',')
    for cons in consequences:
        cons = cons.strip()
        if cons in high_impact or cons in moderate_impact:
            return True, cons
    
    return False, None

def evaluate_sift_score(sift_info):
    """评估SIFT分数"""
    if not sift_info:
        return None, None
    
    # SIFT格式: prediction(score)
    match = re.match(r'(\w+)\(([\d.]+)\)', sift_info)
    if match:
        prediction = match.group(1)
        score = float(match.group(2))
        
        # SIFT: ≤0.05 有害, >0.05 耐受
        is_damaging = score <= 0.05
        return is_damaging, score
    
    return None, None

def evaluate_polyphen_score(polyphen_info):
    """评估PolyPhen分数"""
    if not polyphen_info:
        return None, None
    
    # PolyPhen格式: prediction(score)
    match = re.match(r'(\w+)\(([\d.]+)\)', polyphen_info)
    if match:
        prediction = match.group(1)
        score = float(match.group(2))
        
        # PolyPhen: >0.5 可能有害, >0.85 很可能有害
        is_damaging = score > 0.5
        return is_damaging, score
    
    return None, None

def evaluate_cadd_score(cadd_info):
    """评估CADD分数"""
    if not cadd_info:
        return None, None
    
    try:
        score = float(cadd_info)
        # CADD: ≥15 有害 (前1%有害变异)
        is_damaging = score >= 15.0
        return is_damaging, score
    except ValueError:
        return None, None

def get_allele_frequency(extra_info):
    """获取等位基因频率"""
    frequencies = {}
    
    # 常见的频率字段
    freq_fields = {
        'gnomAD_AF': 'gnomAD',
        'AF': '1000G',
        'ExAC_AF': 'ExAC',
        'ESP_AF': 'ESP'
    }
    
    for field, source in freq_fields.items():
        if field in extra_info:
            try:
                freq = float(extra_info[field])
                frequencies[source] = freq
            except ValueError:
                continue
    
    return frequencies

def is_rare_variant(frequencies, threshold=0.01):
    """判断是否为罕见变异"""
    if not frequencies:
        return True  # 没有频率信息，假设为罕见
    
    # 如果任何数据库中频率超过阈值，则不是罕见变异
    for source, freq in frequencies.items():
        if freq > threshold:
            return False
    
    return True

def filter_damaging_variants(vep_file, output_file=None):
    """筛选有害变异"""
    
    damaging_variants = []
    total_variants = 0
    
    try:
        with open(vep_file, 'r') as f:
            # 跳过头部注释
            for line in f:
                if line.startswith('## ENSEMBL VARIANT EFFECT PREDICTOR'):
                    continue
                if line.startswith('#Uploaded_variation'):
                    # 这是列标题行
                    header = line.strip()
                    if output_file:
                        damaging_variants.append(header)
                    continue
                if line.startswith('#'):
                    continue
                
                variant = parse_vep_line(line)
                if variant is None:
                    continue
                
                total_variants += 1
                
                # 解析Extra字段
                extra_info = parse_extra_field(variant['extra'])
                
                # 评估变异是否有害
                is_damaging = False
                damaging_reasons = []
                
                # 1. 检查变异后果
                consequence_damaging, consequence_type = is_damaging_consequence(variant['consequence'])
                if consequence_damaging:
                    is_damaging = True
                    damaging_reasons.append(f"consequence:{consequence_type}")
                
                # 2. 检查SIFT分数
                if 'SIFT' in extra_info:
                    sift_damaging, sift_score = evaluate_sift_score(extra_info['SIFT'])
                    if sift_damaging:
                        is_damaging = True
                        damaging_reasons.append(f"SIFT:{sift_score}")
                
                # 3. 检查PolyPhen分数
                if 'PolyPhen' in extra_info:
                    polyphen_damaging, polyphen_score = evaluate_polyphen_score(extra_info['PolyPhen'])
                    if polyphen_damaging:
                        is_damaging = True
                        damaging_reasons.append(f"PolyPhen:{polyphen_score}")
                
                # 4. 检查CADD分数
                if 'CADD_PHRED' in extra_info:
                    cadd_damaging, cadd_score = evaluate_cadd_score(extra_info['CADD_PHRED'])
                    if cadd_damaging:
                        is_damaging = True
                        damaging_reasons.append(f"CADD:{cadd_score}")
                
                # 5. 检查是否为罕见变异
                frequencies = get_allele_frequency(extra_info)
                is_rare = is_rare_variant(frequencies)
                
                # 只保留有害且罕见的变异
                if is_damaging and is_rare:
                    # 添加筛选原因到Extra字段
                    reasons_str = "|".join(damaging_reasons)
                    freq_str = "|".join([f"{k}:{v}" for k, v in frequencies.items()])
                    
                    enhanced_extra = variant['extra']
                    if enhanced_extra:
                        enhanced_extra += f";DAMAGING_REASONS={reasons_str}"
                        if freq_str:
                            enhanced_extra += f";FREQUENCIES={freq_str}"
                    else:
                        enhanced_extra = f"DAMAGING_REASONS={reasons_str}"
                        if freq_str:
                            enhanced_extra += f";FREQUENCIES={freq_str}"
                    
                    # 重构输出行
                    output_fields = [
                        variant['uploaded_variation'],
                        variant['location'],
                        variant['allele'],
                        variant['gene'],
                        variant['feature'],
                        variant['feature_type'],
                        variant['consequence'],
                        variant['cdna_position'],
                        variant['cds_position'],
                        variant['protein_position'],
                        variant['amino_acids'],
                        variant['codons'],
                        variant['existing_variation'],
                        enhanced_extra
                    ]
                    
                    damaging_variants.append('\t'.join(output_fields))
    
    except FileNotFoundError:
        print(f"错误: 找不到文件 {vep_file}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"错误: 处理文件时出现问题 - {e}", file=sys.stderr)
        return None
    
    # 输出结果
    if output_file:
        with open(output_file, 'w') as f:
            for variant in damaging_variants:
                f.write(variant + '\n')
    else:
        for variant in damaging_variants:
            print(variant)
    
    # 输出统计信息到stderr
    print(f"总变异数量: {total_variants}", file=sys.stderr)
    print(f"有害变异数量: {len(damaging_variants)-1}", file=sys.stderr)  # -1 因为包含了header
    if total_variants > 0:
        percentage = ((len(damaging_variants)-1) / total_variants) * 100
        print(f"有害变异比例: {percentage:.2f}%", file=sys.stderr)
    
    return len(damaging_variants) - 1  # 返回有害变异数量（不包含header）

def main():
    parser = argparse.ArgumentParser(description='筛选VEP注释结果中的有害变异')
    parser.add_argument('vep_file', help='VEP注释结果文件')
    parser.add_argument('--output', '-o', help='输出文件（默认输出到stdout）')
    parser.add_argument('--freq-threshold', type=float, default=0.01, 
                       help='等位基因频率阈值（默认0.01）')
    
    args = parser.parse_args()
    
    # 筛选有害变异
    result = filter_damaging_variants(args.vep_file, args.output)
    
    if result is None:
        sys.exit(1)

if __name__ == "__main__":
    main()