#!/usr/bin/env python3
"""
Python数据结构练习脚本
课程：高通量测序数据分析 - Python编程基础
作者：王运生教授
日期：2025年
"""

def list_operations():
    """列表操作练习"""
    print("=== 列表操作练习 ===")
    
    # 创建基因列表
    genes = ["BRCA1", "BRCA2", "TP53", "EGFR", "MYC"]
    print(f"原始基因列表: {genes}")
    
    # 基本访问操作
    print(f"\n--- 基本访问操作 ---")
    print(f"第一个基因: {genes[0]}")
    print(f"最后一个基因: {genes[-1]}")
    print(f"前三个基因: {genes[:3]}")
    print(f"后两个基因: {genes[-2:]}")
    print(f"中间的基因: {genes[1:4]}")
    
    # 修改列表
    print(f"\n--- 修改列表 ---")
    genes_copy = genes.copy()  # 创建副本以保持原列表不变
    
    genes_copy.append("KRAS")
    print(f"添加KRAS后: {genes_copy}")
    
    genes_copy.insert(2, "PIK3CA")
    print(f"在位置2插入PIK3CA: {genes_copy}")
    
    genes_copy.remove("TP53")
    print(f"删除TP53后: {genes_copy}")
    
    popped_gene = genes_copy.pop()
    print(f"弹出最后一个基因 {popped_gene}: {genes_copy}")
    
    # 列表方法
    print(f"\n--- 列表方法 ---")
    print(f"列表长度: {len(genes)}")
    print(f"BRCA1的索引: {genes.index('BRCA1')}")
    print(f"基因数量: {genes.count('BRCA1')}")
    
    # 排序
    sorted_genes = sorted(genes)
    print(f"按字母排序: {sorted_genes}")
    
    # 反转
    reversed_genes = genes[::-1]
    print(f"反转列表: {reversed_genes}")

def list_comprehensions():
    """列表推导式练习"""
    print("\n=== 列表推导式练习 ===")
    
    # 基础列表推导式
    sequences = ["ATCG", "GCTA", "TTAA", "CCGG", "ATAT"]
    print(f"原始序列: {sequences}")
    
    # 计算每个序列的长度
    lengths = [len(seq) for seq in sequences]
    print(f"序列长度: {lengths}")
    
    # 计算GC含量
    gc_contents = [(seq.count('G') + seq.count('C')) / len(seq) for seq in sequences]
    print(f"GC含量: {[f'{gc:.2f}' for gc in gc_contents]}")
    
    # 条件过滤 - 高GC含量序列
    high_gc_sequences = [seq for seq in sequences if (seq.count('G') + seq.count('C')) / len(seq) > 0.5]
    print(f"高GC含量序列 (>50%): {high_gc_sequences}")
    
    # 转换为大写
    upper_sequences = [seq.upper() for seq in sequences]
    print(f"大写序列: {upper_sequences}")
    
    # 复杂的列表推导式 - 序列分析
    sequence_analysis = [
        {
            'sequence': seq,
            'length': len(seq),
            'gc_content': (seq.count('G') + seq.count('C')) / len(seq),
            'has_start_codon': 'ATG' in seq
        }
        for seq in sequences
    ]
    
    print(f"\n--- 序列分析结果 ---")
    for analysis in sequence_analysis:
        print(f"序列: {analysis['sequence']}, 长度: {analysis['length']}, "
              f"GC含量: {analysis['gc_content']:.2f}, 含起始密码子: {analysis['has_start_codon']}")

def dictionary_operations():
    """字典操作练习"""
    print("\n=== 字典操作练习 ===")
    
    # 创建基因表达字典
    gene_expression = {
        "BRCA1": 5.2,
        "BRCA2": 3.8,
        "TP53": 7.1,
        "EGFR": 4.5,
        "MYC": 6.3
    }
    
    print(f"基因表达数据: {gene_expression}")
    
    # 基本操作
    print(f"\n--- 基本操作 ---")
    print(f"BRCA1的表达水平: {gene_expression['BRCA1']}")
    print(f"所有基因: {list(gene_expression.keys())}")
    print(f"所有表达水平: {list(gene_expression.values())}")
    print(f"键值对数量: {len(gene_expression)}")
    
    # 安全访问
    print(f"KRAS的表达水平: {gene_expression.get('KRAS', '未检测到')}")
    
    # 添加和修改
    gene_expression["KRAS"] = 4.8
    gene_expression["PIK3CA"] = 5.5
    gene_expression["BRCA1"] = 5.5  # 修改现有值
    print(f"添加新基因后: {gene_expression}")
    
    # 删除
    del gene_expression["PIK3CA"]
    print(f"删除PIK3CA后: {gene_expression}")
    
    # 遍历字典
    print(f"\n--- 遍历字典 ---")
    print("基因表达水平:")
    for gene, expression in gene_expression.items():
        print(f"  {gene}: {expression}")
    
    # 字典方法
    print(f"\n--- 字典方法 ---")
    print(f"所有键: {list(gene_expression.keys())}")
    print(f"所有值: {list(gene_expression.values())}")
    print(f"所有项: {list(gene_expression.items())}")

def dictionary_comprehensions():
    """字典推导式练习"""
    print("\n=== 字典推导式练习 ===")
    
    # 基础字典推导式
    genes = ["BRCA1", "BRCA2", "TP53", "EGFR", "MYC"]
    
    # 创建基因长度字典
    gene_lengths = {gene: len(gene) for gene in genes}
    print(f"基因名长度: {gene_lengths}")
    
    # 条件字典推导式 - 筛选高表达基因
    gene_expression = {"BRCA1": 5.2, "BRCA2": 3.8, "TP53": 7.1, "EGFR": 4.5, "MYC": 6.3}
    high_expression = {gene: expr for gene, expr in gene_expression.items() if expr > 5.0}
    print(f"高表达基因 (>5.0): {high_expression}")
    
    # 转换字典值
    log_expression = {gene: round(expr * 2, 1) for gene, expr in gene_expression.items()}
    print(f"表达水平翻倍: {log_expression}")
    
    # 复杂字典推导式
    gene_categories = {
        gene: "高表达" if expr > 6.0 else "中等表达" if expr > 4.0 else "低表达"
        for gene, expr in gene_expression.items()
    }
    print(f"基因表达分类: {gene_categories}")

def tuple_operations():
    """元组操作练习"""
    print("\n=== 元组操作练习 ===")
    
    # 创建元组
    gene_info = ("BRCA1", 17, 43044295, 43125483, "+")
    print(f"基因信息元组: {gene_info}")
    
    # 元组解包
    name, chromosome, start, end, strand = gene_info
    print(f"基因名: {name}")
    print(f"染色体: {chromosome}")
    print(f"起始位置: {start}")
    print(f"结束位置: {end}")
    print(f"链方向: {strand}")
    
    # 元组索引
    print(f"\n--- 元组索引 ---")
    print(f"第一个元素: {gene_info[0]}")
    print(f"最后一个元素: {gene_info[-1]}")
    print(f"前三个元素: {gene_info[:3]}")
    
    # 命名元组
    from collections import namedtuple
    
    Gene = namedtuple('Gene', ['name', 'chromosome', 'start', 'end', 'strand'])
    brca1 = Gene("BRCA1", 17, 43044295, 43125483, "+")
    brca2 = Gene("BRCA2", 13, 32315086, 32400266, "+")
    
    print(f"\n--- 命名元组 ---")
    print(f"BRCA1基因名: {brca1.name}")
    print(f"BRCA1染色体: {brca1.chromosome}")
    print(f"BRCA2起始位置: {brca2.start}")
    
    # 元组列表
    genes_list = [brca1, brca2]
    print(f"\n基因列表:")
    for gene in genes_list:
        length = gene.end - gene.start + 1
        print(f"{gene.name}: chr{gene.chromosome}:{gene.start}-{gene.end} ({length:,} bp)")

def set_operations():
    """集合操作练习"""
    print("\n=== 集合操作练习 ===")
    
    # 创建集合
    nucleotides = {"A", "T", "C", "G"}
    amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
    
    print(f"DNA核苷酸: {nucleotides}")
    print(f"氨基酸 (前10个): {sorted(list(amino_acids))[:10]}...")
    
    # 集合操作
    purines = {"A", "G"}
    pyrimidines = {"C", "T"}
    
    print(f"\n--- 集合运算 ---")
    print(f"嘌呤: {purines}")
    print(f"嘧啶: {pyrimidines}")
    print(f"并集 (所有碱基): {purines | pyrimidines}")
    print(f"交集: {purines & pyrimidines}")  # 空集
    print(f"差集 (嘌呤-嘧啶): {purines - pyrimidines}")
    print(f"对称差集: {purines ^ pyrimidines}")
    
    # 集合方法
    print(f"\n--- 集合方法 ---")
    all_bases = set()
    all_bases.add("A")
    all_bases.add("T")
    all_bases.update(["C", "G"])
    print(f"构建的碱基集合: {all_bases}")
    
    # 去重应用
    sequence_with_duplicates = ["A", "T", "C", "G", "A", "T", "C", "G", "A"]
    unique_bases = set(sequence_with_duplicates)
    print(f"原序列: {sequence_with_duplicates}")
    print(f"去重后: {unique_bases}")

def nested_data_structures():
    """嵌套数据结构练习"""
    print("\n=== 嵌套数据结构练习 ===")
    
    # 基因数据库
    gene_database = {
        "BRCA1": {
            "chromosome": 17,
            "start": 43044295,
            "end": 43125483,
            "strand": "+",
            "type": "protein_coding",
            "expression": [5.2, 4.8, 5.5, 4.9, 5.1],
            "mutations": ["c.68_69delAG", "c.181T>G", "c.5266dupC"]
        },
        "BRCA2": {
            "chromosome": 13,
            "start": 32315086,
            "end": 32400266,
            "strand": "+",
            "type": "protein_coding",
            "expression": [3.8, 4.1, 3.5, 3.9, 4.0],
            "mutations": ["c.5946delT", "c.9097dupA"]
        },
        "TP53": {
            "chromosome": 17,
            "start": 7661779,
            "end": 7687550,
            "strand": "-",
            "type": "protein_coding",
            "expression": [7.1, 6.8, 7.3, 7.0, 6.9],
            "mutations": ["c.524G>A", "c.743G>A", "c.818G>A"]
        }
    }
    
    print("基因数据库信息:")
    for gene_name, gene_info in gene_database.items():
        print(f"\n{gene_name}:")
        print(f"  位置: chr{gene_info['chromosome']}:{gene_info['start']}-{gene_info['end']}")
        print(f"  链方向: {gene_info['strand']}")
        print(f"  平均表达水平: {sum(gene_info['expression']) / len(gene_info['expression']):.2f}")
        print(f"  突变数量: {len(gene_info['mutations'])}")
    
    # 复杂查询
    print(f"\n--- 复杂查询 ---")
    
    # 查找17号染色体上的基因
    chr17_genes = [name for name, info in gene_database.items() if info['chromosome'] == 17]
    print(f"17号染色体上的基因: {chr17_genes}")
    
    # 查找高表达基因
    high_expr_genes = []
    for name, info in gene_database.items():
        avg_expr = sum(info['expression']) / len(info['expression'])
        if avg_expr > 5.0:
            high_expr_genes.append((name, avg_expr))
    
    print(f"高表达基因 (>5.0):")
    for gene, expr in high_expr_genes:
        print(f"  {gene}: {expr:.2f}")
    
    # 统计突变信息
    total_mutations = sum(len(info['mutations']) for info in gene_database.values())
    print(f"总突变数: {total_mutations}")

def data_analysis_example():
    """数据分析示例"""
    print("\n=== 数据分析示例 ===")
    
    # 模拟基因表达实验数据
    experiment_data = {
        "samples": ["Control_1", "Control_2", "Control_3", "Treatment_1", "Treatment_2", "Treatment_3"],
        "genes": {
            "BRCA1": [5.2, 5.0, 5.1, 7.8, 8.1, 7.9],
            "BRCA2": [3.8, 3.9, 3.7, 4.2, 4.5, 4.3],
            "TP53": [7.1, 7.0, 7.2, 9.5, 9.8, 9.3],
            "EGFR": [4.5, 4.4, 4.6, 4.7, 4.8, 4.6],
            "MYC": [6.3, 6.1, 6.4, 8.9, 9.2, 8.7]
        }
    }
    
    print("实验数据分析:")
    
    # 计算每个基因的统计信息
    for gene, values in experiment_data["genes"].items():
        control_values = values[:3]  # 前3个是对照组
        treatment_values = values[3:]  # 后3个是处理组
        
        control_mean = sum(control_values) / len(control_values)
        treatment_mean = sum(treatment_values) / len(treatment_values)
        fold_change = treatment_mean / control_mean
        
        print(f"\n{gene}:")
        print(f"  对照组平均: {control_mean:.2f}")
        print(f"  处理组平均: {treatment_mean:.2f}")
        print(f"  倍数变化: {fold_change:.2f}")
        
        if fold_change > 1.5:
            print(f"  状态: 上调")
        elif fold_change < 0.67:
            print(f"  状态: 下调")
        else:
            print(f"  状态: 无显著变化")

def exercise_challenges():
    """挑战练习"""
    print("\n=== 挑战练习 ===")
    
    print("挑战1: 序列统计分析器")
    def sequence_statistics(sequences):
        """分析多个序列的统计信息"""
        stats = {
            'total_sequences': len(sequences),
            'total_length': sum(len(seq) for seq in sequences),
            'average_length': sum(len(seq) for seq in sequences) / len(sequences) if sequences else 0,
            'nucleotide_counts': {'A': 0, 'T': 0, 'C': 0, 'G': 0},
            'gc_contents': []
        }
        
        for seq in sequences:
            seq = seq.upper()
            for nucleotide in 'ATCG':
                stats['nucleotide_counts'][nucleotide] += seq.count(nucleotide)
            
            gc_content = (seq.count('G') + seq.count('C')) / len(seq) if len(seq) > 0 else 0
            stats['gc_contents'].append(gc_content)
        
        stats['average_gc'] = sum(stats['gc_contents']) / len(stats['gc_contents']) if stats['gc_contents'] else 0
        
        return stats
    
    test_sequences = [
        "ATCGATCGATCG",
        "GCGCGCGCGCGC", 
        "ATATATATATATAT",
        "CCCGGGAAATTT",
        "TTTTAAAACCCCGGGG"
    ]
    
    stats = sequence_statistics(test_sequences)
    print(f"序列统计结果:")
    print(f"  序列数量: {stats['total_sequences']}")
    print(f"  总长度: {stats['total_length']} bp")
    print(f"  平均长度: {stats['average_length']:.1f} bp")
    print(f"  核苷酸计数: {stats['nucleotide_counts']}")
    print(f"  平均GC含量: {stats['average_gc']:.2%}")
    
    print(f"\n挑战2: 基因表达数据透视")
    # 创建基因表达矩阵
    expression_matrix = [
        ["Gene", "Sample1", "Sample2", "Sample3", "Sample4"],
        ["BRCA1", 5.2, 4.8, 5.5, 4.9],
        ["BRCA2", 3.8, 4.1, 3.5, 3.9],
        ["TP53", 7.1, 6.8, 7.3, 7.0],
        ["EGFR", 4.5, 4.2, 4.8, 4.6]
    ]
    
    # 转换为字典格式
    header = expression_matrix[0]
    data_dict = {}
    
    for row in expression_matrix[1:]:
        gene_name = row[0]
        gene_data = {}
        for i, sample in enumerate(header[1:], 1):
            gene_data[sample] = row[i]
        data_dict[gene_name] = gene_data
    
    print("基因表达矩阵转换结果:")
    for gene, samples in data_dict.items():
        values = list(samples.values())
        mean_expr = sum(values) / len(values)
        print(f"{gene}: 平均表达 = {mean_expr:.2f}, 样本数据 = {samples}")

def main():
    """主函数"""
    print("Python数据结构练习")
    print("=" * 50)
    
    # 执行各个练习
    list_operations()
    list_comprehensions()
    dictionary_operations()
    dictionary_comprehensions()
    tuple_operations()
    set_operations()
    nested_data_structures()
    data_analysis_example()
    exercise_challenges()
    
    print("\n" + "=" * 50)
    print("数据结构练习完成！")
    print("=" * 50)

if __name__ == "__main__":
    main()