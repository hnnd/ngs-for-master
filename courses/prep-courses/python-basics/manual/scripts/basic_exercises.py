#!/usr/bin/env python3
"""
Python基础语法练习脚本
课程：高通量测序数据分析 - Python编程基础
作者：王运生教授
日期：2025年
"""

def basic_data_types():
    """基础数据类型练习"""
    print("=== 基础数据类型练习 ===")
    
    # 基础数据类型
    gene_name = "BRCA1"
    chromosome = 17
    gc_content = 0.42
    is_tumor_suppressor = True
    
    print(f"基因名: {gene_name} (类型: {type(gene_name).__name__})")
    print(f"染色体: {chromosome} (类型: {type(chromosome).__name__})")
    print(f"GC含量: {gc_content:.2%} (类型: {type(gc_content).__name__})")
    print(f"是否为抑癌基因: {is_tumor_suppressor} (类型: {type(is_tumor_suppressor).__name__})")
    
    # 类型转换练习
    print("\n--- 类型转换练习 ---")
    str_number = "123"
    str_float = "3.14"
    
    print(f"字符串 '{str_number}' 转为整数: {int(str_number)}")
    print(f"字符串 '{str_float}' 转为浮点数: {float(str_float)}")
    print(f"整数 {chromosome} 转为字符串: '{str(chromosome)}'")
    print(f"浮点数 {gc_content} 转为字符串: '{str(gc_content)}'")

def string_operations():
    """字符串操作练习"""
    print("\n=== 字符串操作练习 ===")
    
    # DNA序列操作
    dna_sequence = "ATCGATCGATCG"
    
    print(f"原始序列: {dna_sequence}")
    print(f"序列长度: {len(dna_sequence)}")
    print(f"大写序列: {dna_sequence.upper()}")
    print(f"小写序列: {dna_sequence.lower()}")
    
    # 核苷酸计数
    print(f"\n--- 核苷酸计数 ---")
    print(f"A的数量: {dna_sequence.count('A')}")
    print(f"T的数量: {dna_sequence.count('T')}")
    print(f"C的数量: {dna_sequence.count('C')}")
    print(f"G的数量: {dna_sequence.count('G')}")
    
    # 字符串方法
    print(f"\n--- 字符串方法 ---")
    print(f"查找'CG'的位置: {dna_sequence.find('CG')}")
    print(f"替换A为T: {dna_sequence.replace('A', 'T')}")
    print(f"按'CG'分割: {dna_sequence.split('CG')}")
    
    # 字符串格式化
    print(f"\n--- 字符串格式化 ---")
    gene_id = "ENSG00000139618"
    start_pos = 43044295
    end_pos = 43125483
    
    # f-string格式化
    location = f"基因{gene_name}位于{chromosome}号染色体{start_pos}-{end_pos}"
    print(f"f-string: {location}")
    
    # format方法
    location2 = "基因{}位于{}号染色体{}-{}".format(gene_name, chromosome, start_pos, end_pos)
    print(f"format方法: {location2}")

def analyze_sequence(sequence):
    """分析DNA序列的基本信息"""
    sequence = sequence.upper()
    length = len(sequence)
    
    # 计算各核苷酸数量
    a_count = sequence.count('A')
    t_count = sequence.count('T')
    c_count = sequence.count('C')
    g_count = sequence.count('G')
    
    # 计算GC含量
    gc_content = (g_count + c_count) / length if length > 0 else 0
    
    # 计算AT含量
    at_content = (a_count + t_count) / length if length > 0 else 0
    
    return {
        'sequence': sequence,
        'length': length,
        'A': a_count,
        'T': t_count,
        'C': c_count,
        'G': g_count,
        'GC_content': gc_content,
        'AT_content': at_content,
        'GC_percentage': gc_content * 100,
        'AT_percentage': at_content * 100
    }

def sequence_analysis_demo():
    """序列分析演示"""
    print("\n=== 序列分析演示 ===")
    
    # 测试序列
    test_sequences = [
        "ATCGATCGATCG",
        "GCGCGCGCGCGC",
        "ATATATATATATAT",
        "CCCGGGAAATTT"
    ]
    
    for i, seq in enumerate(test_sequences, 1):
        print(f"\n--- 序列 {i} ---")
        result = analyze_sequence(seq)
        
        print(f"序列: {result['sequence']}")
        print(f"长度: {result['length']} bp")
        print(f"核苷酸组成: A={result['A']}, T={result['T']}, C={result['C']}, G={result['G']}")
        print(f"GC含量: {result['GC_percentage']:.1f}%")
        print(f"AT含量: {result['AT_percentage']:.1f}%")

def string_formatting_advanced():
    """高级字符串格式化"""
    print("\n=== 高级字符串格式化 ===")
    
    # 生物信息学数据格式化
    genes_data = [
        {"name": "BRCA1", "chr": 17, "start": 43044295, "end": 43125483, "expression": 5.234},
        {"name": "BRCA2", "chr": 13, "start": 32315086, "end": 32400266, "expression": 3.876},
        {"name": "TP53", "chr": 17, "start": 7661779, "end": 7687550, "expression": 7.123}
    ]
    
    print("基因信息表:")
    print(f"{'基因名':<8} {'染色体':<6} {'起始位置':<12} {'结束位置':<12} {'表达水平':<8}")
    print("-" * 55)
    
    for gene in genes_data:
        print(f"{gene['name']:<8} {gene['chr']:<6} {gene['start']:<12} {gene['end']:<12} {gene['expression']:<8.2f}")
    
    # 序列格式化输出
    print(f"\n--- FASTA格式输出 ---")
    sequence = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    
    # 模拟FASTA格式
    fasta_header = f">gene1|BRCA1|length={len(sequence)}"
    print(fasta_header)
    
    # 每行60个字符
    line_length = 60
    for i in range(0, len(sequence), line_length):
        print(sequence[i:i+line_length])

def exercise_problems():
    """练习题"""
    print("\n=== 练习题 ===")
    
    print("练习1: 计算反向互补序列")
    def reverse_complement(sequence):
        """计算DNA序列的反向互补序列"""
        complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        
        # 方法1：使用循环
        reverse_seq = sequence[::-1]  # 反向
        complement_seq = ""
        for base in reverse_seq:
            complement_seq += complement_map.get(base, base)
        
        return complement_seq
    
    test_seq = "ATCGATCG"
    rev_comp = reverse_complement(test_seq)
    print(f"原序列: {test_seq}")
    print(f"反向互补: {rev_comp}")
    
    print(f"\n练习2: 寻找序列中的模式")
    def find_pattern(sequence, pattern):
        """在序列中寻找特定模式"""
        positions = []
        start = 0
        
        while True:
            pos = sequence.find(pattern, start)
            if pos == -1:
                break
            positions.append(pos)
            start = pos + 1
        
        return positions
    
    long_sequence = "ATCGATCGATCGAAATCGATCGCGATCG"
    pattern = "CG"
    positions = find_pattern(long_sequence, pattern)
    print(f"序列: {long_sequence}")
    print(f"模式 '{pattern}' 出现的位置: {positions}")
    
    print(f"\n练习3: 翻译密码子")
    def translate_codon(codon):
        """翻译单个密码子为氨基酸"""
        # 简化的遗传密码表
        codon_table = {
            'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
            'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
            'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
            'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
            'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
            'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
            'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
            'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }
        
        return codon_table.get(codon.upper(), 'X')  # X表示未知氨基酸
    
    test_codons = ["ATG", "TTC", "AAA", "TAA", "GGG"]
    print("密码子翻译:")
    for codon in test_codons:
        amino_acid = translate_codon(codon)
        print(f"{codon} -> {amino_acid}")

def main():
    """主函数"""
    print("Python基础语法练习")
    print("=" * 50)
    
    # 执行各个练习
    basic_data_types()
    string_operations()
    sequence_analysis_demo()
    string_formatting_advanced()
    exercise_problems()
    
    print("\n" + "=" * 50)
    print("基础语法练习完成！")
    print("=" * 50)

if __name__ == "__main__":
    main()