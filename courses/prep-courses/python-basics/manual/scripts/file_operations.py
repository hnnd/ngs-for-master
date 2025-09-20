#!/usr/bin/env python3
"""
Python文件操作练习脚本
课程：高通量测序数据分析 - Python编程基础
作者：王运生教授
日期：2025年
"""

import os
import csv
import json
from pathlib import Path

def create_sample_files():
    """创建示例文件用于练习"""
    print("=== 创建示例文件 ===")
    
    # 确保data目录存在
    data_dir = Path("data")
    data_dir.mkdir(exist_ok=True)
    
    # 创建FASTA文件
    fasta_content = """>gene1|BRCA1|breast_cancer_gene_1
ATGAAATTTCCCGGGAAATAGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
CCCGGGAAATTTAAATTTCCCGGGAAATAGATCGATCGATCGATCGATCGAT
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
>gene2|TP53|tumor_suppressor_protein_53
ATGCCCGGGAAATTTAAATTTCCCGGGAAATAGATCGATCGATCGATCGATC
GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
AAATTTCCCGGGAAATAGATCGATCGATCGATCGATCGATCGATCGATCGAT
CTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
>gene3|EGFR|epidermal_growth_factor_receptor
ATGAAATTTCCCGGGAAATAGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
TTTCCCGGGAAATAGATCGATCGATCGATCGATCGATCGATCGATCGATCGA
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
>gene4|MYC|myelocytomatosis_oncogene
ATGCCCGGGAAATTTAAATTTCCCGGGAAATAGATCGATCGATCGATCGATC
GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
CCCGGGAAATTTAAATTTCCCGGGAAATAGATCGATCGATCGATCGATCGAT
CTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
>gene5|KRAS|kirsten_rat_sarcoma_oncogene
ATGAAATTTCCCGGGAAATAGATCGATCGATCGATCGATCGATCGATCGATCG
CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT
AAATTTCCCGGGAAATAGATCGATCGATCGATCGATCGATCGATCGATCGAT
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT"""
    
    fasta_file = data_dir / "sample_sequences.fasta"
    with open(fasta_file, "w") as f:
        f.write(fasta_content)
    print(f"✓ 创建文件: {fasta_file}")
    
    # 创建基因表达CSV文件
    expression_data = [
        ["gene_name", "sample1", "sample2", "sample3", "sample4", "sample5", "gene_type"],
        ["BRCA1", "5.2", "4.8", "5.5", "4.9", "5.1", "tumor_suppressor"],
        ["BRCA2", "3.8", "4.1", "3.5", "3.9", "4.0", "tumor_suppressor"],
        ["TP53", "7.1", "6.8", "7.3", "7.0", "6.9", "tumor_suppressor"],
        ["EGFR", "4.5", "4.2", "4.8", "4.6", "4.4", "growth_factor"],
        ["MYC", "6.3", "6.0", "6.5", "6.1", "6.2", "oncogene"],
        ["KRAS", "4.8", "5.1", "4.6", "4.9", "5.0", "oncogene"],
        ["PIK3CA", "5.5", "5.2", "5.8", "5.4", "5.6", "kinase"]
    ]
    
    csv_file = data_dir / "gene_expression.csv"
    with open(csv_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(expression_data)
    print(f"✓ 创建文件: {csv_file}")
    
    # 创建质量分数文件
    quality_content = """# 测序质量分数数据
# 每行代表一个读段的质量分数
# 格式：位置1 位置2 位置3 ... 位置N
35 38 40 42 38 35 30 28 25 22 20 18 15 12 10 8 5
40 42 45 43 40 38 35 32 30 28 25 23 20 18 15 12 8
38 40 42 45 42 40 38 35 33 30 28 25 22 20 17 14 10
42 45 47 45 43 40 38 35 32 30 27 25 22 19 16 13 9
36 38 40 42 40 38 35 32 30 27 25 22 20 17 15 12 8
39 41 43 45 43 41 38 36 33 31 28 26 23 21 18 15 11
37 39 41 43 41 39 36 34 31 29 26 24 21 19 16 13 9
41 43 45 47 45 43 40 38 35 33 30 28 25 23 20 17 13
38 40 42 44 42 40 37 35 32 30 27 25 22 20 17 14 10
40 42 44 46 44 42 39 37 34 32 29 27 24 22 19 16 12"""
    
    quality_file = data_dir / "quality_scores.txt"
    with open(quality_file, "w") as f:
        f.write(quality_content)
    print(f"✓ 创建文件: {quality_file}")
    
    # 创建基因注释JSON文件
    annotation_data = {
        "BRCA1": {
            "full_name": "BRCA1 DNA repair associated",
            "chromosome": "17",
            "start": 43044295,
            "end": 43125483,
            "strand": "+",
            "gene_type": "protein_coding",
            "description": "Tumor suppressor gene involved in DNA repair"
        },
        "BRCA2": {
            "full_name": "BRCA2 DNA repair associated", 
            "chromosome": "13",
            "start": 32315086,
            "end": 32400266,
            "strand": "+",
            "gene_type": "protein_coding",
            "description": "Tumor suppressor gene involved in homologous recombination"
        },
        "TP53": {
            "full_name": "tumor protein p53",
            "chromosome": "17", 
            "start": 7661779,
            "end": 7687550,
            "strand": "-",
            "gene_type": "protein_coding",
            "description": "Guardian of the genome, cell cycle checkpoint"
        }
    }
    
    json_file = data_dir / "gene_annotations.json"
    with open(json_file, "w") as f:
        json.dump(annotation_data, f, indent=2)
    print(f"✓ 创建文件: {json_file}")

def basic_file_operations():
    """基础文件操作练习"""
    print("\n=== 基础文件操作练习 ===")
    
    # 读取整个文件
    print("--- 读取整个文件 ---")
    try:
        with open("data/quality_scores.txt", "r") as f:
            content = f.read()
            lines = content.split('\n')
            non_comment_lines = [line for line in lines if not line.startswith('#') and line.strip()]
            print(f"文件总行数: {len(lines)}")
            print(f"数据行数: {len(non_comment_lines)}")
            print(f"前两行数据:")
            for i, line in enumerate(non_comment_lines[:2]):
                print(f"  行{i+1}: {line}")
    except FileNotFoundError:
        print("文件不存在，请先运行create_sample_files()")
    
    # 逐行读取
    print(f"\n--- 逐行读取文件 ---")
    try:
        with open("data/quality_scores.txt", "r") as f:
            line_count = 0
            data_lines = 0
            for line in f:
                line_count += 1
                line = line.strip()
                if line and not line.startswith('#'):
                    data_lines += 1
                    if data_lines <= 3:  # 只显示前3行数据
                        scores = line.split()
                        avg_score = sum(int(score) for score in scores) / len(scores)
                        print(f"  数据行{data_lines}: 平均质量分数 = {avg_score:.1f}")
            
            print(f"总行数: {line_count}, 数据行数: {data_lines}")
    except FileNotFoundError:
        print("文件不存在")
    
    # 写入文件
    print(f"\n--- 写入文件 ---")
    output_file = "data/analysis_results.txt"
    
    results = [
        "基因表达分析结果",
        "=" * 30,
        "高表达基因: BRCA1, TP53, MYC",
        "中等表达基因: EGFR, KRAS, PIK3CA", 
        "低表达基因: BRCA2",
        "",
        "分析完成时间: 2025-01-01 12:00:00"
    ]
    
    with open(output_file, "w") as f:
        for result in results:
            f.write(result + "\n")
    
    print(f"✓ 结果已写入: {output_file}")
    
    # 追加写入
    with open(output_file, "a") as f:
        f.write("\n--- 追加信息 ---\n")
        f.write("样本数量: 5\n")
        f.write("基因数量: 7\n")
    
    print(f"✓ 追加信息已写入: {output_file}")

def fasta_file_processing():
    """FASTA文件处理"""
    print("\n=== FASTA文件处理 ===")
    
    def read_fasta_simple(filename):
        """简单的FASTA文件读取器"""
        sequences = {}
        current_id = None
        
        try:
            with open(filename, "r") as f:
                for line in f:
                    line = line.strip()
                    if line.startswith(">"):
                        # 解析序列ID和描述
                        header_parts = line[1:].split("|")
                        current_id = header_parts[0] if header_parts else line[1:]
                        sequences[current_id] = {
                            'sequence': '',
                            'description': line[1:],
                            'gene_name': header_parts[1] if len(header_parts) > 1 else '',
                            'full_name': header_parts[2] if len(header_parts) > 2 else ''
                        }
                    elif current_id:
                        sequences[current_id]['sequence'] += line
            
            return sequences
        except FileNotFoundError:
            print(f"错误: 文件 {filename} 不存在")
            return {}
    
    def analyze_sequence(sequence):
        """分析单个序列"""
        sequence = sequence.upper()
        length = len(sequence)
        
        if length == 0:
            return None
        
        # 计算核苷酸组成
        composition = {
            'A': sequence.count('A'),
            'T': sequence.count('T'),
            'C': sequence.count('C'),
            'G': sequence.count('G')
        }
        
        # 计算GC含量
        gc_content = (composition['G'] + composition['C']) / length
        
        # 寻找起始密码子
        start_codons = []
        for i in range(0, length - 2, 3):
            codon = sequence[i:i+3]
            if codon == 'ATG':
                start_codons.append(i)
        
        return {
            'length': length,
            'composition': composition,
            'gc_content': gc_content,
            'start_codons': start_codons
        }
    
    # 读取和分析FASTA文件
    print("读取FASTA文件...")
    sequences = read_fasta_simple("data/sample_sequences.fasta")
    
    if sequences:
        print(f"成功读取 {len(sequences)} 个序列")
        
        print(f"\n--- 序列分析结果 ---")
        for seq_id, seq_data in sequences.items():
            analysis = analyze_sequence(seq_data['sequence'])
            if analysis:
                print(f"\n序列ID: {seq_id}")
                print(f"基因名: {seq_data['gene_name']}")
                print(f"长度: {analysis['length']} bp")
                print(f"GC含量: {analysis['gc_content']:.2%}")
                print(f"核苷酸组成: A={analysis['composition']['A']}, "
                      f"T={analysis['composition']['T']}, "
                      f"C={analysis['composition']['C']}, "
                      f"G={analysis['composition']['G']}")
                print(f"起始密码子位置: {analysis['start_codons']}")
                print(f"前50个碱基: {seq_data['sequence'][:50]}...")
    
    # 写入分析结果到新的FASTA文件
    print(f"\n--- 生成分析报告 ---")
    report_file = "data/sequence_analysis_report.txt"
    
    with open(report_file, "w") as f:
        f.write("FASTA序列分析报告\n")
        f.write("=" * 50 + "\n\n")
        
        for seq_id, seq_data in sequences.items():
            analysis = analyze_sequence(seq_data['sequence'])
            if analysis:
                f.write(f"序列ID: {seq_id}\n")
                f.write(f"基因名: {seq_data['gene_name']}\n")
                f.write(f"描述: {seq_data['full_name']}\n")
                f.write(f"长度: {analysis['length']} bp\n")
                f.write(f"GC含量: {analysis['gc_content']:.2%}\n")
                f.write(f"起始密码子数量: {len(analysis['start_codons'])}\n")
                f.write("-" * 30 + "\n")
    
    print(f"✓ 分析报告已保存: {report_file}")

def csv_file_processing():
    """CSV文件处理"""
    print("\n=== CSV文件处理 ===")
    
    # 使用csv模块读取
    print("--- 使用csv模块读取 ---")
    try:
        with open("data/gene_expression.csv", "r") as f:
            reader = csv.DictReader(f)
            
            print(f"列名: {reader.fieldnames}")
            
            gene_data = []
            for row in reader:
                # 转换数值列
                expression_values = []
                for i in range(1, 6):  # sample1 到 sample5
                    sample_key = f"sample{i}"
                    if sample_key in row:
                        expression_values.append(float(row[sample_key]))
                
                gene_info = {
                    'name': row['gene_name'],
                    'type': row['gene_type'],
                    'expressions': expression_values,
                    'mean_expression': sum(expression_values) / len(expression_values) if expression_values else 0
                }
                gene_data.append(gene_info)
            
            print(f"\n基因表达数据:")
            for gene in gene_data:
                print(f"{gene['name']:<8} ({gene['type']:<15}): "
                      f"平均表达 = {gene['mean_expression']:.2f}")
    
    except FileNotFoundError:
        print("CSV文件不存在")
        return
    
    # 数据分析和统计
    print(f"\n--- 数据分析 ---")
    
    # 按基因类型分组
    gene_types = {}
    for gene in gene_data:
        gene_type = gene['type']
        if gene_type not in gene_types:
            gene_types[gene_type] = []
        gene_types[gene_type].append(gene)
    
    print("按基因类型分组:")
    for gene_type, genes in gene_types.items():
        expressions = [gene['mean_expression'] for gene in genes]
        avg_expression = sum(expressions) / len(expressions)
        print(f"  {gene_type}: {len(genes)}个基因, 平均表达 = {avg_expression:.2f}")
    
    # 识别高表达和低表达基因
    all_expressions = [gene['mean_expression'] for gene in gene_data]
    overall_mean = sum(all_expressions) / len(all_expressions)
    
    high_expression_genes = [gene for gene in gene_data if gene['mean_expression'] > overall_mean * 1.2]
    low_expression_genes = [gene for gene in gene_data if gene['mean_expression'] < overall_mean * 0.8]
    
    print(f"\n表达水平分类 (基准: {overall_mean:.2f}):")
    print(f"高表达基因 (>{overall_mean*1.2:.2f}): {[g['name'] for g in high_expression_genes]}")
    print(f"低表达基因 (<{overall_mean*0.8:.2f}): {[g['name'] for g in low_expression_genes]}")
    
    # 写入分析结果到新的CSV文件
    print(f"\n--- 保存分析结果 ---")
    
    output_data = []
    for gene in gene_data:
        # 分类表达水平
        if gene['mean_expression'] > overall_mean * 1.2:
            expression_category = "高表达"
        elif gene['mean_expression'] < overall_mean * 0.8:
            expression_category = "低表达"
        else:
            expression_category = "中等表达"
        
        output_data.append({
            'gene_name': gene['name'],
            'gene_type': gene['type'],
            'mean_expression': round(gene['mean_expression'], 2),
            'expression_category': expression_category,
            'sample_count': len(gene['expressions'])
        })
    
    output_file = "data/expression_analysis.csv"
    with open(output_file, "w", newline="") as f:
        fieldnames = ['gene_name', 'gene_type', 'mean_expression', 'expression_category', 'sample_count']
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        
        writer.writeheader()
        writer.writerows(output_data)
    
    print(f"✓ 分析结果已保存: {output_file}")

def json_file_processing():
    """JSON文件处理"""
    print("\n=== JSON文件处理 ===")
    
    # 读取JSON文件
    print("--- 读取JSON注释文件 ---")
    try:
        with open("data/gene_annotations.json", "r") as f:
            annotations = json.load(f)
        
        print(f"成功读取 {len(annotations)} 个基因的注释信息")
        
        for gene_name, info in annotations.items():
            print(f"\n{gene_name}:")
            print(f"  全名: {info['full_name']}")
            print(f"  位置: chr{info['chromosome']}:{info['start']}-{info['end']}")
            print(f"  链方向: {info['strand']}")
            print(f"  描述: {info['description']}")
            
            # 计算基因长度
            gene_length = info['end'] - info['start'] + 1
            print(f"  长度: {gene_length:,} bp")
    
    except FileNotFoundError:
        print("JSON文件不存在")
        return
    
    # 合并多个数据源
    print(f"\n--- 数据整合 ---")
    
    # 读取表达数据
    expression_data = {}
    try:
        with open("data/gene_expression.csv", "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                gene_name = row['gene_name']
                expressions = []
                for i in range(1, 6):
                    sample_key = f"sample{i}"
                    if sample_key in row:
                        expressions.append(float(row[sample_key]))
                
                expression_data[gene_name] = {
                    'mean_expression': sum(expressions) / len(expressions) if expressions else 0,
                    'gene_type': row['gene_type']
                }
    except FileNotFoundError:
        expression_data = {}
    
    # 整合数据
    integrated_data = {}
    for gene_name in annotations.keys():
        integrated_data[gene_name] = {
            'annotation': annotations[gene_name],
            'expression': expression_data.get(gene_name, {'mean_expression': 0, 'gene_type': 'unknown'})
        }
    
    # 保存整合后的数据
    output_file = "data/integrated_gene_data.json"
    with open(output_file, "w") as f:
        json.dump(integrated_data, f, indent=2)
    
    print(f"✓ 整合数据已保存: {output_file}")
    
    # 生成汇总报告
    print(f"\n--- 生成汇总报告 ---")
    
    report_data = {
        'summary': {
            'total_genes': len(integrated_data),
            'chromosomes': list(set(data['annotation']['chromosome'] for data in integrated_data.values())),
            'gene_types': list(set(data['expression']['gene_type'] for data in integrated_data.values()))
        },
        'statistics': {
            'average_expression': sum(data['expression']['mean_expression'] for data in integrated_data.values()) / len(integrated_data),
            'total_genomic_span': sum(data['annotation']['end'] - data['annotation']['start'] + 1 for data in integrated_data.values())
        }
    }
    
    summary_file = "data/analysis_summary.json"
    with open(summary_file, "w") as f:
        json.dump(report_data, f, indent=2)
    
    print(f"✓ 汇总报告已保存: {summary_file}")
    print(f"  总基因数: {report_data['summary']['total_genes']}")
    print(f"  涉及染色体: {', '.join(report_data['summary']['chromosomes'])}")
    print(f"  平均表达水平: {report_data['statistics']['average_expression']:.2f}")

def file_path_operations():
    """文件路径操作"""
    print("\n=== 文件路径操作 ===")
    
    # 使用pathlib
    print("--- 使用pathlib进行路径操作 ---")
    
    data_dir = Path("data")
    print(f"数据目录: {data_dir}")
    print(f"绝对路径: {data_dir.absolute()}")
    print(f"目录是否存在: {data_dir.exists()}")
    
    # 列出目录中的文件
    if data_dir.exists():
        print(f"\n目录中的文件:")
        for file_path in data_dir.iterdir():
            if file_path.is_file():
                size = file_path.stat().st_size
                print(f"  {file_path.name}: {size} bytes")
    
    # 文件操作
    print(f"\n--- 文件信息 ---")
    
    files_to_check = [
        "data/sample_sequences.fasta",
        "data/gene_expression.csv", 
        "data/quality_scores.txt"
    ]
    
    for file_path_str in files_to_check:
        file_path = Path(file_path_str)
        if file_path.exists():
            stat = file_path.stat()
            print(f"{file_path.name}:")
            print(f"  大小: {stat.st_size} bytes")
            print(f"  扩展名: {file_path.suffix}")
            print(f"  文件名(无扩展名): {file_path.stem}")
        else:
            print(f"{file_path.name}: 文件不存在")

def error_handling_demo():
    """文件操作错误处理演示"""
    print("\n=== 错误处理演示 ===")
    
    def safe_read_file(filename):
        """安全地读取文件"""
        try:
            with open(filename, "r") as f:
                content = f.read()
                return content
        except FileNotFoundError:
            print(f"错误: 文件 {filename} 不存在")
            return None
        except PermissionError:
            print(f"错误: 没有权限读取文件 {filename}")
            return None
        except UnicodeDecodeError:
            print(f"错误: 文件 {filename} 编码格式不支持")
            return None
        except Exception as e:
            print(f"未知错误: {e}")
            return None
    
    def safe_write_file(filename, content):
        """安全地写入文件"""
        try:
            # 确保目录存在
            file_path = Path(filename)
            file_path.parent.mkdir(parents=True, exist_ok=True)
            
            with open(filename, "w") as f:
                f.write(content)
            return True
        except PermissionError:
            print(f"错误: 没有权限写入文件 {filename}")
            return False
        except OSError as e:
            print(f"系统错误: {e}")
            return False
        except Exception as e:
            print(f"未知错误: {e}")
            return False
    
    # 测试错误处理
    print("测试文件读取错误处理:")
    
    # 测试不存在的文件
    content = safe_read_file("nonexistent_file.txt")
    print(f"读取不存在文件的结果: {content}")
    
    # 测试正常文件
    content = safe_read_file("data/quality_scores.txt")
    if content:
        print(f"成功读取文件，内容长度: {len(content)} 字符")
    
    # 测试文件写入
    print(f"\n测试文件写入:")
    test_content = "这是一个测试文件\n包含多行内容\n用于演示错误处理"
    
    success = safe_write_file("data/test_output.txt", test_content)
    if success:
        print("✓ 文件写入成功")
    else:
        print("✗ 文件写入失败")

def batch_file_processing():
    """批量文件处理示例"""
    print("\n=== 批量文件处理示例 ===")
    
    # 创建多个测试文件
    print("--- 创建测试文件 ---")
    
    test_dir = Path("data/batch_test")
    test_dir.mkdir(exist_ok=True)
    
    # 创建多个小的FASTA文件
    sequences_data = {
        "gene_set_1.fasta": [
            (">gene1", "ATCGATCGATCG"),
            (">gene2", "GCTAGCTAGCTA")
        ],
        "gene_set_2.fasta": [
            (">gene3", "TTAACCGGTTAA"),
            (">gene4", "CGATCGATCGAT")
        ],
        "gene_set_3.fasta": [
            (">gene5", "AAATTTCCCGGG"),
            (">gene6", "GGGCCCAAATTT")
        ]
    }
    
    for filename, sequences in sequences_data.items():
        file_path = test_dir / filename
        with open(file_path, "w") as f:
            for header, sequence in sequences:
                f.write(f"{header}\n{sequence}\n")
        print(f"✓ 创建文件: {file_path}")
    
    # 批量处理文件
    print(f"\n--- 批量处理FASTA文件 ---")
    
    def process_fasta_file(file_path):
        """处理单个FASTA文件"""
        sequences = {}
        current_id = None
        
        with open(file_path, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    current_id = line[1:]
                    sequences[current_id] = ""
                elif current_id:
                    sequences[current_id] += line
        
        # 分析序列
        results = []
        for seq_id, sequence in sequences.items():
            gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence) if sequence else 0
            results.append({
                'id': seq_id,
                'length': len(sequence),
                'gc_content': gc_content
            })
        
        return results
    
    # 处理所有FASTA文件
    all_results = []
    fasta_files = list(test_dir.glob("*.fasta"))
    
    for fasta_file in fasta_files:
        print(f"处理文件: {fasta_file.name}")
        results = process_fasta_file(fasta_file)
        
        for result in results:
            result['source_file'] = fasta_file.name
            all_results.append(result)
            print(f"  {result['id']}: 长度={result['length']}, GC含量={result['gc_content']:.2%}")
    
    # 保存汇总结果
    summary_file = test_dir / "batch_analysis_summary.csv"
    with open(summary_file, "w", newline="") as f:
        fieldnames = ['source_file', 'id', 'length', 'gc_content']
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        
        writer.writeheader()
        writer.writerows(all_results)
    
    print(f"\n✓ 批量处理完成，结果保存至: {summary_file}")
    print(f"总共处理了 {len(all_results)} 个序列")

def main():
    """主函数"""
    print("Python文件操作练习")
    print("=" * 50)
    
    # 创建示例文件
    create_sample_files()
    
    # 执行各个练习
    basic_file_operations()
    fasta_file_processing()
    csv_file_processing()
    json_file_processing()
    file_path_operations()
    error_handling_demo()
    batch_file_processing()
    
    print("\n" + "=" * 50)
    print("文件操作练习完成！")
    print("=" * 50)

if __name__ == "__main__":
    main()