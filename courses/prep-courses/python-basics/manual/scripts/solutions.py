#!/usr/bin/env python3
"""
Python编程基础练习题解答
课程：高通量测序数据分析 - Python编程基础
作者：王运生教授
日期：2025年
"""

def exercise_1_sequence_analysis():
    """练习1: 序列分析工具"""
    print("=== 练习1: 序列分析工具 ===")
    
    class DNASequenceAnalyzer:
        """DNA序列分析工具类"""
        
        def __init__(self, sequence):
            """初始化序列分析器"""
            self.sequence = sequence.upper().replace(' ', '').replace('\n', '')
            self.valid_bases = set('ATCG')
            
            # 验证序列
            if not all(base in self.valid_bases for base in self.sequence):
                raise ValueError("序列包含无效字符，只允许A、T、C、G")
        
        def gc_content(self):
            """计算GC含量"""
            if len(self.sequence) == 0:
                return 0.0
            
            gc_count = self.sequence.count('G') + self.sequence.count('C')
            return gc_count / len(self.sequence)
        
        def nucleotide_composition(self):
            """计算核苷酸组成"""
            composition = {}
            for base in 'ATCG':
                composition[base] = self.sequence.count(base)
            
            # 计算百分比
            total = len(self.sequence)
            if total > 0:
                composition_percent = {base: count/total for base, count in composition.items()}
            else:
                composition_percent = {base: 0.0 for base in 'ATCG'}
            
            return composition, composition_percent
        
        def reverse_complement(self):
            """计算反向互补序列"""
            complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
            
            # 先互补，再反向
            complement = ''.join(complement_map[base] for base in self.sequence)
            return complement[::-1]
        
        def find_orfs(self, min_length=100):
            """寻找开放阅读框（ORF）"""
            start_codon = "ATG"
            stop_codons = {"TAA", "TAG", "TGA"}
            orfs = []
            
            # 搜索三个阅读框
            for frame in range(3):
                i = frame
                while i < len(self.sequence) - 2:
                    codon = self.sequence[i:i+3]
                    
                    if codon == start_codon:
                        # 寻找终止密码子
                        j = i + 3
                        while j < len(self.sequence) - 2:
                            stop_codon = self.sequence[j:j+3]
                            if stop_codon in stop_codons:
                                orf_length = j - i + 3
                                if orf_length >= min_length:
                                    orfs.append({
                                        'start': i + 1,  # 1-based position
                                        'end': j + 3,
                                        'length': orf_length,
                                        'frame': frame + 1,
                                        'sequence': self.sequence[i:j+3]
                                    })
                                break
                            j += 3
                        
                        # 如果没有找到终止密码子，检查是否到达序列末尾
                        if j >= len(self.sequence) - 2:
                            orf_length = len(self.sequence) - i
                            if orf_length >= min_length:
                                orfs.append({
                                    'start': i + 1,
                                    'end': len(self.sequence),
                                    'length': orf_length,
                                    'frame': frame + 1,
                                    'sequence': self.sequence[i:],
                                    'incomplete': True
                                })
                    
                    i += 3
            
            return orfs
        
        def codon_usage(self):
            """分析密码子使用频率"""
            codon_counts = {}
            
            # 统计所有可能的密码子
            for i in range(0, len(self.sequence) - 2, 3):
                codon = self.sequence[i:i+3]
                if len(codon) == 3:  # 确保是完整的密码子
                    codon_counts[codon] = codon_counts.get(codon, 0) + 1
            
            # 计算频率
            total_codons = sum(codon_counts.values())
            codon_frequencies = {}
            if total_codons > 0:
                codon_frequencies = {codon: count/total_codons 
                                   for codon, count in codon_counts.items()}
            
            return codon_counts, codon_frequencies
        
        def translate(self, frame=1):
            """翻译DNA序列为蛋白质序列"""
            # 标准遗传密码表
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
            
            protein = ""
            start_pos = frame - 1  # 转换为0-based索引
            
            for i in range(start_pos, len(self.sequence) - 2, 3):
                codon = self.sequence[i:i+3]
                if len(codon) == 3:
                    amino_acid = codon_table.get(codon, 'X')  # X表示未知氨基酸
                    protein += amino_acid
                    if amino_acid == '*':  # 遇到终止密码子
                        break
            
            return protein
        
        def generate_report(self):
            """生成完整的分析报告"""
            composition, composition_percent = self.nucleotide_composition()
            orfs = self.find_orfs(min_length=50)  # 降低最小长度以便演示
            codon_counts, codon_frequencies = self.codon_usage()
            
            report = f"""
DNA序列分析报告
===============

基本信息:
- 序列长度: {len(self.sequence)} bp
- GC含量: {self.gc_content():.2%}

核苷酸组成:
- A: {composition['A']} ({composition_percent['A']:.1%})
- T: {composition['T']} ({composition_percent['T']:.1%})
- C: {composition['C']} ({composition_percent['C']:.1%})
- G: {composition['G']} ({composition_percent['G']:.1%})

开放阅读框 (ORF):
- 找到 {len(orfs)} 个ORF (最小长度50bp)
"""
            
            for i, orf in enumerate(orfs, 1):
                incomplete = " (不完整)" if orf.get('incomplete', False) else ""
                report += f"  ORF{i}: 位置{orf['start']}-{orf['end']}, 长度{orf['length']}bp, 阅读框{orf['frame']}{incomplete}\n"
            
            # 显示最常见的密码子
            if codon_frequencies:
                top_codons = sorted(codon_frequencies.items(), key=lambda x: x[1], reverse=True)[:5]
                report += f"\n最常见的密码子:\n"
                for codon, freq in top_codons:
                    report += f"  {codon}: {freq:.2%}\n"
            
            # 三个阅读框的翻译
            report += f"\n翻译结果:\n"
            for frame in range(1, 4):
                protein = self.translate(frame)
                if len(protein) > 0:
                    display_protein = protein[:50] + "..." if len(protein) > 50 else protein
                    report += f"  阅读框{frame}: {display_protein}\n"
            
            return report
    
    # 测试序列分析器
    print("测试DNA序列分析器:")
    
    # 测试序列（包含多个ORF）
    test_sequence = """
    ATGAAATTTCCCGGGAAATAGATCGATCGATCGATCGATCGATCGATCGATCG
    ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
    CCCGGGAAATTTAAATTTCCCGGGAAATAGATCGATCGATCGATCGATCGAT
    GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
    ATGCCCGGGAAATTTAAATTTCCCGGGAAATAGATCGATCGATCGATCGATC
    GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
    AAATTTCCCGGGAAATAGATCGATCGATCGATCGATCGATCGATCGATCGAT
    CTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
    TAAATGAAATTTCCCGGGAAATAGATCGATCGATCGATCGATCGATCGATCG
    """
    
    try:
        analyzer = DNASequenceAnalyzer(test_sequence)
        print(analyzer.generate_report())
        
        # 测试反向互补
        print(f"原序列前50bp: {analyzer.sequence[:50]}")
        print(f"反向互补前50bp: {analyzer.reverse_complement()[:50]}")
        
    except ValueError as e:
        print(f"错误: {e}")

def exercise_2_expression_analysis():
    """练习2: 基因表达数据分析"""
    print("\n=== 练习2: 基因表达数据分析 ===")
    
    import numpy as np
    import pandas as pd
    from scipy import stats
    import matplotlib.pyplot as plt
    
    class ExpressionAnalyzer:
        """基因表达数据分析器"""
        
        def __init__(self, expression_file, sample_info_file=None):
            """初始化表达分析器"""
            self.expression_df = pd.read_csv(expression_file, index_col=0)
            
            if sample_info_file:
                self.sample_info = pd.read_csv(sample_info_file, index_col=0)
            else:
                # 如果没有样本信息，从样本名推断
                self.sample_info = self._infer_sample_info()
        
        def _infer_sample_info(self):
            """从样本名推断样本信息"""
            sample_info = pd.DataFrame(index=self.expression_df.index)
            
            # 简单的推断逻辑
            sample_info['condition'] = ['Control' if 'control' in sample.lower() or 'ctrl' in sample.lower() 
                                      else 'Treatment' for sample in self.expression_df.index]
            
            return sample_info
        
        def quality_control(self):
            """质量控制分析"""
            print("--- 质量控制分析 ---")
            
            # 基本统计
            print("表达数据基本统计:")
            print(self.expression_df.describe())
            
            # 检查缺失值
            missing_count = self.expression_df.isnull().sum().sum()
            print(f"\n缺失值总数: {missing_count}")
            
            if missing_count > 0:
                print("每个基因的缺失值数量:")
                missing_per_gene = self.expression_df.isnull().sum()
                print(missing_per_gene[missing_per_gene > 0])
            
            # 样本质量评估
            sample_stats = pd.DataFrame({
                'mean_expression': self.expression_df.mean(axis=1),
                'median_expression': self.expression_df.median(axis=1),
                'std_expression': self.expression_df.std(axis=1),
                'expressed_genes': (self.expression_df > 0).sum(axis=1)
            })
            
            print(f"\n样本质量统计:")
            print(sample_stats.describe())
            
            # 识别异常样本
            mean_threshold = sample_stats['mean_expression'].mean() + 2 * sample_stats['mean_expression'].std()
            outlier_samples = sample_stats[sample_stats['mean_expression'] > mean_threshold].index
            
            if len(outlier_samples) > 0:
                print(f"\n可能的异常样本: {list(outlier_samples)}")
            
            return sample_stats
        
        def correlation_analysis(self):
            """相关性分析"""
            print("\n--- 相关性分析 ---")
            
            # 样本间相关性
            sample_correlation = self.expression_df.T.corr()
            
            print("样本间相关性统计:")
            correlation_values = sample_correlation.values
            # 去除对角线元素
            correlation_values = correlation_values[np.triu_indices_from(correlation_values, k=1)]
            
            print(f"平均相关性: {np.mean(correlation_values):.3f}")
            print(f"相关性范围: {np.min(correlation_values):.3f} - {np.max(correlation_values):.3f}")
            
            # 基因间相关性
            gene_correlation = self.expression_df.corr()
            
            # 找到高度相关的基因对
            high_corr_pairs = []
            for i in range(len(gene_correlation.columns)):
                for j in range(i+1, len(gene_correlation.columns)):
                    corr_value = gene_correlation.iloc[i, j]
                    if abs(corr_value) > 0.8:  # 高相关性阈值
                        high_corr_pairs.append({
                            'gene1': gene_correlation.columns[i],
                            'gene2': gene_correlation.columns[j],
                            'correlation': corr_value
                        })
            
            if high_corr_pairs:
                print(f"\n高度相关的基因对 (|r| > 0.8):")
                for pair in high_corr_pairs:
                    print(f"  {pair['gene1']} - {pair['gene2']}: r = {pair['correlation']:.3f}")
            
            return sample_correlation, gene_correlation
        
        def differential_expression_analysis(self, control_condition='Control', treatment_condition='Treatment'):
            """差异表达分析"""
            print(f"\n--- 差异表达分析 ({control_condition} vs {treatment_condition}) ---")
            
            # 获取样本分组
            control_samples = self.sample_info[self.sample_info['condition'] == control_condition].index
            treatment_samples = self.sample_info[self.sample_info['condition'] == treatment_condition].index
            
            print(f"对照组样本数: {len(control_samples)}")
            print(f"处理组样本数: {len(treatment_samples)}")
            
            if len(control_samples) == 0 or len(treatment_samples) == 0:
                print("错误: 缺少对照组或处理组样本")
                return None
            
            # 差异表达分析
            results = []
            
            for gene in self.expression_df.columns:
                control_values = self.expression_df.loc[control_samples, gene].values
                treatment_values = self.expression_df.loc[treatment_samples, gene].values
                
                # 基本统计
                control_mean = np.mean(control_values)
                treatment_mean = np.mean(treatment_values)
                control_std = np.std(control_values)
                treatment_std = np.std(treatment_values)
                
                # 倍数变化
                fold_change = treatment_mean / control_mean if control_mean > 0 else float('inf')
                log2_fold_change = np.log2(fold_change) if fold_change > 0 and fold_change != float('inf') else 0
                
                # t检验
                try:
                    t_stat, p_value = stats.ttest_ind(control_values, treatment_values)
                except:
                    t_stat, p_value = 0, 1
                
                results.append({
                    'gene': gene,
                    'control_mean': control_mean,
                    'control_std': control_std,
                    'treatment_mean': treatment_mean,
                    'treatment_std': treatment_std,
                    'fold_change': fold_change,
                    'log2_fold_change': log2_fold_change,
                    't_statistic': t_stat,
                    'p_value': p_value
                })
            
            # 创建结果DataFrame
            diff_df = pd.DataFrame(results)
            
            # 多重检验校正 (Bonferroni)
            diff_df['p_adjusted'] = diff_df['p_value'] * len(diff_df)
            diff_df['p_adjusted'] = diff_df['p_adjusted'].clip(upper=1.0)
            
            # 分类基因
            def classify_gene(row):
                if row['p_adjusted'] < 0.05:
                    if row['fold_change'] > 1.5:
                        return 'Upregulated'
                    elif row['fold_change'] < 0.67:
                        return 'Downregulated'
                return 'Not Significant'
            
            diff_df['classification'] = diff_df.apply(classify_gene, axis=1)
            
            # 统计结果
            classification_counts = diff_df['classification'].value_counts()
            print(f"\n差异表达基因统计:")
            for classification, count in classification_counts.items():
                print(f"  {classification}: {count}")
            
            # 显示显著差异的基因
            significant_genes = diff_df[diff_df['classification'] != 'Not Significant'].sort_values('p_adjusted')
            
            if len(significant_genes) > 0:
                print(f"\n显著差异表达基因 (前10个):")
                display_cols = ['gene', 'fold_change', 'log2_fold_change', 'p_value', 'p_adjusted', 'classification']
                print(significant_genes[display_cols].head(10).to_string(index=False))
            
            return diff_df
        
        def clustering_analysis(self):
            """聚类分析"""
            print(f"\n--- 聚类分析 ---")
            
            from sklearn.cluster import KMeans
            from sklearn.preprocessing import StandardScaler
            
            # 标准化数据
            scaler = StandardScaler()
            expression_scaled = scaler.fit_transform(self.expression_df)
            
            # K-means聚类 (样本聚类)
            n_clusters = len(self.sample_info['condition'].unique())
            kmeans = KMeans(n_clusters=n_clusters, random_state=42)
            cluster_labels = kmeans.fit_predict(expression_scaled)
            
            # 添加聚类结果到样本信息
            cluster_df = pd.DataFrame({
                'sample': self.expression_df.index,
                'cluster': cluster_labels,
                'condition': self.sample_info['condition'].values
            })
            
            print("样本聚类结果:")
            print(cluster_df)
            
            # 评估聚类质量
            from sklearn.metrics import adjusted_rand_score
            
            # 将条件转换为数值标签
            condition_labels = pd.Categorical(self.sample_info['condition']).codes
            ari_score = adjusted_rand_score(condition_labels, cluster_labels)
            
            print(f"\n聚类质量评估 (调整兰德指数): {ari_score:.3f}")
            
            return cluster_df
        
        def generate_visualizations(self, output_dir="data/expression_plots"):
            """生成可视化图表"""
            print(f"\n--- 生成可视化图表 ---")
            
            import os
            os.makedirs(output_dir, exist_ok=True)
            
            # 1. 表达分布箱线图
            plt.figure(figsize=(12, 6))
            
            conditions = self.sample_info['condition'].unique()
            plot_data = []
            plot_labels = []
            
            for condition in conditions:
                condition_samples = self.sample_info[self.sample_info['condition'] == condition].index
                condition_data = self.expression_df.loc[condition_samples]
                
                for sample in condition_samples:
                    plot_data.append(condition_data.loc[sample].values)
                    plot_labels.append(f"{sample}\n({condition})")
            
            plt.boxplot(plot_data, labels=plot_labels)
            plt.title('Gene Expression Distribution by Sample')
            plt.xlabel('Samples')
            plt.ylabel('Expression Level')
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.savefig(f"{output_dir}/expression_boxplot.png", dpi=300, bbox_inches='tight')
            plt.show()
            
            # 2. 相关性热图
            sample_corr, gene_corr = self.correlation_analysis()
            
            plt.figure(figsize=(10, 8))
            import seaborn as sns
            sns.heatmap(sample_corr, annot=True, fmt='.3f', cmap='coolwarm', center=0)
            plt.title('Sample Correlation Heatmap')
            plt.tight_layout()
            plt.savefig(f"{output_dir}/sample_correlation_heatmap.png", dpi=300, bbox_inches='tight')
            plt.show()
            
            print(f"✓ 图表已保存到: {output_dir}")
    
    # 创建示例数据并测试
    print("创建示例表达数据...")
    
    # 生成示例数据
    np.random.seed(42)
    genes = ["BRCA1", "BRCA2", "TP53", "EGFR", "MYC", "KRAS", "PIK3CA", "AKT1"]
    samples = ["Control_1", "Control_2", "Control_3", "Treatment_1", "Treatment_2", "Treatment_3"]
    
    # 创建表达矩阵
    expression_data = np.random.normal(5.0, 1.5, (len(samples), len(genes)))
    expression_data = np.abs(expression_data)
    
    # 为处理组添加差异表达
    treatment_indices = [3, 4, 5]
    upregulated_genes = [0, 2, 4]  # BRCA1, TP53, MYC
    
    for gene_idx in upregulated_genes:
        for sample_idx in treatment_indices:
            expression_data[sample_idx, gene_idx] *= 2.0
    
    # 创建DataFrame并保存
    expression_df = pd.DataFrame(expression_data, index=samples, columns=genes)
    expression_df.to_csv("data/test_expression.csv")
    
    # 创建样本信息
    sample_info = pd.DataFrame({
        'condition': ['Control', 'Control', 'Control', 'Treatment', 'Treatment', 'Treatment']
    }, index=samples)
    sample_info.to_csv("data/test_sample_info.csv")
    
    # 测试分析器
    analyzer = ExpressionAnalyzer("data/test_expression.csv", "data/test_sample_info.csv")
    
    # 执行分析
    analyzer.quality_control()
    analyzer.correlation_analysis()
    diff_results = analyzer.differential_expression_analysis()
    analyzer.clustering_analysis()
    analyzer.generate_visualizations()

def exercise_3_file_converter():
    """练习3: 文件格式转换工具"""
    print("\n=== 练习3: 文件格式转换工具 ===")
    
    import csv
    import json
    from pathlib import Path
    
    class FileConverter:
        """生物信息学文件格式转换工具"""
        
        def __init__(self):
            self.supported_formats = ['fasta', 'csv', 'json', 'txt']
        
        def fasta_to_csv(self, fasta_file, csv_file):
            """FASTA格式转CSV格式"""
            sequences = []
            current_id = None
            current_seq = ""
            
            try:
                with open(fasta_file, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if line.startswith('>'):
                            # 保存前一个序列
                            if current_id:
                                sequences.append({
                                    'sequence_id': current_id,
                                    'sequence': current_seq,
                                    'length': len(current_seq),
                                    'gc_content': self._calculate_gc_content(current_seq)
                                })
                            
                            # 开始新序列
                            current_id = line[1:]  # 去掉>符号
                            current_seq = ""
                        else:
                            current_seq += line
                    
                    # 保存最后一个序列
                    if current_id:
                        sequences.append({
                            'sequence_id': current_id,
                            'sequence': current_seq,
                            'length': len(current_seq),
                            'gc_content': self._calculate_gc_content(current_seq)
                        })
                
                # 写入CSV文件
                with open(csv_file, 'w', newline='') as f:
                    if sequences:
                        fieldnames = sequences[0].keys()
                        writer = csv.DictWriter(f, fieldnames=fieldnames)
                        writer.writeheader()
                        writer.writerows(sequences)
                
                print(f"✓ 成功转换 {len(sequences)} 个序列: {fasta_file} -> {csv_file}")
                return True
                
            except Exception as e:
                print(f"✗ 转换失败: {e}")
                return False
        
        def csv_to_fasta(self, csv_file, fasta_file, id_column='sequence_id', seq_column='sequence'):
            """CSV格式转FASTA格式"""
            try:
                sequences_written = 0
                
                with open(csv_file, 'r') as f_in, open(fasta_file, 'w') as f_out:
                    reader = csv.DictReader(f_in)
                    
                    for row in reader:
                        seq_id = row.get(id_column, f"seq_{sequences_written + 1}")
                        sequence = row.get(seq_column, "")
                        
                        if sequence:
                            f_out.write(f">{seq_id}\n")
                            # 每行60个字符
                            for i in range(0, len(sequence), 60):
                                f_out.write(sequence[i:i+60] + "\n")
                            sequences_written += 1
                
                print(f"✓ 成功转换 {sequences_written} 个序列: {csv_file} -> {fasta_file}")
                return True
                
            except Exception as e:
                print(f"✗ 转换失败: {e}")
                return False
        
        def expression_to_json(self, csv_file, json_file):
            """表达数据CSV转JSON格式"""
            try:
                data = {}
                
                with open(csv_file, 'r') as f:
                    reader = csv.DictReader(f)
                    
                    for row in reader:
                        sample_id = row.get('sample_id') or row.get('Sample') or list(row.keys())[0]
                        
                        # 提取数值列
                        expression_data = {}
                        for key, value in row.items():
                            if key != sample_id:
                                try:
                                    expression_data[key] = float(value)
                                except ValueError:
                                    expression_data[key] = value
                        
                        data[sample_id] = expression_data
                
                # 写入JSON文件
                with open(json_file, 'w') as f:
                    json.dump(data, f, indent=2)
                
                print(f"✓ 成功转换表达数据: {csv_file} -> {json_file}")
                return True
                
            except Exception as e:
                print(f"✗ 转换失败: {e}")
                return False
        
        def batch_convert(self, input_dir, output_dir, conversion_type):
            """批量文件转换"""
            input_path = Path(input_dir)
            output_path = Path(output_dir)
            output_path.mkdir(exist_ok=True)
            
            conversion_map = {
                'fasta_to_csv': ('.fasta', '.csv', self.fasta_to_csv),
                'csv_to_fasta': ('.csv', '.fasta', self.csv_to_fasta),
                'csv_to_json': ('.csv', '.json', self.expression_to_json)
            }
            
            if conversion_type not in conversion_map:
                print(f"不支持的转换类型: {conversion_type}")
                return False
            
            input_ext, output_ext, convert_func = conversion_map[conversion_type]
            
            # 查找输入文件
            input_files = list(input_path.glob(f"*{input_ext}"))
            
            if not input_files:
                print(f"在 {input_dir} 中没有找到 {input_ext} 文件")
                return False
            
            success_count = 0
            
            for input_file in input_files:
                output_file = output_path / (input_file.stem + output_ext)
                
                print(f"转换: {input_file.name} -> {output_file.name}")
                
                if convert_func(str(input_file), str(output_file)):
                    success_count += 1
            
            print(f"\n批量转换完成: {success_count}/{len(input_files)} 个文件成功转换")
            return success_count == len(input_files)
        
        def validate_file(self, file_path, file_type):
            """验证文件格式"""
            try:
                if file_type == 'fasta':
                    return self._validate_fasta(file_path)
                elif file_type == 'csv':
                    return self._validate_csv(file_path)
                elif file_type == 'json':
                    return self._validate_json(file_path)
                else:
                    print(f"不支持的文件类型: {file_type}")
                    return False
            except Exception as e:
                print(f"验证文件时出错: {e}")
                return False
        
        def _validate_fasta(self, file_path):
            """验证FASTA文件格式"""
            with open(file_path, 'r') as f:
                lines = f.readlines()
            
            if not lines:
                print("文件为空")
                return False
            
            if not lines[0].startswith('>'):
                print("FASTA文件必须以'>'开头")
                return False
            
            sequence_count = 0
            for line in lines:
                if line.startswith('>'):
                    sequence_count += 1
            
            print(f"✓ FASTA文件验证通过，包含 {sequence_count} 个序列")
            return True
        
        def _validate_csv(self, file_path):
            """验证CSV文件格式"""
            with open(file_path, 'r') as f:
                reader = csv.reader(f)
                rows = list(reader)
            
            if len(rows) < 2:
                print("CSV文件至少需要包含标题行和一行数据")
                return False
            
            header_length = len(rows[0])
            for i, row in enumerate(rows[1:], 2):
                if len(row) != header_length:
                    print(f"第{i}行的列数与标题行不匹配")
                    return False
            
            print(f"✓ CSV文件验证通过，{len(rows)-1} 行数据，{header_length} 列")
            return True
        
        def _validate_json(self, file_path):
            """验证JSON文件格式"""
            with open(file_path, 'r') as f:
                data = json.load(f)
            
            if not isinstance(data, dict):
                print("JSON文件应该包含字典格式的数据")
                return False
            
            print(f"✓ JSON文件验证通过，包含 {len(data)} 个条目")
            return True
        
        def _calculate_gc_content(self, sequence):
            """计算GC含量"""
            if not sequence:
                return 0.0
            
            sequence = sequence.upper()
            gc_count = sequence.count('G') + sequence.count('C')
            return gc_count / len(sequence)
    
    # 测试文件转换器
    print("测试文件格式转换器...")
    
    converter = FileConverter()
    
    # 创建测试目录
    test_dir = Path("data/conversion_test")
    test_dir.mkdir(exist_ok=True)
    
    # 创建测试FASTA文件
    test_fasta = test_dir / "test_sequences.fasta"
    fasta_content = """>gene1|BRCA1
ATGAAATTTCCCGGGAAATAGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>gene2|TP53
ATGCCCGGGAAATTTAAATTTCCCGGGAAATAGATCGATCGATCGATCGATC
GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
>gene3|EGFR
ATGAAATTTCCCGGGAAATAGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"""
    
    with open(test_fasta, 'w') as f:
        f.write(fasta_content)
    
    # 测试FASTA到CSV转换
    test_csv = test_dir / "converted_sequences.csv"
    converter.fasta_to_csv(str(test_fasta), str(test_csv))
    
    # 测试CSV到FASTA转换
    test_fasta_back = test_dir / "converted_back.fasta"
    converter.csv_to_fasta(str(test_csv), str(test_fasta_back))
    
    # 验证文件
    print("\n文件验证结果:")
    converter.validate_file(str(test_fasta), 'fasta')
    converter.validate_file(str(test_csv), 'csv')
    
    print(f"\n✓ 转换测试完成，文件保存在: {test_dir}")

def main():
    """主函数 - 运行所有练习题解答"""
    print("Python编程基础练习题解答")
    print("=" * 50)
    
    try:
        # 创建数据目录
        Path("data").mkdir(exist_ok=True)
        
        # 运行练习题
        exercise_1_sequence_analysis()
        exercise_2_expression_analysis()
        exercise_3_file_converter()
        
        print("\n" + "=" * 50)
        print("所有练习题解答完成！")
        print("=" * 50)
        
        print("\n学习建议:")
        print("1. 仔细阅读每个解答的代码实现")
        print("2. 尝试修改参数和输入数据")
        print("3. 思考如何优化和扩展这些功能")
        print("4. 将这些工具应用到实际的生物信息学问题中")
        
    except Exception as e:
        print(f"运行过程中出现错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()