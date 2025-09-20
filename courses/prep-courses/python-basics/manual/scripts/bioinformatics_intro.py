#!/usr/bin/env python3
"""
生物信息学Python库介绍脚本
课程：高通量测序数据分析 - Python编程基础
作者：王运生教授
日期：2025年
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# 设置matplotlib中文显示
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

def numpy_for_bioinformatics():
    """NumPy在生物信息学中的应用"""
    print("=== NumPy在生物信息学中的应用 ===")
    
    # 创建基因表达矩阵
    print("--- 基因表达矩阵操作 ---")
    
    # 模拟基因表达数据
    np.random.seed(42)  # 设置随机种子以获得可重现结果
    
    genes = ["BRCA1", "BRCA2", "TP53", "EGFR", "MYC", "KRAS", "PIK3CA"]
    samples = ["Control_1", "Control_2", "Control_3", "Treatment_1", "Treatment_2", "Treatment_3"]
    
    # 创建表达矩阵 (基因 x 样本)
    expression_matrix = np.random.normal(5.0, 1.5, (len(genes), len(samples)))
    expression_matrix = np.abs(expression_matrix)  # 确保表达值为正
    
    # 为处理组添加一些上调的基因
    treatment_indices = [3, 4, 5]  # Treatment样本的索引
    upregulated_genes = [0, 2, 4]  # BRCA1, TP53, MYC的索引
    
    for gene_idx in upregulated_genes:
        for sample_idx in treatment_indices:
            expression_matrix[gene_idx, sample_idx] *= 1.5  # 上调1.5倍
    
    print(f"表达矩阵形状: {expression_matrix.shape}")
    print(f"基因数: {len(genes)}, 样本数: {len(samples)}")
    
    # 基本统计
    print(f"\n--- 基本统计分析 ---")
    print(f"整体平均表达水平: {np.mean(expression_matrix):.2f}")
    print(f"表达水平标准差: {np.std(expression_matrix):.2f}")
    print(f"最高表达水平: {np.max(expression_matrix):.2f}")
    print(f"最低表达水平: {np.min(expression_matrix):.2f}")
    
    # 每个基因的统计
    gene_means = np.mean(expression_matrix, axis=1)
    gene_stds = np.std(expression_matrix, axis=1)
    
    print(f"\n每个基因的平均表达水平:")
    for i, gene in enumerate(genes):
        print(f"{gene:<8}: {gene_means[i]:.2f} ± {gene_stds[i]:.2f}")
    
    # 样本间相关性
    print(f"\n--- 样本间相关性分析 ---")
    correlation_matrix = np.corrcoef(expression_matrix.T)  # 转置以计算样本间相关性
    
    print(f"样本间相关性矩阵:")
    print(f"{'样本':<12}", end="")
    for sample in samples:
        print(f"{sample:>10}", end="")
    print()
    
    for i, sample1 in enumerate(samples):
        print(f"{sample1:<12}", end="")
        for j, sample2 in enumerate(samples):
            print(f"{correlation_matrix[i, j]:>10.3f}", end="")
        print()
    
    # 差异表达分析
    print(f"\n--- 差异表达分析 ---")
    control_samples = expression_matrix[:, :3]  # 前3个样本是对照组
    treatment_samples = expression_matrix[:, 3:]  # 后3个样本是处理组
    
    control_means = np.mean(control_samples, axis=1)
    treatment_means = np.mean(treatment_samples, axis=1)
    fold_changes = treatment_means / control_means
    
    print(f"差异表达基因分析:")
    for i, gene in enumerate(genes):
        fc = fold_changes[i]
        status = "上调" if fc > 1.5 else "下调" if fc < 0.67 else "无显著变化"
        print(f"{gene:<8}: 倍数变化 = {fc:.2f} ({status})")
    
    return expression_matrix, genes, samples

def pandas_for_bioinformatics():
    """Pandas在生物信息学中的应用"""
    print("\n=== Pandas在生物信息学中的应用 ===")
    
    # 创建基因注释DataFrame
    print("--- 基因注释数据处理 ---")
    
    gene_annotations = pd.DataFrame({
        'gene_symbol': ['BRCA1', 'BRCA2', 'TP53', 'EGFR', 'MYC', 'KRAS', 'PIK3CA'],
        'chromosome': ['17', '13', '17', '7', '8', '12', '3'],
        'start_position': [43044295, 32315086, 7661779, 55019017, 127735434, 25205246, 179148114],
        'end_position': [43125483, 32400266, 7687550, 55211628, 127742951, 25250929, 179218294],
        'gene_type': ['tumor_suppressor', 'tumor_suppressor', 'tumor_suppressor', 
                     'growth_factor_receptor', 'oncogene', 'oncogene', 'kinase'],
        'pathway': ['DNA_repair', 'DNA_repair', 'cell_cycle', 'EGFR_signaling', 
                   'MYC_signaling', 'RAS_signaling', 'PI3K_signaling']
    })
    
    # 计算基因长度
    gene_annotations['gene_length'] = gene_annotations['end_position'] - gene_annotations['start_position'] + 1
    
    print("基因注释信息:")
    print(gene_annotations)
    
    # 基本统计
    print(f"\n--- 基本统计信息 ---")
    print(gene_annotations.describe())
    
    # 按基因类型分组
    print(f"\n--- 按基因类型分组统计 ---")
    type_stats = gene_annotations.groupby('gene_type').agg({
        'gene_length': ['count', 'mean', 'std'],
        'chromosome': lambda x: len(set(x))  # 涉及的染色体数
    })
    
    print(type_stats)
    
    # 创建表达数据DataFrame
    print(f"\n--- 表达数据处理 ---")
    
    # 使用之前生成的表达矩阵
    expression_matrix, genes, samples = numpy_for_bioinformatics()
    
    # 创建表达数据DataFrame
    expression_df = pd.DataFrame(
        expression_matrix.T,  # 转置：行为样本，列为基因
        index=samples,
        columns=genes
    )
    
    print("表达数据DataFrame:")
    print(expression_df.head())
    
    # 添加样本信息
    sample_info = pd.DataFrame({
        'sample_id': samples,
        'condition': ['Control', 'Control', 'Control', 'Treatment', 'Treatment', 'Treatment'],
        'batch': ['Batch1', 'Batch1', 'Batch2', 'Batch1', 'Batch2', 'Batch2']
    })
    
    sample_info.set_index('sample_id', inplace=True)
    print(f"\n样本信息:")
    print(sample_info)
    
    # 合并数据进行分析
    print(f"\n--- 条件间比较分析 ---")
    
    # 按条件分组计算统计
    condition_stats = []
    for condition in ['Control', 'Treatment']:
        condition_samples = sample_info[sample_info['condition'] == condition].index
        condition_data = expression_df.loc[condition_samples]
        
        stats = {
            'condition': condition,
            'sample_count': len(condition_data),
            'mean_expression': condition_data.mean().mean(),
            'std_expression': condition_data.std().mean()
        }
        condition_stats.append(stats)
    
    condition_df = pd.DataFrame(condition_stats)
    print(condition_df)
    
    # 差异表达分析
    print(f"\n--- 差异表达基因识别 ---")
    
    control_samples = sample_info[sample_info['condition'] == 'Control'].index
    treatment_samples = sample_info[sample_info['condition'] == 'Treatment'].index
    
    control_means = expression_df.loc[control_samples].mean()
    treatment_means = expression_df.loc[treatment_samples].mean()
    
    diff_expr = pd.DataFrame({
        'gene': genes,
        'control_mean': control_means.values,
        'treatment_mean': treatment_means.values
    })
    
    diff_expr['fold_change'] = diff_expr['treatment_mean'] / diff_expr['control_mean']
    diff_expr['log2_fold_change'] = np.log2(diff_expr['fold_change'])
    
    # 分类差异表达基因
    diff_expr['regulation'] = diff_expr['fold_change'].apply(
        lambda x: 'Upregulated' if x > 1.5 else 'Downregulated' if x < 0.67 else 'No change'
    )
    
    print("差异表达分析结果:")
    print(diff_expr)
    
    return expression_df, gene_annotations, sample_info, diff_expr

def matplotlib_visualization():
    """使用Matplotlib进行生物信息学数据可视化"""
    print("\n=== Matplotlib数据可视化 ===")
    
    # 获取数据
    expression_df, gene_annotations, sample_info, diff_expr = pandas_for_bioinformatics()
    
    # 创建图形目录
    output_dir = Path("data/plots")
    output_dir.mkdir(exist_ok=True)
    
    # 1. 基因表达水平柱状图
    print("--- 生成基因表达柱状图 ---")
    
    plt.figure(figsize=(12, 6))
    
    # 计算每个基因的平均表达
    gene_means = expression_df.mean()
    gene_stds = expression_df.std()
    
    bars = plt.bar(range(len(gene_means)), gene_means.values, 
                   yerr=gene_stds.values, capsize=5,
                   color=['red' if reg == 'Upregulated' else 'blue' if reg == 'Downregulated' else 'gray' 
                         for reg in diff_expr['regulation']])
    
    plt.xlabel('Genes')
    plt.ylabel('Expression Level (log2 FPKM)')
    plt.title('Gene Expression Levels Across All Samples')
    plt.xticks(range(len(gene_means)), gene_means.index, rotation=45)
    
    # 添加图例
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='red', label='Upregulated'),
        Patch(facecolor='blue', label='Downregulated'),
        Patch(facecolor='gray', label='No change')
    ]
    plt.legend(handles=legend_elements)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'gene_expression_barplot.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 2. 样本间表达分布箱线图
    print("--- 生成表达分布箱线图 ---")
    
    plt.figure(figsize=(10, 6))
    
    # 准备数据
    plot_data = []
    plot_labels = []
    plot_colors = []
    
    for sample in expression_df.index:
        plot_data.append(expression_df.loc[sample].values)
        plot_labels.append(sample)
        condition = sample_info.loc[sample, 'condition']
        plot_colors.append('lightblue' if condition == 'Control' else 'lightcoral')
    
    box_plot = plt.boxplot(plot_data, labels=plot_labels, patch_artist=True)
    
    # 设置颜色
    for patch, color in zip(box_plot['boxes'], plot_colors):
        patch.set_facecolor(color)
    
    plt.xlabel('Samples')
    plt.ylabel('Expression Level')
    plt.title('Expression Distribution Across Samples')
    plt.xticks(rotation=45)
    
    # 添加图例
    legend_elements = [
        Patch(facecolor='lightblue', label='Control'),
        Patch(facecolor='lightcoral', label='Treatment')
    ]
    plt.legend(handles=legend_elements)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'expression_distribution_boxplot.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 3. 火山图 (Volcano Plot)
    print("--- 生成火山图 ---")
    
    plt.figure(figsize=(10, 8))
    
    # 计算p值的模拟 (实际应用中需要统计检验)
    np.random.seed(42)
    diff_expr['p_value'] = np.random.uniform(0.001, 0.1, len(diff_expr))
    diff_expr['neg_log10_p'] = -np.log10(diff_expr['p_value'])
    
    # 绘制散点图
    colors = []
    for _, row in diff_expr.iterrows():
        if row['fold_change'] > 1.5 and row['p_value'] < 0.05:
            colors.append('red')  # 显著上调
        elif row['fold_change'] < 0.67 and row['p_value'] < 0.05:
            colors.append('blue')  # 显著下调
        else:
            colors.append('gray')  # 无显著变化
    
    plt.scatter(diff_expr['log2_fold_change'], diff_expr['neg_log10_p'], 
                c=colors, alpha=0.7, s=100)
    
    # 添加基因标签
    for _, row in diff_expr.iterrows():
        if abs(row['log2_fold_change']) > 0.5 and row['neg_log10_p'] > 1:
            plt.annotate(row['gene'], 
                        (row['log2_fold_change'], row['neg_log10_p']),
                        xytext=(5, 5), textcoords='offset points',
                        fontsize=9, alpha=0.8)
    
    # 添加阈值线
    plt.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.5, label='p=0.05')
    plt.axvline(x=np.log2(1.5), color='black', linestyle='--', alpha=0.5, label='FC=1.5')
    plt.axvline(x=np.log2(0.67), color='black', linestyle='--', alpha=0.5, label='FC=0.67')
    
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-Log10 P-value')
    plt.title('Volcano Plot: Differential Gene Expression')
    plt.grid(True, alpha=0.3)
    
    # 添加图例
    legend_elements = [
        Patch(facecolor='red', label='Upregulated'),
        Patch(facecolor='blue', label='Downregulated'),
        Patch(facecolor='gray', label='No change')
    ]
    plt.legend(handles=legend_elements)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'volcano_plot.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"✓ 所有图表已保存到: {output_dir}")

def seaborn_advanced_visualization():
    """使用Seaborn进行高级可视化"""
    print("\n=== Seaborn高级可视化 ===")
    
    # 获取数据
    expression_df, gene_annotations, sample_info, diff_expr = pandas_for_bioinformatics()
    
    # 创建图形目录
    output_dir = Path("data/plots")
    output_dir.mkdir(exist_ok=True)
    
    # 设置seaborn样式
    sns.set_style("whitegrid")
    sns.set_palette("husl")
    
    # 1. 表达数据热图
    print("--- 生成表达热图 ---")
    
    plt.figure(figsize=(10, 8))
    
    # 标准化数据用于热图
    expression_normalized = (expression_df - expression_df.mean()) / expression_df.std()
    
    # 创建注释信息
    condition_colors = sample_info['condition'].map({'Control': 'lightblue', 'Treatment': 'lightcoral'})
    
    # 绘制热图
    sns.heatmap(expression_normalized.T, 
                annot=True, 
                fmt='.2f',
                cmap='RdBu_r',
                center=0,
                cbar_kws={'label': 'Normalized Expression'},
                xticklabels=True,
                yticklabels=True)
    
    plt.title('Gene Expression Heatmap (Normalized)')
    plt.xlabel('Samples')
    plt.ylabel('Genes')
    
    # 添加条件颜色条
    ax = plt.gca()
    for i, (sample, color) in enumerate(condition_colors.items()):
        ax.add_patch(plt.Rectangle((i, len(expression_df.columns)), 1, 0.1, 
                                  facecolor=color, clip_on=False))
    
    plt.tight_layout()
    plt.savefig(output_dir / 'expression_heatmap.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 2. 相关性热图
    print("--- 生成相关性热图 ---")
    
    plt.figure(figsize=(8, 6))
    
    # 计算基因间相关性
    gene_correlation = expression_df.corr()
    
    # 绘制相关性热图
    mask = np.triu(np.ones_like(gene_correlation, dtype=bool))  # 只显示下三角
    
    sns.heatmap(gene_correlation,
                mask=mask,
                annot=True,
                fmt='.3f',
                cmap='coolwarm',
                center=0,
                square=True,
                cbar_kws={'label': 'Correlation Coefficient'})
    
    plt.title('Gene Expression Correlation Matrix')
    plt.tight_layout()
    plt.savefig(output_dir / 'correlation_heatmap.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 3. 分组比较小提琴图
    print("--- 生成分组比较图 ---")
    
    # 准备长格式数据
    melted_data = expression_df.reset_index().melt(
        id_vars=['index'], 
        var_name='gene', 
        value_name='expression'
    )
    melted_data = melted_data.rename(columns={'index': 'sample'})
    
    # 添加条件信息
    melted_data['condition'] = melted_data['sample'].map(
        sample_info['condition'].to_dict()
    )
    
    # 选择几个有代表性的基因
    selected_genes = ['BRCA1', 'TP53', 'MYC', 'EGFR']
    selected_data = melted_data[melted_data['gene'].isin(selected_genes)]
    
    plt.figure(figsize=(12, 6))
    
    sns.violinplot(data=selected_data, 
                   x='gene', 
                   y='expression', 
                   hue='condition',
                   split=True,
                   inner='quart')
    
    plt.title('Gene Expression Distribution by Condition')
    plt.xlabel('Genes')
    plt.ylabel('Expression Level')
    plt.legend(title='Condition')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'violin_plot_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 4. 主成分分析 (PCA) 图
    print("--- 生成PCA图 ---")
    
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    
    # 标准化数据
    scaler = StandardScaler()
    expression_scaled = scaler.fit_transform(expression_df)
    
    # 执行PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(expression_scaled)
    
    # 创建PCA DataFrame
    pca_df = pd.DataFrame({
        'PC1': pca_result[:, 0],
        'PC2': pca_result[:, 1],
        'sample': expression_df.index,
        'condition': [sample_info.loc[sample, 'condition'] for sample in expression_df.index]
    })
    
    plt.figure(figsize=(10, 8))
    
    sns.scatterplot(data=pca_df, 
                    x='PC1', 
                    y='PC2', 
                    hue='condition',
                    s=100,
                    alpha=0.8)
    
    # 添加样本标签
    for _, row in pca_df.iterrows():
        plt.annotate(row['sample'], 
                    (row['PC1'], row['PC2']),
                    xytext=(5, 5), 
                    textcoords='offset points',
                    fontsize=9)
    
    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
    plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
    plt.title('Principal Component Analysis of Gene Expression')
    plt.legend(title='Condition')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'pca_plot.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"PCA解释的方差比例: PC1={pca.explained_variance_ratio_[0]:.1%}, PC2={pca.explained_variance_ratio_[1]:.1%}")

def biopython_demo():
    """BioPython库演示"""
    print("\n=== BioPython库演示 ===")
    
    try:
        from Bio.Seq import Seq
        from Bio.SeqUtils import GC, molecular_weight
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
        from Bio import SeqIO
        from Bio.SeqRecord import SeqRecord
        
        print("✓ BioPython导入成功")
        
        # 1. 序列对象操作
        print("\n--- 序列对象操作 ---")
        
        # DNA序列
        dna_seq = Seq("ATGAAATTTCCCGGGAAATAGATCGATCGATCGATCGATCGATCGATCGATCG")
        
        print(f"DNA序列: {dna_seq}")
        print(f"序列长度: {len(dna_seq)} bp")
        print(f"GC含量: {GC(dna_seq):.2f}%")
        
        # 序列操作
        print(f"反向序列: {dna_seq[::-1]}")
        print(f"互补序列: {dna_seq.complement()}")
        print(f"反向互补: {dna_seq.reverse_complement()}")
        
        # 转录和翻译
        rna_seq = dna_seq.transcribe()
        protein_seq = dna_seq.translate()
        
        print(f"转录RNA: {rna_seq}")
        print(f"翻译蛋白质: {protein_seq}")
        
        # 2. 蛋白质序列分析
        print(f"\n--- 蛋白质序列分析 ---")
        
        # 移除终止密码子
        protein_clean = str(protein_seq).replace('*', '')
        
        if protein_clean:
            protein_analysis = ProteinAnalysis(protein_clean)
            
            print(f"蛋白质序列: {protein_clean}")
            print(f"分子量: {protein_analysis.molecular_weight():.2f} Da")
            print(f"等电点: {protein_analysis.isoelectric_point():.2f}")
            
            # 氨基酸组成
            aa_composition = protein_analysis.get_amino_acids_percent()
            print(f"氨基酸组成 (前5个):")
            for aa, percent in sorted(aa_composition.items(), key=lambda x: x[1], reverse=True)[:5]:
                print(f"  {aa}: {percent:.2%}")
        
        # 3. 创建和写入FASTA文件
        print(f"\n--- FASTA文件操作 ---")
        
        # 创建序列记录
        sequences = [
            SeqRecord(Seq("ATGAAATTTCCCGGGAAATAGATCGATCGATCG"), 
                     id="gene1", 
                     description="BRCA1 partial sequence"),
            SeqRecord(Seq("ATGCCCGGGAAATTTAAATTTCCCGGGAAATAG"), 
                     id="gene2", 
                     description="TP53 partial sequence"),
            SeqRecord(Seq("ATGAAATTTCCCGGGAAATAGATCGATCGATCG"), 
                     id="gene3", 
                     description="EGFR partial sequence")
        ]
        
        # 写入FASTA文件
        output_file = "data/biopython_sequences.fasta"
        with open(output_file, "w") as f:
            SeqIO.write(sequences, f, "fasta")
        
        print(f"✓ 序列已写入: {output_file}")
        
        # 读取FASTA文件
        print(f"\n读取FASTA文件:")
        for record in SeqIO.parse(output_file, "fasta"):
            gc_content = GC(record.seq)
            print(f"ID: {record.id}")
            print(f"描述: {record.description}")
            print(f"长度: {len(record.seq)} bp")
            print(f"GC含量: {gc_content:.2f}%")
            print(f"序列: {record.seq}")
            print("-" * 40)
        
        # 4. 序列搜索和模式匹配
        print(f"\n--- 序列模式搜索 ---")
        
        def find_orfs(sequence, min_length=30):
            """寻找开放阅读框"""
            start_codon = "ATG"
            stop_codons = ["TAA", "TAG", "TGA"]
            orfs = []
            
            sequence_str = str(sequence).upper()
            
            # 搜索三个阅读框
            for frame in range(3):
                for i in range(frame, len(sequence_str) - 2, 3):
                    codon = sequence_str[i:i+3]
                    if codon == start_codon:
                        # 寻找终止密码子
                        for j in range(i + 3, len(sequence_str) - 2, 3):
                            stop_codon = sequence_str[j:j+3]
                            if stop_codon in stop_codons:
                                orf_length = j - i + 3
                                if orf_length >= min_length:
                                    orfs.append({
                                        'start': i,
                                        'end': j + 3,
                                        'length': orf_length,
                                        'frame': frame + 1,
                                        'sequence': sequence_str[i:j+3]
                                    })
                                break
            
            return orfs
        
        # 测试ORF查找
        test_sequence = Seq("ATGAAATTTCCCGGGAAATAGATCGATCGATCGATCGATCGATCGATCGATCGTAG")
        orfs = find_orfs(test_sequence)
        
        print(f"在序列中找到 {len(orfs)} 个ORF:")
        for i, orf in enumerate(orfs, 1):
            print(f"ORF {i}: 位置 {orf['start']}-{orf['end']}, "
                  f"长度 {orf['length']} bp, 阅读框 {orf['frame']}")
    
    except ImportError:
        print("✗ BioPython未安装")
        print("可以使用以下命令安装: pip install biopython")
        
        # 提供简化版本的序列分析
        print("\n使用内置函数进行序列分析:")
        
        def simple_gc_content(sequence):
            """简单的GC含量计算"""
            sequence = sequence.upper()
            gc_count = sequence.count('G') + sequence.count('C')
            return gc_count / len(sequence) if len(sequence) > 0 else 0
        
        def simple_reverse_complement(sequence):
            """简单的反向互补"""
            complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
            sequence = sequence.upper()
            complement = ''.join(complement_map.get(base, base) for base in sequence)
            return complement[::-1]
        
        # 测试简化函数
        test_seq = "ATCGATCGATCG"
        print(f"测试序列: {test_seq}")
        print(f"GC含量: {simple_gc_content(test_seq):.2%}")
        print(f"反向互补: {simple_reverse_complement(test_seq)}")

def practical_example():
    """实际应用示例：基因表达数据分析流水线"""
    print("\n=== 实际应用示例：基因表达分析流水线 ===")
    
    print("--- 步骤1: 数据加载和预处理 ---")
    
    # 模拟从多个文件加载数据
    np.random.seed(42)
    
    # 基因信息
    genes = ["BRCA1", "BRCA2", "TP53", "EGFR", "MYC", "KRAS", "PIK3CA", "AKT1", "PTEN", "RB1"]
    samples = [f"Sample_{i:02d}" for i in range(1, 21)]  # 20个样本
    conditions = ["Control"] * 10 + ["Treatment"] * 10
    
    # 生成表达数据
    expression_data = np.random.normal(5.0, 1.5, (len(genes), len(samples)))
    expression_data = np.abs(expression_data)
    
    # 为处理组的某些基因添加差异表达
    treatment_indices = list(range(10, 20))
    upregulated_genes = [0, 2, 4, 6]  # BRCA1, TP53, MYC, PIK3CA
    downregulated_genes = [1, 3, 9]   # BRCA2, EGFR, RB1
    
    for gene_idx in upregulated_genes:
        for sample_idx in treatment_indices:
            expression_data[gene_idx, sample_idx] *= np.random.uniform(1.5, 2.5)
    
    for gene_idx in downregulated_genes:
        for sample_idx in treatment_indices:
            expression_data[gene_idx, sample_idx] *= np.random.uniform(0.3, 0.6)
    
    # 创建DataFrame
    expression_df = pd.DataFrame(
        expression_data.T,
        index=samples,
        columns=genes
    )
    
    sample_info = pd.DataFrame({
        'sample_id': samples,
        'condition': conditions,
        'batch': [f"Batch{(i//5) + 1}" for i in range(20)]
    }).set_index('sample_id')
    
    print(f"✓ 加载了 {len(genes)} 个基因，{len(samples)} 个样本的表达数据")
    
    print("\n--- 步骤2: 质量控制 ---")
    
    # 基本统计
    print("表达数据基本统计:")
    print(expression_df.describe())
    
    # 检查缺失值
    missing_values = expression_df.isnull().sum().sum()
    print(f"缺失值数量: {missing_values}")
    
    # 样本质量评估
    sample_means = expression_df.mean(axis=1)
    sample_stds = expression_df.std(axis=1)
    
    print(f"样本表达水平范围: {sample_means.min():.2f} - {sample_means.max():.2f}")
    print(f"样本标准差范围: {sample_stds.min():.2f} - {sample_stds.max():.2f}")
    
    print("\n--- 步骤3: 差异表达分析 ---")
    
    # 分组计算
    control_samples = sample_info[sample_info['condition'] == 'Control'].index
    treatment_samples = sample_info[sample_info['condition'] == 'Treatment'].index
    
    control_data = expression_df.loc[control_samples]
    treatment_data = expression_df.loc[treatment_samples]
    
    # 统计检验 (简化版t检验)
    from scipy import stats
    
    diff_results = []
    for gene in genes:
        control_values = control_data[gene].values
        treatment_values = treatment_data[gene].values
        
        # t检验
        t_stat, p_value = stats.ttest_ind(control_values, treatment_values)
        
        # 计算倍数变化
        control_mean = np.mean(control_values)
        treatment_mean = np.mean(treatment_values)
        fold_change = treatment_mean / control_mean if control_mean > 0 else 0
        
        diff_results.append({
            'gene': gene,
            'control_mean': control_mean,
            'treatment_mean': treatment_mean,
            'fold_change': fold_change,
            'log2_fold_change': np.log2(fold_change) if fold_change > 0 else 0,
            'p_value': p_value,
            't_statistic': t_stat
        })
    
    diff_df = pd.DataFrame(diff_results)
    
    # 多重检验校正 (Bonferroni)
    diff_df['p_adjusted'] = diff_df['p_value'] * len(diff_df)
    diff_df['p_adjusted'] = diff_df['p_adjusted'].clip(upper=1.0)
    
    # 分类基因
    def classify_gene(row):
        if row['p_adjusted'] < 0.05:
            if row['fold_change'] > 1.5:
                return 'Significantly Upregulated'
            elif row['fold_change'] < 0.67:
                return 'Significantly Downregulated'
        return 'Not Significant'
    
    diff_df['classification'] = diff_df.apply(classify_gene, axis=1)
    
    print("差异表达分析结果:")
    print(diff_df[['gene', 'fold_change', 'p_value', 'p_adjusted', 'classification']])
    
    # 统计结果
    classification_counts = diff_df['classification'].value_counts()
    print(f"\n分类统计:")
    for classification, count in classification_counts.items():
        print(f"  {classification}: {count} 个基因")
    
    print("\n--- 步骤4: 结果可视化和保存 ---")
    
    # 保存结果
    output_dir = Path("data/analysis_results")
    output_dir.mkdir(exist_ok=True)
    
    # 保存差异表达结果
    diff_df.to_csv(output_dir / "differential_expression_results.csv", index=False)
    
    # 保存表达矩阵
    expression_df.to_csv(output_dir / "expression_matrix.csv")
    
    # 保存样本信息
    sample_info.to_csv(output_dir / "sample_information.csv")
    
    print(f"✓ 分析结果已保存到: {output_dir}")
    
    # 生成简单的汇总报告
    report_content = f"""
基因表达差异分析报告
=====================

分析概况:
- 基因数量: {len(genes)}
- 样本数量: {len(samples)}
- 对照组样本: {len(control_samples)}
- 处理组样本: {len(treatment_samples)}

差异表达基因统计:
- 显著上调基因: {classification_counts.get('Significantly Upregulated', 0)}
- 显著下调基因: {classification_counts.get('Significantly Downregulated', 0)}
- 无显著变化基因: {classification_counts.get('Not Significant', 0)}

显著差异表达基因列表:
{diff_df[diff_df['classification'] != 'Not Significant'][['gene', 'fold_change', 'p_adjusted']].to_string(index=False)}

分析参数:
- 显著性阈值: p < 0.05 (Bonferroni校正)
- 倍数变化阈值: FC > 1.5 或 FC < 0.67
- 统计方法: 双样本t检验
"""
    
    with open(output_dir / "analysis_report.txt", "w") as f:
        f.write(report_content)
    
    print(f"✓ 分析报告已生成: {output_dir / 'analysis_report.txt'}")

def main():
    """主函数"""
    print("生物信息学Python库介绍")
    print("=" * 50)
    
    # 执行各个演示
    numpy_for_bioinformatics()
    pandas_for_bioinformatics()
    matplotlib_visualization()
    seaborn_advanced_visualization()
    biopython_demo()
    practical_example()
    
    print("\n" + "=" * 50)
    print("生物信息学Python库介绍完成！")
    print("=" * 50)
    
    print("\n推荐进一步学习的库:")
    print("- scikit-learn: 机器学习")
    print("- scipy: 科学计算")
    print("- plotly: 交互式可视化")
    print("- networkx: 网络分析")
    print("- pysam: SAM/BAM文件处理")
    print("- pyvcf: VCF文件处理")

if __name__ == "__main__":
    main()