#!/usr/bin/env python3
"""
测序平台选择决策分析脚本
基于多因素分析帮助选择最适合的测序平台
作者：王运生
日期：2025年
"""

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import argparse

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 平台技术参数
PLATFORM_SPECS = {
    'Illumina': {
        'read_length': 300,
        'accuracy': 99.5,
        'throughput_gb': 6000,
        'cost_per_gb': 25,
        'run_time_hours': 48,
        'complexity': 3,  # 1-5, 5最复杂
        'availability': 5,  # 1-5, 5最容易获得
        'maturity': 5,  # 1-5, 5最成熟
    },
    'PacBio': {
        'read_length': 15000,
        'accuracy': 95.0,
        'throughput_gb': 160,
        'cost_per_gb': 200,
        'run_time_hours': 24,
        'complexity': 4,
        'availability': 3,
        'maturity': 4,
    },
    'Nanopore': {
        'read_length': 30000,
        'accuracy': 92.0,
        'throughput_gb': 50,
        'cost_per_gb': 100,
        'run_time_hours': 48,
        'complexity': 3,
        'availability': 4,
        'maturity': 3,
    }
}

# 应用场景权重配置
APPLICATION_WEIGHTS = {
    'genome_assembly': {
        'read_length': 0.4,
        'accuracy': 0.2,
        'throughput': 0.1,
        'cost': 0.2,
        'time': 0.1
    },
    'rna_seq': {
        'read_length': 0.1,
        'accuracy': 0.3,
        'throughput': 0.3,
        'cost': 0.2,
        'time': 0.1
    },
    'variant_calling': {
        'read_length': 0.2,
        'accuracy': 0.4,
        'throughput': 0.2,
        'cost': 0.1,
        'time': 0.1
    },
    'metagenomics': {
        'read_length': 0.3,
        'accuracy': 0.2,
        'throughput': 0.2,
        'cost': 0.2,
        'time': 0.1
    },
    'structural_variation': {
        'read_length': 0.5,
        'accuracy': 0.2,
        'throughput': 0.1,
        'cost': 0.1,
        'time': 0.1
    }
}

def normalize_scores(platforms_dict, metric, higher_better=True):
    """标准化评分到0-100范围"""
    values = [platforms_dict[p][metric] for p in platforms_dict.keys()]
    min_val, max_val = min(values), max(values)
    
    normalized = {}
    for platform in platforms_dict.keys():
        value = platforms_dict[platform][metric]
        if higher_better:
            score = ((value - min_val) / (max_val - min_val)) * 100
        else:
            score = ((max_val - value) / (max_val - min_val)) * 100
        normalized[platform] = score
    
    return normalized

def calculate_platform_scores(application='genome_assembly'):
    """计算各平台在特定应用场景下的综合评分"""
    if application not in APPLICATION_WEIGHTS:
        raise ValueError(f"不支持的应用场景: {application}")
    
    weights = APPLICATION_WEIGHTS[application]
    
    # 标准化各项指标
    read_length_scores = normalize_scores(PLATFORM_SPECS, 'read_length', True)
    accuracy_scores = normalize_scores(PLATFORM_SPECS, 'accuracy', True)
    throughput_scores = normalize_scores(PLATFORM_SPECS, 'throughput_gb', True)
    cost_scores = normalize_scores(PLATFORM_SPECS, 'cost_per_gb', False)  # 成本越低越好
    time_scores = normalize_scores(PLATFORM_SPECS, 'run_time_hours', False)  # 时间越短越好
    
    # 计算加权综合评分
    platform_scores = {}
    detailed_scores = {}
    
    for platform in PLATFORM_SPECS.keys():
        score = (
            read_length_scores[platform] * weights['read_length'] +
            accuracy_scores[platform] * weights['accuracy'] +
            throughput_scores[platform] * weights['throughput'] +
            cost_scores[platform] * weights['cost'] +
            time_scores[platform] * weights['time']
        )
        
        platform_scores[platform] = score
        detailed_scores[platform] = {
            'read_length': read_length_scores[platform],
            'accuracy': accuracy_scores[platform],
            'throughput': throughput_scores[platform],
            'cost': cost_scores[platform],
            'time': time_scores[platform],
            'total': score
        }
    
    return platform_scores, detailed_scores

def create_decision_matrix_plot(detailed_scores, application):
    """创建决策矩阵可视化"""
    # 准备数据
    platforms = list(detailed_scores.keys())
    metrics = ['read_length', 'accuracy', 'throughput', 'cost', 'time']
    metric_names = ['读长', '准确性', '通量', '成本', '时间']
    
    scores_matrix = []
    for platform in platforms:
        scores_matrix.append([detailed_scores[platform][metric] for metric in metrics])
    
    scores_array = np.array(scores_matrix)
    
    # 创建热图
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # 子图1：决策矩阵热图
    im = ax1.imshow(scores_array, cmap='RdYlGn', aspect='auto', vmin=0, vmax=100)
    
    ax1.set_xticks(range(len(metric_names)))
    ax1.set_xticklabels(metric_names)
    ax1.set_yticks(range(len(platforms)))
    ax1.set_yticklabels(platforms)
    ax1.set_title(f'平台评分矩阵 - {application}')
    
    # 添加数值标签
    for i in range(len(platforms)):
        for j in range(len(metrics)):
            text = ax1.text(j, i, f'{scores_array[i, j]:.0f}',
                           ha="center", va="center", color="black", fontweight='bold')
    
    # 添加颜色条
    cbar = plt.colorbar(im, ax=ax1)
    cbar.set_label('评分 (0-100)')
    
    # 子图2：综合评分柱状图
    total_scores = [detailed_scores[p]['total'] for p in platforms]
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
    
    bars = ax2.bar(platforms, total_scores, color=colors, alpha=0.8)
    ax2.set_ylabel('综合评分')
    ax2.set_title(f'平台综合评分 - {application}')
    ax2.grid(True, alpha=0.3)
    
    # 添加数值标签
    for bar, score in zip(bars, total_scores):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                f'{score:.1f}', ha='center', va='bottom', fontweight='bold')
    
    plt.tight_layout()
    
    # 保存图片
    os.makedirs('../results', exist_ok=True)
    plt.savefig(f'../results/decision_matrix_{application}.png', dpi=300, bbox_inches='tight')
    print(f"决策矩阵图已保存到: ../results/decision_matrix_{application}.png")
    
    return fig

def create_radar_chart(detailed_scores):
    """创建雷达图比较各平台"""
    metrics = ['read_length', 'accuracy', 'throughput', 'cost', 'time']
    metric_names = ['读长', '准确性', '通量', '成本', '时间']
    
    # 计算角度
    angles = np.linspace(0, 2 * np.pi, len(metrics), endpoint=False).tolist()
    angles += angles[:1]  # 闭合图形
    
    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(projection='polar'))
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
    
    for i, (platform, scores) in enumerate(detailed_scores.items()):
        values = [scores[metric] for metric in metrics]
        values += values[:1]  # 闭合图形
        
        ax.plot(angles, values, 'o-', linewidth=2, label=platform, color=colors[i])
        ax.fill(angles, values, alpha=0.25, color=colors[i])
    
    # 设置标签
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(metric_names)
    ax.set_ylim(0, 100)
    ax.set_title('平台技术特征雷达图', size=16, fontweight='bold', pad=20)
    ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.0))
    ax.grid(True)
    
    plt.tight_layout()
    plt.savefig('../results/platform_radar_chart.png', dpi=300, bbox_inches='tight')
    print("雷达图已保存到: ../results/platform_radar_chart.png")
    
    return fig

def generate_recommendation_report(platform_scores, detailed_scores, application):
    """生成推荐报告"""
    # 排序平台
    sorted_platforms = sorted(platform_scores.items(), key=lambda x: x[1], reverse=True)
    
    report = f"""
# 测序平台选择推荐报告

## 应用场景: {application}

## 综合评分排名:
"""
    
    for i, (platform, score) in enumerate(sorted_platforms, 1):
        report += f"{i}. **{platform}**: {score:.1f}分\n"
    
    report += "\n## 详细分析:\n\n"
    
    for platform, score in sorted_platforms:
        scores = detailed_scores[platform]
        specs = PLATFORM_SPECS[platform]
        
        report += f"### {platform}\n"
        report += f"- **综合评分**: {score:.1f}/100\n"
        report += f"- **技术参数**:\n"
        report += f"  - 读长: {specs['read_length']:,} bp\n"
        report += f"  - 准确率: {specs['accuracy']:.1f}%\n"
        report += f"  - 通量: {specs['throughput_gb']:,} GB/run\n"
        report += f"  - 成本: ¥{specs['cost_per_gb']}/GB\n"
        report += f"  - 运行时间: {specs['run_time_hours']} 小时\n"
        report += f"- **各项评分**:\n"
        report += f"  - 读长: {scores['read_length']:.1f}/100\n"
        report += f"  - 准确性: {scores['accuracy']:.1f}/100\n"
        report += f"  - 通量: {scores['throughput']:.1f}/100\n"
        report += f"  - 成本: {scores['cost']:.1f}/100\n"
        report += f"  - 时间: {scores['time']:.1f}/100\n\n"
    
    # 推荐建议
    best_platform = sorted_platforms[0][0]
    report += f"## 推荐建议:\n\n"
    report += f"基于当前应用场景({application})的需求分析，推荐使用 **{best_platform}** 平台。\n\n"
    
    # 添加应用场景特定的建议
    if application == 'genome_assembly':
        report += "**基因组组装项目建议**:\n"
        report += "- 优先考虑长读长技术以跨越重复区域\n"
        report += "- 可考虑短读长+长读长的混合策略\n"
        report += "- 注意平衡成本和组装质量\n"
    elif application == 'rna_seq':
        report += "**转录组分析项目建议**:\n"
        report += "- 短读长技术适合基因表达定量\n"
        report += "- 长读长技术适合转录本发现\n"
        report += "- 考虑样本数量对成本的影响\n"
    
    # 保存报告
    os.makedirs('../results', exist_ok=True)
    with open(f'../results/platform_recommendation_{application}.md', 'w', encoding='utf-8') as f:
        f.write(report)
    
    print(f"推荐报告已保存到: ../results/platform_recommendation_{application}.md")
    
    return report

def compare_all_applications():
    """比较所有应用场景下的平台选择"""
    all_results = {}
    
    for application in APPLICATION_WEIGHTS.keys():
        platform_scores, detailed_scores = calculate_platform_scores(application)
        all_results[application] = platform_scores
    
    # 创建应用场景比较图
    df_comparison = pd.DataFrame(all_results).T
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # 创建分组柱状图
    x = np.arange(len(df_comparison.index))
    width = 0.25
    
    platforms = df_comparison.columns
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
    
    for i, platform in enumerate(platforms):
        ax.bar(x + i * width, df_comparison[platform], width, 
               label=platform, color=colors[i], alpha=0.8)
    
    ax.set_xlabel('应用场景')
    ax.set_ylabel('综合评分')
    ax.set_title('不同应用场景下的平台评分比较')
    ax.set_xticks(x + width)
    ax.set_xticklabels(df_comparison.index, rotation=45)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('../results/application_comparison.png', dpi=300, bbox_inches='tight')
    print("应用场景比较图已保存到: ../results/application_comparison.png")
    
    # 保存比较结果
    df_comparison.round(1).to_csv('../results/application_comparison.csv', encoding='utf-8-sig')
    print("应用场景比较结果已保存到: ../results/application_comparison.csv")
    
    return df_comparison

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='测序平台选择决策分析')
    parser.add_argument('--application', type=str, default='genome_assembly',
                       choices=list(APPLICATION_WEIGHTS.keys()),
                       help='应用场景')
    
    args = parser.parse_args()
    
    print("测序平台选择决策分析")
    print("=" * 50)
    print(f"应用场景: {args.application}")
    print()
    
    # 计算平台评分
    platform_scores, detailed_scores = calculate_platform_scores(args.application)
    
    # 显示结果
    print("平台综合评分:")
    sorted_scores = sorted(platform_scores.items(), key=lambda x: x[1], reverse=True)
    for platform, score in sorted_scores:
        print(f"  {platform}: {score:.1f}分")
    
    print(f"\n推荐平台: {sorted_scores[0][0]}")
    
    # 创建可视化
    create_decision_matrix_plot(detailed_scores, args.application)
    create_radar_chart(detailed_scores)
    
    # 生成推荐报告
    generate_recommendation_report(platform_scores, detailed_scores, args.application)
    
    # 比较所有应用场景
    print("\n进行全应用场景比较分析...")
    compare_all_applications()
    
    print("\n" + "=" * 50)
    print("决策分析完成！")
    print("结果文件:")
    print(f"  - 决策矩阵: ../results/decision_matrix_{args.application}.png")
    print("  - 雷达图: ../results/platform_radar_chart.png")
    print(f"  - 推荐报告: ../results/platform_recommendation_{args.application}.md")
    print("  - 应用比较: ../results/application_comparison.png")

if __name__ == "__main__":
    main()