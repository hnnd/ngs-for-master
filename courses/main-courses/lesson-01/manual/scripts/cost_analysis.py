#!/usr/bin/env python3
"""
测序成本分析脚本
分析不同测序平台的成本效益
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

# 平台成本参数（2025年估算值）
PLATFORM_COSTS = {
    'Illumina': {
        'equipment_cost': 750000,  # 设备成本（元）
        'per_gb_cost': 25,         # 每GB成本（元）
        'setup_cost': 1000,        # 建库成本（元）
        'maintenance_yearly': 50000, # 年维护成本（元）
        'throughput_per_run': 6000,  # 单次运行产出（GB）
        'run_time_hours': 48,        # 运行时间（小时）
        'accuracy': 0.999            # 准确率
    },
    'PacBio': {
        'equipment_cost': 1200000,
        'per_gb_cost': 200,
        'setup_cost': 3000,
        'maintenance_yearly': 80000,
        'throughput_per_run': 160,
        'run_time_hours': 24,
        'accuracy': 0.95
    },
    'Nanopore': {
        'equipment_cost': 150000,
        'per_gb_cost': 100,
        'setup_cost': 1500,
        'maintenance_yearly': 20000,
        'throughput_per_run': 50,
        'run_time_hours': 48,
        'accuracy': 0.92
    }
}

def calculate_project_cost(platform, data_size_gb, coverage=30):
    """计算项目总成本"""
    params = PLATFORM_COSTS[platform]
    
    # 基础成本计算
    sequencing_cost = data_size_gb * params['per_gb_cost']
    setup_cost = params['setup_cost']
    
    # 运行次数计算
    runs_needed = max(1, np.ceil(data_size_gb / params['throughput_per_run']))
    
    # 时间成本（假设人工成本200元/小时）
    total_time = runs_needed * params['run_time_hours']
    labor_cost = total_time * 200
    
    # 设备折旧（假设5年折旧）
    equipment_depreciation = (params['equipment_cost'] / 5) * (total_time / (365 * 24))
    
    # 维护成本分摊
    maintenance_cost = params['maintenance_yearly'] * (total_time / (365 * 24))
    
    total_cost = sequencing_cost + setup_cost + labor_cost + equipment_depreciation + maintenance_cost
    
    return {
        'platform': platform,
        'data_size_gb': data_size_gb,
        'coverage': coverage,
        'sequencing_cost': sequencing_cost,
        'setup_cost': setup_cost,
        'labor_cost': labor_cost,
        'equipment_cost': equipment_depreciation,
        'maintenance_cost': maintenance_cost,
        'total_cost': total_cost,
        'cost_per_gb': total_cost / data_size_gb,
        'runs_needed': runs_needed,
        'total_time_hours': total_time,
        'accuracy': params['accuracy']
    }

def analyze_cost_scenarios():
    """分析不同场景下的成本"""
    scenarios = [
        {'name': '小型项目', 'size_gb': 10, 'coverage': 30},
        {'name': '中型项目', 'size_gb': 100, 'coverage': 30},
        {'name': '大型项目', 'size_gb': 1000, 'coverage': 30},
        {'name': '人类基因组', 'size_gb': 90, 'coverage': 30},
        {'name': '植物基因组', 'size_gb': 300, 'coverage': 20}
    ]
    
    results = []
    
    for scenario in scenarios:
        print(f"\n分析场景: {scenario['name']}")
        print(f"数据量: {scenario['size_gb']} GB, 覆盖度: {scenario['coverage']}x")
        print("-" * 40)
        
        for platform in PLATFORM_COSTS.keys():
            cost_result = calculate_project_cost(
                platform, scenario['size_gb'], scenario['coverage']
            )
            cost_result['scenario'] = scenario['name']
            results.append(cost_result)
            
            print(f"{platform:10s}: ¥{cost_result['total_cost']:,.0f} "
                  f"(¥{cost_result['cost_per_gb']:.0f}/GB)")
    
    return results

def create_cost_comparison_plots(results):
    """创建成本比较图表"""
    df = pd.DataFrame(results)
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('测序平台成本分析', fontsize=16, fontweight='bold')
    
    # 子图1：不同项目规模的总成本比较
    ax1 = axes[0, 0]
    pivot_total = df.pivot(index='scenario', columns='platform', values='total_cost')
    pivot_total.plot(kind='bar', ax=ax1, width=0.8)
    ax1.set_title('项目总成本比较')
    ax1.set_ylabel('总成本 (元)')
    ax1.set_xlabel('项目类型')
    ax1.legend(title='测序平台')
    ax1.tick_params(axis='x', rotation=45)
    
    # 子图2：每GB成本比较
    ax2 = axes[0, 1]
    pivot_per_gb = df.pivot(index='scenario', columns='platform', values='cost_per_gb')
    pivot_per_gb.plot(kind='bar', ax=ax2, width=0.8)
    ax2.set_title('每GB成本比较')
    ax2.set_ylabel('成本 (元/GB)')
    ax2.set_xlabel('项目类型')
    ax2.legend(title='测序平台')
    ax2.tick_params(axis='x', rotation=45)
    
    # 子图3：成本构成分析（以中型项目为例）
    ax3 = axes[1, 0]
    medium_project = df[df['scenario'] == '中型项目']
    
    cost_components = ['sequencing_cost', 'setup_cost', 'labor_cost', 
                      'equipment_cost', 'maintenance_cost']
    component_names = ['测序费用', '建库费用', '人工费用', '设备折旧', '维护费用']
    
    bottom = np.zeros(len(medium_project))
    colors = plt.cm.Set3(np.linspace(0, 1, len(cost_components)))
    
    for i, (component, name) in enumerate(zip(cost_components, component_names)):
        values = medium_project[component].values
        ax3.bar(medium_project['platform'], values, bottom=bottom, 
               label=name, color=colors[i])
        bottom += values
    
    ax3.set_title('成本构成分析（中型项目）')
    ax3.set_ylabel('成本 (元)')
    ax3.legend()
    
    # 子图4：成本效益分析（考虑准确率）
    ax4 = axes[1, 1]
    
    # 计算质量调整成本（成本/准确率）
    df['quality_adjusted_cost'] = df['cost_per_gb'] / df['accuracy']
    
    pivot_quality = df.pivot(index='scenario', columns='platform', values='quality_adjusted_cost')
    pivot_quality.plot(kind='bar', ax=ax4, width=0.8)
    ax4.set_title('质量调整成本比较')
    ax4.set_ylabel('质量调整成本 (元/GB)')
    ax4.set_xlabel('项目类型')
    ax4.legend(title='测序平台')
    ax4.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    
    # 保存图片
    os.makedirs('../results', exist_ok=True)
    plt.savefig('../results/cost_analysis.png', dpi=300, bbox_inches='tight')
    print("成本分析图已保存到: ../results/cost_analysis.png")
    
    return fig

def create_cost_sensitivity_analysis():
    """成本敏感性分析"""
    data_sizes = [10, 30, 50, 100, 200, 500, 1000]  # GB
    
    sensitivity_results = []
    
    for size in data_sizes:
        for platform in PLATFORM_COSTS.keys():
            result = calculate_project_cost(platform, size)
            sensitivity_results.append({
                'data_size': size,
                'platform': platform,
                'total_cost': result['total_cost'],
                'cost_per_gb': result['cost_per_gb']
            })
    
    df_sens = pd.DataFrame(sensitivity_results)
    
    # 创建敏感性分析图
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # 总成本随数据量变化
    for platform in PLATFORM_COSTS.keys():
        platform_data = df_sens[df_sens['platform'] == platform]
        ax1.plot(platform_data['data_size'], platform_data['total_cost'], 
                marker='o', label=platform, linewidth=2)
    
    ax1.set_xlabel('数据量 (GB)')
    ax1.set_ylabel('总成本 (元)')
    ax1.set_title('总成本随数据量变化')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    
    # 每GB成本随数据量变化
    for platform in PLATFORM_COSTS.keys():
        platform_data = df_sens[df_sens['platform'] == platform]
        ax2.plot(platform_data['data_size'], platform_data['cost_per_gb'], 
                marker='s', label=platform, linewidth=2)
    
    ax2.set_xlabel('数据量 (GB)')
    ax2.set_ylabel('每GB成本 (元)')
    ax2.set_title('每GB成本随数据量变化')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_xscale('log')
    
    plt.tight_layout()
    plt.savefig('../results/cost_sensitivity.png', dpi=300, bbox_inches='tight')
    print("敏感性分析图已保存到: ../results/cost_sensitivity.png")
    
    return df_sens

def save_cost_analysis_results(results):
    """保存成本分析结果"""
    df = pd.DataFrame(results)
    
    # 重新排列和重命名列
    df_output = df[['scenario', 'platform', 'data_size_gb', 'total_cost', 
                   'cost_per_gb', 'runs_needed', 'total_time_hours', 'accuracy']].copy()
    
    df_output.columns = ['项目类型', '测序平台', '数据量(GB)', '总成本(元)', 
                        '每GB成本(元)', '运行次数', '总时间(小时)', '准确率']
    
    # 数值格式化
    df_output['总成本(元)'] = df_output['总成本(元)'].round(0).astype(int)
    df_output['每GB成本(元)'] = df_output['每GB成本(元)'].round(0).astype(int)
    df_output['总时间(小时)'] = df_output['总时间(小时)'].round(1)
    
    os.makedirs('../results', exist_ok=True)
    df_output.to_csv('../results/cost_analysis_results.csv', index=False, encoding='utf-8-sig')
    print("成本分析结果已保存到: ../results/cost_analysis_results.csv")
    
    return df_output

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='测序成本分析工具')
    parser.add_argument('--project_size', type=float, default=100, 
                       help='项目数据量 (GB)')
    parser.add_argument('--coverage', type=int, default=30, 
                       help='测序覆盖度')
    
    args = parser.parse_args()
    
    print("测序平台成本分析")
    print("=" * 50)
    
    # 如果指定了特定项目参数，进行单项目分析
    if len(sys.argv) > 1:
        print(f"分析项目: {args.project_size} GB, {args.coverage}x 覆盖度")
        print("-" * 30)
        
        for platform in PLATFORM_COSTS.keys():
            result = calculate_project_cost(platform, args.project_size, args.coverage)
            print(f"{platform}:")
            print(f"  总成本: ¥{result['total_cost']:,.0f}")
            print(f"  每GB成本: ¥{result['cost_per_gb']:.0f}")
            print(f"  运行次数: {result['runs_needed']:.0f}")
            print(f"  总时间: {result['total_time_hours']:.1f} 小时")
            print()
    
    # 进行多场景分析
    print("多场景成本分析...")
    results = analyze_cost_scenarios()
    
    # 创建成本比较图表
    create_cost_comparison_plots(results)
    
    # 敏感性分析
    print("\n进行敏感性分析...")
    create_cost_sensitivity_analysis()
    
    # 保存结果
    save_cost_analysis_results(results)
    
    print("\n" + "=" * 50)
    print("成本分析完成！")
    print("结果文件:")
    print("  - 成本比较图: ../results/cost_analysis.png")
    print("  - 敏感性分析: ../results/cost_sensitivity.png")
    print("  - 详细结果: ../results/cost_analysis_results.csv")

if __name__ == "__main__":
    main()