#!/usr/bin/env python3

"""
MultiQC分析脚本
课程：高通量测序数据分析 - 测序数据质量控制与预处理
作者：王运生
日期：2025-01-20
用法：python run_multiqc.py [input_dir] [output_dir] [report_name]
"""

import os
import sys
import subprocess
import argparse
import json
from datetime import datetime
from pathlib import Path

def check_multiqc_installation():
    """检查MultiQC是否已安装"""
    try:
        result = subprocess.run(['multiqc', '--version'], 
                              capture_output=True, text=True, check=True)
        version = result.stdout.strip()
        print(f"✓ MultiQC已安装: {version}")
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("✗ MultiQC未安装或不在PATH中")
        print("请运行: pip install multiqc")
        return False

def find_fastqc_results(input_dir):
    """查找FastQC结果文件"""
    input_path = Path(input_dir)
    
    # 查找FastQC的ZIP文件
    zip_files = list(input_path.glob("*_fastqc.zip"))
    html_files = list(input_path.glob("*_fastqc.html"))
    
    print(f"在 {input_dir} 中找到:")
    print(f"  FastQC ZIP文件: {len(zip_files)} 个")
    print(f"  FastQC HTML文件: {len(html_files)} 个")
    
    if len(zip_files) == 0 and len(html_files) == 0:
        print("警告：未找到FastQC结果文件")
        print("请确保已运行FastQC分析")
        return False
    
    return True

def create_multiqc_config(output_dir, report_name):
    """创建MultiQC配置文件"""
    config = {
        "title": f"质量控制报告 - {report_name}",
        "subtitle": "高通量测序数据分析课程",
        "intro_text": "本报告由MultiQC生成，整合了FastQC等工具的分析结果。",
        "report_comment": f"生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "show_analysis_paths": False,
        "show_analysis_time": True,
        "custom_logo": None,
        "custom_logo_url": None,
        "custom_logo_title": "NGS数据分析课程",
        "report_header_info": [
            {"课程": "高通量测序数据分析"},
            {"讲师": "王运生教授"},
            {"邮箱": "wangys@hunau.edu.cn"},
            {"生成时间": datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
        ],
        "module_order": [
            "fastqc",
            "trimmomatic",
            "cutadapt"
        ],
        "fastqc": {
            "general_stats_target_coverage": [20, 40, 70, 80, 90, 95, 99]
        }
    }
    
    config_file = Path(output_dir) / "multiqc_config.yaml"
    
    # 将配置写入YAML文件
    import yaml
    try:
        with open(config_file, 'w', encoding='utf-8') as f:
            yaml.dump(config, f, default_flow_style=False, allow_unicode=True)
        print(f"✓ MultiQC配置文件已创建: {config_file}")
        return str(config_file)
    except ImportError:
        print("警告：PyYAML未安装，使用默认配置")
        return None

def run_multiqc(input_dir, output_dir, report_name, config_file=None):
    """运行MultiQC分析"""
    
    # 构建MultiQC命令
    cmd = [
        'multiqc',
        input_dir,
        '-o', output_dir,
        '-n', f'{report_name}.html',
        '--force',  # 覆盖已存在的报告
        '--interactive',  # 生成交互式图表
        '--export',  # 导出数据表格
        '--zip-data-dir',  # 压缩数据目录
    ]
    
    # 如果有配置文件，添加到命令中
    if config_file:
        cmd.extend(['-c', config_file])
    
    print("运行MultiQC命令:")
    print(" ".join(cmd))
    print("=" * 50)
    
    try:
        # 运行MultiQC
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        print("✓ MultiQC分析成功完成")
        print("\nMultiQC输出:")
        print(result.stdout)
        
        if result.stderr:
            print("\n警告信息:")
            print(result.stderr)
            
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"✗ MultiQC分析失败 (退出代码: {e.returncode})")
        print(f"错误信息: {e.stderr}")
        return False

def analyze_multiqc_results(output_dir, report_name):
    """分析MultiQC结果并生成摘要"""
    
    output_path = Path(output_dir)
    
    # 检查生成的文件
    html_report = output_path / f"{report_name}.html"
    data_dir = output_path / "multiqc_data"
    
    print("\n检查生成的文件:")
    
    if html_report.exists():
        size = html_report.stat().st_size / 1024  # KB
        print(f"✓ HTML报告: {html_report} ({size:.1f} KB)")
    else:
        print(f"✗ HTML报告未生成: {html_report}")
        return False
    
    if data_dir.exists():
        data_files = list(data_dir.glob("*.txt"))
        print(f"✓ 数据目录: {data_dir} ({len(data_files)} 个数据文件)")
        
        # 列出主要数据文件
        for data_file in data_files:
            if data_file.name.startswith('multiqc_'):
                print(f"  - {data_file.name}")
    else:
        print(f"✗ 数据目录未生成: {data_dir}")
    
    # 尝试解析general stats
    general_stats_file = data_dir / "multiqc_general_stats.txt"
    if general_stats_file.exists():
        print(f"\n解析样本统计信息:")
        try:
            with open(general_stats_file, 'r') as f:
                lines = f.readlines()
                if len(lines) > 1:  # 有数据行
                    header = lines[0].strip().split('\t')
                    print(f"  样本数量: {len(lines) - 1}")
                    print(f"  统计指标: {len(header) - 1}")  # 减去Sample列
                    
                    # 显示前几个样本的信息
                    for i, line in enumerate(lines[1:4], 1):  # 显示前3个样本
                        if line.strip():
                            sample_data = line.strip().split('\t')
                            sample_name = sample_data[0]
                            print(f"  样本{i}: {sample_name}")
        except Exception as e:
            print(f"  解析统计文件时出错: {e}")
    
    return True

def create_summary_report(output_dir, report_name, input_dir):
    """创建分析摘要报告"""
    
    summary_file = Path(output_dir) / "analysis_summary.txt"
    
    with open(summary_file, 'w', encoding='utf-8') as f:
        f.write("MultiQC分析摘要报告\n")
        f.write("=" * 30 + "\n\n")
        f.write(f"分析时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"输入目录: {input_dir}\n")
        f.write(f"输出目录: {output_dir}\n")
        f.write(f"报告名称: {report_name}\n\n")
        
        # 统计输入文件
        input_path = Path(input_dir)
        zip_files = list(input_path.glob("*_fastqc.zip"))
        html_files = list(input_path.glob("*_fastqc.html"))
        
        f.write("输入文件统计:\n")
        f.write(f"- FastQC ZIP文件: {len(zip_files)} 个\n")
        f.write(f"- FastQC HTML文件: {len(html_files)} 个\n\n")
        
        f.write("输入文件列表:\n")
        for zip_file in zip_files:
            f.write(f"- {zip_file.name}\n")
        
        # 输出文件信息
        output_path = Path(output_dir)
        html_report = output_path / f"{report_name}.html"
        data_dir = output_path / "multiqc_data"
        
        f.write(f"\n输出文件:\n")
        if html_report.exists():
            size = html_report.stat().st_size / 1024
            f.write(f"- HTML报告: {html_report.name} ({size:.1f} KB)\n")
        
        if data_dir.exists():
            data_files = list(data_dir.glob("*.txt"))
            f.write(f"- 数据文件: {len(data_files)} 个\n")
            for data_file in data_files:
                f.write(f"  - {data_file.name}\n")
        
        f.write(f"\n使用说明:\n")
        f.write(f"1. 在浏览器中打开 {html_report.name} 查看完整报告\n")
        f.write(f"2. multiqc_data/ 目录包含原始数据，可用于进一步分析\n")
        f.write(f"3. 重点关注样本间的质量差异和异常样本\n")
    
    print(f"✓ 分析摘要已保存: {summary_file}")

def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description="MultiQC分析脚本 - 整合FastQC等工具的质量控制报告",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:
  python run_multiqc.py                                    # 使用默认参数
  python run_multiqc.py results/fastqc_raw               # 指定输入目录
  python run_multiqc.py results/fastqc_raw results/multiqc raw_data_report  # 完整参数
        """
    )
    
    parser.add_argument('input_dir', nargs='?', default='results/fastqc_raw',
                       help='FastQC结果目录 (默认: results/fastqc_raw)')
    parser.add_argument('output_dir', nargs='?', default='results/multiqc',
                       help='MultiQC输出目录 (默认: results/multiqc)')
    parser.add_argument('report_name', nargs='?', default='multiqc_report',
                       help='报告名称 (默认: multiqc_report)')
    parser.add_argument('--config', help='MultiQC配置文件路径')
    parser.add_argument('--no-config', action='store_true', help='不创建配置文件')
    
    args = parser.parse_args()
    
    print("=" * 50)
    print("MultiQC分析脚本")
    print("=" * 50)
    print(f"输入目录: {args.input_dir}")
    print(f"输出目录: {args.output_dir}")
    print(f"报告名称: {args.report_name}")
    print("=" * 50)
    
    # 检查MultiQC安装
    if not check_multiqc_installation():
        sys.exit(1)
    
    # 检查输入目录
    if not os.path.exists(args.input_dir):
        print(f"错误：输入目录不存在: {args.input_dir}")
        sys.exit(1)
    
    # 查找FastQC结果
    if not find_fastqc_results(args.input_dir):
        print("请先运行FastQC分析生成结果文件")
        sys.exit(1)
    
    # 创建输出目录
    os.makedirs(args.output_dir, exist_ok=True)
    
    # 创建配置文件（如果需要）
    config_file = None
    if not args.no_config and not args.config:
        try:
            config_file = create_multiqc_config(args.output_dir, args.report_name)
        except Exception as e:
            print(f"警告：创建配置文件失败: {e}")
    elif args.config:
        config_file = args.config
    
    # 运行MultiQC
    success = run_multiqc(args.input_dir, args.output_dir, args.report_name, config_file)
    
    if success:
        # 分析结果
        analyze_multiqc_results(args.output_dir, args.report_name)
        
        # 创建摘要报告
        create_summary_report(args.output_dir, args.report_name, args.input_dir)
        
        print("\n" + "=" * 50)
        print("MultiQC分析完成！")
        print("=" * 50)
        print(f"HTML报告: {args.output_dir}/{args.report_name}.html")
        print(f"数据目录: {args.output_dir}/multiqc_data/")
        print(f"分析摘要: {args.output_dir}/analysis_summary.txt")
        print("\n建议下一步:")
        print("1. 在浏览器中打开HTML报告")
        print("2. 检查样本质量和异常情况")
        print("3. 根据结果制定数据清洗策略")
        print("=" * 50)
        
    else:
        print("\nMultiQC分析失败，请检查错误信息")
        sys.exit(1)

if __name__ == "__main__":
    main()