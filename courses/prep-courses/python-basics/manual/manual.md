# Python编程基础实践操作手册

## 课程信息
- **课程名称**：高通量测序数据分析 - Python编程基础
- **主讲教师**：王运生教授
- **联系邮箱**：wangys@hunau.edu.cn
- **办公室**：16教420室
- **上课地点**：105机房

## 实验目标
通过本次实践课程，学生将能够：
1. 配置Python开发环境并熟悉基本操作
2. 掌握Python基础语法和数据结构
3. 编写简单的生物信息学数据处理脚本
4. 使用Python处理生物学数据文件
5. 了解并使用生物信息学相关Python库

## 环境要求

### 软件环境
- **Python版本**：3.8或更高版本
- **推荐发行版**：Anaconda或Miniconda
- **代码编辑器**：VS Code、PyCharm或Jupyter Notebook
- **必需Python包**：
  ```bash
  pip install numpy pandas matplotlib seaborn biopython jupyter
  ```

### 硬件要求
- **内存**：至少4GB RAM
- **存储空间**：至少2GB可用空间
- **网络**：用于下载包和数据

### 数据准备
本次实验将使用以下示例数据文件：
- `sample_sequences.fasta`：示例DNA序列文件
- `gene_expression.csv`：基因表达数据
- `quality_scores.txt`：测序质量分数数据

## 操作步骤

### 步骤1：Python环境设置和验证

#### 1.1 检查Python安装
```bash
# 检查Python版本
python3 --version
# 或者
python --version

# 检查pip版本
pip --version
```

**预期结果**：显示Python 3.8+版本信息

#### 1.2 创建虚拟环境（推荐）
```bash
# 创建虚拟环境
python3 -m venv bioinfo_env

# 激活虚拟环境
# Linux/Mac:
source bioinfo_env/bin/activate
# Windows:
bioinfo_env\Scripts\activate

# 验证环境激活
which python
```

#### 1.3 安装必需的包
```bash
# 升级pip
pip install --upgrade pip

# 安装基础科学计算包
pip install numpy pandas matplotlib seaborn

# 安装生物信息学包
pip install biopython

# 安装Jupyter Notebook
pip install jupyter

# 验证安装
python -c "import numpy, pandas, matplotlib, Bio; print('所有包安装成功！')"
```

**运行脚本**：`scripts/setup.py`

### 步骤2：Python基础语法练习

#### 2.1 启动Python交互环境
```bash
# 启动Python解释器
python3

# 或启动Jupyter Notebook
jupyter notebook
```

#### 2.2 基础数据类型练习
在Python解释器中执行以下代码：

```python
# 基础数据类型
gene_name = "BRCA1"
chromosome = 17
gc_content = 0.42
is_tumor_suppressor = True

print(f"基因名: {gene_name}")
print(f"染色体: {chromosome}")
print(f"GC含量: {gc_content:.2%}")
print(f"是否为抑癌基因: {is_tumor_suppressor}")

# 类型检查
print(f"gene_name的类型: {type(gene_name)}")
print(f"chromosome的类型: {type(chromosome)}")
```

#### 2.3 字符串操作练习
```python
# DNA序列操作
dna_sequence = "ATCGATCGATCG"

# 基本操作
print(f"序列长度: {len(dna_sequence)}")
print(f"大写序列: {dna_sequence.upper()}")
print(f"A的数量: {dna_sequence.count('A')}")

# 序列分析
def analyze_sequence(seq):
    """分析DNA序列的基本信息"""
    seq = seq.upper()
    length = len(seq)
    a_count = seq.count('A')
    t_count = seq.count('T')
    c_count = seq.count('C')
    g_count = seq.count('G')
    gc_content = (g_count + c_count) / length if length > 0 else 0
    
    return {
        'length': length,
        'A': a_count,
        'T': t_count,
        'C': c_count,
        'G': g_count,
        'GC_content': gc_content
    }

# 测试函数
result = analyze_sequence(dna_sequence)
print("序列分析结果:")
for key, value in result.items():
    if key == 'GC_content':
        print(f"{key}: {value:.2%}")
    else:
        print(f"{key}: {value}")
```

**运行脚本**：`scripts/basic_exercises.py`

### 步骤3：数据结构操作

#### 3.1 列表操作练习
```python
# 基因列表操作
genes = ["BRCA1", "BRCA2", "TP53", "EGFR", "MYC"]

# 基本操作
print(f"基因列表: {genes}")
print(f"第一个基因: {genes[0]}")
print(f"最后一个基因: {genes[-1]}")
print(f"前三个基因: {genes[:3]}")

# 添加和删除
genes.append("KRAS")
genes.insert(2, "PIK3CA")
print(f"添加后的列表: {genes}")

# 列表推导式 - 筛选包含'BR'的基因
br_genes = [gene for gene in genes if 'BR' in gene]
print(f"包含'BR'的基因: {br_genes}")

# 计算每个基因名的长度
gene_lengths = [len(gene) for gene in genes]
print(f"基因名长度: {gene_lengths}")
```

#### 3.2 字典操作练习
```python
# 基因表达数据
gene_expression = {
    "BRCA1": 5.2,
    "BRCA2": 3.8,
    "TP53": 7.1,
    "EGFR": 4.5,
    "MYC": 6.3
}

# 基本操作
print(f"BRCA1的表达水平: {gene_expression['BRCA1']}")
print(f"所有基因: {list(gene_expression.keys())}")
print(f"所有表达水平: {list(gene_expression.values())}")

# 添加新数据
gene_expression["KRAS"] = 4.8
gene_expression["PIK3CA"] = 5.5

# 筛选高表达基因（>5.0）
high_expression = {gene: expr for gene, expr in gene_expression.items() if expr > 5.0}
print(f"高表达基因: {high_expression}")

# 计算统计信息
expressions = list(gene_expression.values())
mean_expr = sum(expressions) / len(expressions)
max_expr = max(expressions)
min_expr = min(expressions)

print(f"平均表达水平: {mean_expr:.2f}")
print(f"最高表达水平: {max_expr}")
print(f"最低表达水平: {min_expr}")
```

**运行脚本**：`scripts/data_structures.py`

### 步骤4：控制流程和函数

#### 4.1 条件语句练习
```python
def classify_expression(expression_level):
    """根据表达水平对基因进行分类"""
    if expression_level > 6.0:
        return "高表达"
    elif expression_level > 3.0:
        return "中等表达"
    elif expression_level > 1.0:
        return "低表达"
    else:
        return "几乎不表达"

# 测试分类函数
test_levels = [7.5, 4.2, 2.1, 0.5]
for level in test_levels:
    category = classify_expression(level)
    print(f"表达水平 {level}: {category}")
```

#### 4.2 循环练习
```python
# 批量处理序列
sequences = [
    "ATCGATCG",
    "GCTAGCTA", 
    "TTAACCGG",
    "CGATCGAT"
]

print("序列分析结果:")
for i, seq in enumerate(sequences, 1):
    analysis = analyze_sequence(seq)
    print(f"序列 {i}: 长度={analysis['length']}, GC含量={analysis['GC_content']:.2%}")

# 寻找特定模式
pattern = "CG"
print(f"\n寻找模式 '{pattern}':")
for i, seq in enumerate(sequences, 1):
    count = seq.count(pattern)
    if count > 0:
        print(f"序列 {i}: 找到 {count} 个 '{pattern}' 模式")
```

#### 4.3 函数定义练习
```python
def reverse_complement(sequence):
    """计算DNA序列的反向互补序列"""
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    # 反向并互补
    reverse_seq = sequence[::-1]  # 反向
    complement_seq = ''.join(complement_map.get(base, base) for base in reverse_seq)
    
    return complement_seq

def translate_dna(sequence, start_pos=0):
    """将DNA序列翻译为蛋白质序列"""
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
    
    protein = ""
    for i in range(start_pos, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        amino_acid = codon_table.get(codon, 'X')  # X表示未知氨基酸
        protein += amino_acid
        if amino_acid == '*':  # 遇到终止密码子
            break
    
    return protein

# 测试函数
test_sequence = "ATGAAATTTCCCGGGAAATAG"
print(f"原序列: {test_sequence}")
print(f"反向互补: {reverse_complement(test_sequence)}")
print(f"翻译结果: {translate_dna(test_sequence)}")
```

### 步骤5：文件操作

#### 5.1 创建示例数据文件
首先创建一些示例数据文件用于练习：

```python
# 创建示例FASTA文件
fasta_content = """>gene1|BRCA1
ATGAAATTTCCCGGGAAATAGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>gene2|TP53
ATGCCCGGGAAATTTAAATTTCCCGGGAAATAGATCGATCGATCGATCGATC
GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
>gene3|EGFR
ATGAAATTTCCCGGGAAATAGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"""

with open("sample_sequences.fasta", "w") as f:
    f.write(fasta_content)

# 创建基因表达CSV文件
import csv

expression_data = [
    ["gene_name", "sample1", "sample2", "sample3", "sample4"],
    ["BRCA1", "5.2", "4.8", "5.5", "4.9"],
    ["BRCA2", "3.8", "4.1", "3.5", "3.9"],
    ["TP53", "7.1", "6.8", "7.3", "7.0"],
    ["EGFR", "4.5", "4.2", "4.8", "4.6"],
    ["MYC", "6.3", "6.0", "6.5", "6.1"]
]

with open("gene_expression.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerows(expression_data)

print("示例数据文件创建完成")
```

#### 5.2 读取和处理FASTA文件
```python
def read_fasta_simple(filename):
    """简单的FASTA文件读取器"""
    sequences = {}
    current_id = None
    
    try:
        with open(filename, "r") as file:
            for line in file:
                line = line.strip()
                if line.startswith(">"):
                    # 解析序列ID
                    current_id = line[1:].split("|")[0]  # 取>后面第一部分
                    sequences[current_id] = ""
                elif current_id:
                    sequences[current_id] += line
        
        return sequences
    except FileNotFoundError:
        print(f"错误：文件 {filename} 不存在")
        return {}

# 读取并分析FASTA文件
sequences = read_fasta_simple("sample_sequences.fasta")
print("FASTA文件分析结果:")
for seq_id, sequence in sequences.items():
    analysis = analyze_sequence(sequence)
    print(f"{seq_id}: 长度={analysis['length']}, GC含量={analysis['GC_content']:.2%}")
```

#### 5.3 处理CSV文件
```python
import csv

def read_expression_csv(filename):
    """读取基因表达CSV文件"""
    data = {}
    
    try:
        with open(filename, "r") as file:
            reader = csv.DictReader(file)
            for row in reader:
                gene_name = row["gene_name"]
                # 提取数值列（除了gene_name）
                expressions = []
                for key, value in row.items():
                    if key != "gene_name":
                        try:
                            expressions.append(float(value))
                        except ValueError:
                            continue
                data[gene_name] = expressions
        
        return data
    except FileNotFoundError:
        print(f"错误：文件 {filename} 不存在")
        return {}

# 读取并分析表达数据
expression_data = read_expression_csv("gene_expression.csv")
print("\n基因表达数据分析:")
for gene, values in expression_data.items():
    mean_expr = sum(values) / len(values)
    max_expr = max(values)
    min_expr = min(values)
    print(f"{gene}: 平均={mean_expr:.2f}, 最大={max_expr}, 最小={min_expr}")
```

**运行脚本**：`scripts/file_operations.py`

### 步骤6：使用生物信息学库

#### 6.1 NumPy数组操作
```python
import numpy as np

# 创建表达矩阵
genes = ["BRCA1", "BRCA2", "TP53", "EGFR", "MYC"]
samples = ["Sample1", "Sample2", "Sample3", "Sample4"]

# 模拟表达数据
np.random.seed(42)  # 设置随机种子以获得可重现结果
expression_matrix = np.random.normal(5.0, 1.5, (len(genes), len(samples)))
expression_matrix = np.abs(expression_matrix)  # 确保表达值为正

print("基因表达矩阵:")
print(f"{'基因名':<8}", end="")
for sample in samples:
    print(f"{sample:>10}", end="")
print()

for i, gene in enumerate(genes):
    print(f"{gene:<8}", end="")
    for j in range(len(samples)):
        print(f"{expression_matrix[i, j]:>10.2f}", end="")
    print()

# 计算统计信息
print(f"\n每个基因的平均表达水平:")
gene_means = np.mean(expression_matrix, axis=1)
for i, gene in enumerate(genes):
    print(f"{gene}: {gene_means[i]:.2f}")

print(f"\n每个样本的平均表达水平:")
sample_means = np.mean(expression_matrix, axis=0)
for i, sample in enumerate(samples):
    print(f"{sample}: {sample_means[i]:.2f}")
```

#### 6.2 Pandas数据分析
```python
import pandas as pd

# 创建DataFrame
df = pd.DataFrame(expression_matrix, 
                  index=genes, 
                  columns=samples)

print("使用Pandas分析基因表达数据:")
print(df)

print("\n基本统计信息:")
print(df.describe())

print("\n每个基因的统计信息:")
gene_stats = df.T.describe().T  # 转置以获得基因的统计信息
print(gene_stats)

# 筛选高表达基因（平均表达>5.0）
high_expression_genes = df[df.mean(axis=1) > 5.0]
print(f"\n高表达基因（平均>5.0）:")
print(high_expression_genes)

# 相关性分析
correlation_matrix = df.T.corr()  # 样本间相关性
print(f"\n样本间相关性:")
print(correlation_matrix)
```

#### 6.3 BioPython序列分析
```python
try:
    from Bio.Seq import Seq
    from Bio.SeqUtils import GC
    from Bio import SeqIO
    
    # 序列对象操作
    dna_seq = Seq("ATGAAATTTCCCGGGAAATAGATCGATCGATCG")
    
    print("BioPython序列分析:")
    print(f"DNA序列: {dna_seq}")
    print(f"长度: {len(dna_seq)}")
    print(f"GC含量: {GC(dna_seq):.2f}%")
    print(f"反向互补: {dna_seq.reverse_complement()}")
    print(f"转录为RNA: {dna_seq.transcribe()}")
    print(f"翻译为蛋白质: {dna_seq.translate()}")
    
    # 处理FASTA文件
    print(f"\n使用BioPython读取FASTA文件:")
    for record in SeqIO.parse("sample_sequences.fasta", "fasta"):
        print(f"ID: {record.id}")
        print(f"描述: {record.description}")
        print(f"长度: {len(record.seq)}")
        print(f"GC含量: {GC(record.seq):.2f}%")
        print(f"前50个碱基: {record.seq[:50]}")
        print("-" * 50)

except ImportError:
    print("BioPython未安装，跳过BioPython练习")
    print("可以使用以下命令安装：pip install biopython")
```

**运行脚本**：`scripts/bioinformatics_intro.py`

### 步骤7：数据可视化

#### 7.1 基础图表绘制
```python
import matplotlib.pyplot as plt
import numpy as np

# 设置中文字体（如果需要）
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 基因表达柱状图
genes = ["BRCA1", "BRCA2", "TP53", "EGFR", "MYC"]
expression = [5.2, 3.8, 7.1, 4.5, 6.3]

plt.figure(figsize=(10, 6))
bars = plt.bar(genes, expression, color=['red', 'blue', 'green', 'orange', 'purple'])
plt.title('Gene Expression Levels')
plt.xlabel('Gene Names')
plt.ylabel('Expression Level (log2 FPKM)')
plt.xticks(rotation=45)

# 添加数值标签
for bar, value in zip(bars, expression):
    plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1, 
             f'{value}', ha='center', va='bottom')

plt.tight_layout()
plt.savefig('gene_expression.png', dpi=300, bbox_inches='tight')
plt.show()

# GC含量分布直方图
np.random.seed(42)
gc_contents = np.random.normal(0.45, 0.1, 1000)
gc_contents = np.clip(gc_contents, 0, 1)  # 限制在0-1范围内

plt.figure(figsize=(8, 6))
plt.hist(gc_contents, bins=30, alpha=0.7, color='skyblue', edgecolor='black')
plt.title('GC Content Distribution')
plt.xlabel('GC Content')
plt.ylabel('Frequency')
plt.axvline(np.mean(gc_contents), color='red', linestyle='--', 
            label=f'Mean: {np.mean(gc_contents):.3f}')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('gc_distribution.png', dpi=300, bbox_inches='tight')
plt.show()
```

## 预期结果

### 步骤1预期结果
- Python环境成功配置
- 所有必需包安装完成
- 虚拟环境正常工作

### 步骤2预期结果
- 能够执行基本Python语句
- 正确显示变量类型和值
- 字符串操作功能正常

### 步骤3预期结果
- 列表和字典操作正确执行
- 数据筛选和统计计算准确
- 推导式语法使用正确

### 步骤4预期结果
- 条件语句正确分类表达水平
- 循环能够批量处理序列
- 自定义函数正常工作

### 步骤5预期结果
- 成功创建和读取文件
- FASTA和CSV文件解析正确
- 文件操作异常处理有效

### 步骤6预期结果
- NumPy数组操作和统计计算正确
- Pandas数据分析功能正常
- BioPython序列分析准确

### 步骤7预期结果
- 生成清晰的基因表达柱状图
- GC含量分布直方图显示正常
- 图片文件成功保存

## 故障排除

### 常见问题1：包安装失败
**症状**：pip install命令报错
**解决方案**：
```bash
# 升级pip
python -m pip install --upgrade pip

# 使用国内镜像源
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple/ package_name

# 检查网络连接
ping pypi.org
```

### 常见问题2：中文字符显示异常
**症状**：图表中中文显示为方块
**解决方案**：
```python
# 安装中文字体
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
```

### 常见问题3：文件路径错误
**症状**：FileNotFoundError
**解决方案**：
```python
import os
print("当前工作目录:", os.getcwd())
print("文件是否存在:", os.path.exists("filename.txt"))

# 使用绝对路径
import os.path
filepath = os.path.abspath("filename.txt")
```

### 常见问题4：内存不足
**症状**：处理大文件时程序崩溃
**解决方案**：
```python
# 分块读取大文件
def read_large_file(filename, chunk_size=1000):
    with open(filename, 'r') as f:
        while True:
            chunk = f.readlines(chunk_size)
            if not chunk:
                break
            yield chunk
```

## 扩展练习

### 练习1：序列分析工具
编写一个完整的DNA序列分析工具，包括：
- GC含量计算
- 开放阅读框（ORF）查找
- 反向互补序列生成
- 密码子使用频率分析

### 练习2：基因表达数据分析
使用提供的基因表达数据：
- 计算基因间相关性
- 识别差异表达基因
- 绘制热图和散点图
- 进行聚类分析

### 练习3：文件格式转换
编写脚本实现：
- FASTA到CSV格式转换
- 多个文件批量处理
- 数据质量检查和报告
- 错误处理和日志记录

### 练习4：生物信息学流水线
设计一个简单的分析流水线：
- 读取多个FASTA文件
- 批量计算序列统计信息
- 生成分析报告
- 创建可视化图表

## 参考资料

### 官方文档
- [Python官方教程](https://docs.python.org/3/tutorial/)
- [NumPy用户指南](https://numpy.org/doc/stable/user/)
- [Pandas用户指南](https://pandas.pydata.org/docs/user_guide/)
- [Matplotlib教程](https://matplotlib.org/stable/tutorials/index.html)
- [BioPython教程](https://biopython.org/wiki/Documentation)

### 推荐书籍
- 《Python生物信息学编程》
- 《利用Python进行数据分析》
- 《Python科学计算》

### 在线资源
- [Python.org](https://www.python.org/)
- [Anaconda文档](https://docs.anaconda.com/)
- [Jupyter Notebook文档](https://jupyter-notebook.readthedocs.io/)

## 课后作业

1. **基础练习**：完成所有步骤中的代码练习
2. **编程作业**：编写一个DNA序列分析脚本，要求包含至少5个分析功能
3. **数据分析**：使用提供的基因表达数据进行探索性数据分析
4. **可视化作业**：创建至少3种不同类型的生物信息学相关图表

**提交要求**：
- 所有Python脚本文件（.py格式）
- 分析结果和图表
- 简短的分析报告（500字以内）

**截止时间**：下次课前

---

**注意事项**：
- 保存好所有练习代码，后续课程会用到
- 遇到问题及时询问，不要积累疑问
- 多练习，熟能生巧
- 关注代码规范和注释习惯