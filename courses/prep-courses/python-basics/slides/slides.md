---
marp: true
theme: ngs-course
paginate: true
header: '高通量测序数据分析 | 王运生教授'
footer: 'wangys@hunau.edu.cn | 16教420室'
---

<!-- 
Python编程基础课程幻灯片
课程名称：高通量测序数据分析 - 预备课程
主讲教师：王运生教授
联系邮箱：wangys@hunau.edu.cn
办公室：16教420室
上课地点：105机房
-->

<!-- _class: title -->
# Python编程基础
## 高通量测序数据分析 - 预备课程

**主讲教师：** 王运生教授  
**联系邮箱：** wangys@hunau.edu.cn  
**办公室：** 16教420室  
**上课地点：** 105机房  

---

<!-- _class: toc -->
# 本次课程内容

1. Python简介与环境配置
2. 基础语法与数据类型
3. 数据结构详解
4. 控制流程与函数
5. 文件操作与异常处理
6. 生物信息学Python库

---

<!-- _class: content -->
# Python简介
## 学习目标
**学习目标：**
- 掌握Python基础编程技能
- 熟悉生物信息学常用数据结构
- 能够处理生物学数据文件
- 了解生物信息学Python生态

---

<!-- _class: content -->
# Python简介

## 什么是Python？

- **高级编程语言**：简洁易读的语法
- **解释型语言**：无需编译，直接运行
- **跨平台**：Windows、Linux、macOS都支持
- **丰富的库生态**：科学计算、数据分析、机器学习

---

<!-- _class: content -->
# Python环境配置

## 安装方式
## 为什么生物信息学选择Python？

- 语法简单，学习曲线平缓
- 强大的数据处理能力
- 丰富的生物信息学库（BioPython、pandas等）
- 活跃的科学计算社区

---

<!-- _class: multi-column -->
# Python在生物信息学中的应用

<div class="columns">
<div class="column">

## 数据处理
- **序列分析**：DNA/RNA/蛋白质序列处理
- **文件格式转换**：FASTA、FASTQ、VCF等
- **数据清洗**：质量控制、过滤
- **统计分析**：描述性统计、假设检验

</div>

---

<div class="column">

## 可视化分析
- **图表绘制**：matplotlib、seaborn
- **交互式可视化**：plotly、bokeh
- **基因组浏览**：pyGenomeTracks
- **网络分析**：networkx

</div>
</div>

---

<!-- _class: content -->
# Python环境配置

## 安装方式

### 1. Anaconda（推荐）
```bash
# 下载并安装Anaconda
wget https://repo.anaconda.com/archive/Anaconda3-latest-Linux-x86_64.sh
bash Anaconda3-latest-Linux-x86_64.sh

# 创建生物信息学环境
conda create -n bioinfo python=3.9
conda activate bioinfo
```

---

<!-- _class: code -->
# Python环境配置
### 2. 系统Python + pip
```bash
# Ubuntu/Debian
sudo apt install python3 python3-pip

# 安装包管理工具
pip install --upgrade pip
```

---

<!-- _class: code -->
# 第一个Python程序

```python
# hello_world.py
print("Hello, 生物信息学世界！")

# 运行方式
# 1. 交互式解释器
python3
>>> print("Hello World")

# 2. 脚本文件
python3 hello_world.py

# 3. Jupyter Notebook（推荐用于数据分析）
jupyter notebook
```
---

<!-- _class: content -->
# Python基础语法
## Python解释器特点
- **交互式**：适合测试和学习
- **脚本模式**：适合编写程序
- **Jupyter**：适合数据分析和可视化

---

<!-- _class: content -->
# 基础语法规则

## 代码风格特点

- **缩进敏感**：使用4个空格表示代码块
- **大小写敏感**：`Variable`和`variable`是不同的
- **注释**：`#`单行注释，`"""`多行注释

---

<!-- _class: code -->
# 基础语法示例
```python
# 这是单行注释
"""
这是多行注释
可以写多行内容
"""

# 缩进示例
if True:
    print("这行有4个空格缩进")
    if True:
        print("这行有8个空格缩进")
```

---

<!-- _class: multi-column -->
# 数据类型

<div class="columns">
<div class="column">

## 基本数据类型
```python
# 整数
age = 25
chromosome = 23

# 浮点数
gc_content = 0.42
temperature = 37.5

# 字符串
gene_name = "BRCA1"
sequence = "ATCGATCG"

# 布尔值
is_coding = True
has_mutation = False
```

</div>

---

<div class="column">

## 类型检查和转换
```python
# 查看类型
type(42)          # <class 'int'>
type(3.14)        # <class 'float'>
type("DNA")       # <class 'str'>

# 类型转换
int("123")        # 123
float("3.14")     # 3.14
str(42)           # "42"
bool(1)           # True
```

</div>
</div>

---

<!-- _class: content -->
# 字符串操作

## 字符串基础操作

```python
# 字符串定义
dna_seq = "ATCGATCGATCG"
gene_id = 'ENSG00000139618'

# 字符串拼接
full_id = gene_id + "_" + "BRCA2"
formatted = f"基因ID: {gene_id}, 长度: {len(dna_seq)}"

# 字符串方法
dna_seq.upper()           # 转大写
dna_seq.lower()           # 转小写
dna_seq.count("AT")       # 计数子串
dna_seq.replace("A", "T") # 替换
dna_seq.split("CG")       # 分割
```

---

<!-- _class: code -->
# 字符串格式化

```python
# 生物信息学中的字符串格式化
gene_name = "BRCA1"
chromosome = 17
start_pos = 43044295
end_pos = 43125483

# f-string格式化（Python 3.6+，推荐）
location = f"{gene_name}位于{chromosome}号染色体{start_pos}-{end_pos}"

# format方法
location = "{}位于{}号染色体{}-{}".format(gene_name, chromosome, start_pos, end_pos)

# 百分号格式化
location = "%s位于%d号染色体%d-%d" % (gene_name, chromosome, start_pos, end_pos)

print(location)
# 输出：BRCA1位于17号染色体43044295-43125483
```

---

<!-- _class: content -->
# 列表（List）

## 列表基础操作

```python
# 创建列表
genes = ["BRCA1", "BRCA2", "TP53", "EGFR"]
scores = [0.95, 0.87, 0.92, 0.78]
mixed = ["gene", 123, True, 3.14]

# 访问元素
first_gene = genes[0]        # "BRCA1"
last_gene = genes[-1]        # "EGFR"
subset = genes[1:3]          # ["BRCA2", "TP53"]

# 修改列表
genes.append("MYC")          # 添加元素
genes.insert(1, "KRAS")      # 插入元素
genes.remove("TP53")         # 删除元素
popped = genes.pop()         # 弹出最后一个元素
```

---

<!-- _class: code -->
# 列表高级操作

```python
# 列表推导式（List Comprehension）
sequences = ["ATCG", "GCTA", "TTAA", "CCGG"]

# 传统方法
gc_contents = []
for seq in sequences:
    gc_count = seq.count("G") + seq.count("C")
    gc_content = gc_count / len(seq)
    gc_contents.append(gc_content)

# 列表推导式（更Pythonic）
gc_contents = [(seq.count("G") + seq.count("C")) / len(seq) for seq in sequences]

# 条件过滤
high_gc = [seq for seq in sequences if (seq.count("G") + seq.count("C")) / len(seq) > 0.5]

print(gc_contents)  # [0.5, 0.5, 0.0, 1.0]
print(high_gc)      # ['ATCG', 'GCTA', 'CCGG']
```

---

<!-- _class: content -->
# 字典（Dictionary）

## 字典基础操作

```python
# 创建字典
gene_info = {
    "name": "BRCA1",
    "chromosome": 17,
    "start": 43044295,
    "end": 43125483,
    "strand": "+",
    "type": "protein_coding"
}

# 访问和修改
gene_name = gene_info["name"]           # 获取值
gene_info["length"] = 81188             # 添加键值对
gene_info["chromosome"] = "chr17"       # 修改值

# 安全访问
strand = gene_info.get("strand", "unknown")  # 如果键不存在返回默认值
```

---

<!-- _class: code -->
# 字典高级操作

```python
# 遍历字典
gene_expression = {
    "BRCA1": 5.2,
    "BRCA2": 3.8,
    "TP53": 7.1,
    "EGFR": 4.5
}

# 遍历键
for gene in gene_expression:
    print(f"基因: {gene}")

# 遍历键值对
for gene, expression in gene_expression.items():
    print(f"{gene}: {expression}")

# 字典推导式
high_expression = {gene: expr for gene, expr in gene_expression.items() if expr > 5.0}
print(high_expression)  # {'BRCA1': 5.2, 'TP53': 7.1}
```

---

<!-- _class: content -->
# 元组（Tuple）和集合（Set）

## 元组 - 不可变序列

```python
# 创建元组
coordinates = (17, 43044295, 43125483)  # 染色体位置
gene_info = ("BRCA1", "tumor_suppressor", 17)

# 元组解包
chromosome, start, end = coordinates
name, function, chr_num = gene_info

# 命名元组（更清晰的数据结构）
from collections import namedtuple
Gene = namedtuple('Gene', ['name', 'chromosome', 'start', 'end'])
brca1 = Gene("BRCA1", 17, 43044295, 43125483)
print(brca1.name)  # BRCA1
```
---

<!-- _class: code -->
## 集合 - 唯一元素集合

```python
# 创建集合
nucleotides = {"A", "T", "C", "G"}
amino_acids = set("ACDEFGHIKLMNPQRSTVWY")

# 集合操作
purines = {"A", "G"}
pyrimidines = {"C", "T"}
all_bases = purines | pyrimidines      # 并集
common = purines & pyrimidines         # 交集（空集）
```

---

<!-- _class: content -->
# 条件语句

## if-elif-else结构

```python
# 基因表达水平分类
expression_level = 6.5

if expression_level > 8.0:
    category = "高表达"
elif expression_level > 4.0:
    category = "中等表达"
elif expression_level > 1.0:
    category = "低表达"
else:
    category = "几乎不表达"

print(f"表达水平: {expression_level}, 分类: {category}")

# 序列质量评估
def assess_sequence_quality(gc_content, length):
    if 0.4 <= gc_content <= 0.6 and length >= 100:
        return "高质量序列"
    elif 0.3 <= gc_content <= 0.7 and length >= 50:
        return "中等质量序列"
    else:
        return "低质量序列"
```

---

<!-- _class: content -->
# 循环结构

## for循环

```python
# 遍历序列
dna_sequence = "ATCGATCGATCG"
for nucleotide in dna_sequence:
    print(f"核苷酸: {nucleotide}")

# 遍历列表
genes = ["BRCA1", "BRCA2", "TP53"]
for i, gene in enumerate(genes):
    print(f"第{i+1}个基因: {gene}")

# 遍历字典
expression_data = {"BRCA1": 5.2, "TP53": 7.1}
for gene, level in expression_data.items():
    print(f"{gene}的表达水平: {level}")
```
---

<!-- _class: code -->
## while循环

```python
# 序列搜索
sequence = "ATCGATCGATCGAAATCG"
pattern = "AAA"
position = 0

while position < len(sequence):
    if sequence[position:position+3] == pattern:
        print(f"在位置{position}找到模式{pattern}")
        break
    position += 1
```

---

<!-- _class: content -->
# 函数定义

## 基础函数

```python
def calculate_gc_content(sequence):
    """计算DNA序列的GC含量"""
    gc_count = sequence.count('G') + sequence.count('C')
    total_count = len(sequence)
    return gc_count / total_count if total_count > 0 else 0

# 使用函数
dna = "ATCGATCGATCG"
gc_percent = calculate_gc_content(dna)
print(f"GC含量: {gc_percent:.2%}")

def reverse_complement(sequence):
    """计算DNA序列的反向互补序列"""
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(complement.get(base, base) for base in sequence[::-1])

# 测试
original = "ATCGATCG"
rev_comp = reverse_complement(original)
print(f"原序列: {original}")
print(f"反向互补: {rev_comp}")
```

---

<!-- _class: code -->
# 函数高级特性

```python
# 默认参数
def translate_dna(sequence, genetic_code="standard"):
    """翻译DNA序列为蛋白质序列"""
    codon_table = {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        # ... 更多密码子
        "TAA": "*", "TAG": "*", "TGA": "*"  # 终止密码子
    }
    
    protein = ""
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        amino_acid = codon_table.get(codon, "X")  # X表示未知氨基酸
        protein += amino_acid
        if amino_acid == "*":  # 遇到终止密码子停止翻译
            break
    
    return protein

# 可变参数
def calculate_statistics(*values):
    """计算多个数值的统计信息"""
    if not values:
        return None
    
    return {
        "count": len(values),
        "sum": sum(values),
        "mean": sum(values) / len(values),
        "min": min(values),
        "max": max(values)
    }

# 使用示例
stats = calculate_statistics(5.2, 3.8, 7.1, 4.5, 6.3)
print(stats)
```

---

<!-- _class: content -->
# 文件操作

## 读取文件

```python
# 基本文件读取
with open("sequences.txt", "r") as file:
    content = file.read()
    print(content)

# 逐行读取
with open("gene_list.txt", "r") as file:
    for line in file:
        gene_name = line.strip()  # 去除换行符
        print(f"处理基因: {gene_name}")

# 读取FASTA文件示例
def read_fasta(filename):
    """简单的FASTA文件读取器"""
    sequences = {}
    current_id = None
    
    with open(filename, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                current_id = line[1:]  # 去掉>符号
                sequences[current_id] = ""
            elif current_id:
                sequences[current_id] += line
    
    return sequences
```

---

<!-- _class: code -->
# 文件写入和CSV处理

```python
# 写入文件
results = [
    ("BRCA1", 5.2, "高表达"),
    ("BRCA2", 3.8, "中等表达"),
    ("TP53", 7.1, "高表达")
]

with open("expression_results.txt", "w") as file:
    file.write("基因名\t表达水平\t分类\n")
    for gene, level, category in results:
        file.write(f"{gene}\t{level}\t{category}\n")

# CSV文件处理
import csv

# 读取CSV
def read_expression_data(filename):
    """读取基因表达数据CSV文件"""
    data = {}
    with open(filename, "r") as file:
        reader = csv.DictReader(file)
        for row in reader:
            gene = row["gene_name"]
            expression = float(row["expression_level"])
            data[gene] = expression
    return data

# 写入CSV
def write_results_csv(data, filename):
    """将结果写入CSV文件"""
    with open(filename, "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["基因名", "表达水平", "分类"])
        for gene, level, category in data:
            writer.writerow([gene, level, category])
```

---

<!-- _class: content -->
# 异常处理

## try-except结构

```python
def safe_divide(a, b):
    """安全的除法运算"""
    try:
        result = a / b
        return result
    except ZeroDivisionError:
        print("错误：除数不能为零")
        return None
    except TypeError:
        print("错误：参数必须是数字")
        return None

def read_sequence_file(filename):
    """安全地读取序列文件"""
    try:
        with open(filename, "r") as file:
            sequences = []
            for line in file:
                if not line.startswith(">"):
                    sequences.append(line.strip())
            return "".join(sequences)
    except FileNotFoundError:
        print(f"错误：文件 {filename} 不存在")
        return None
    except PermissionError:
        print(f"错误：没有权限读取文件 {filename}")
        return None
    except Exception as e:
        print(f"未知错误：{e}")
        return None
```

---

<!-- _class: content -->
# 生物信息学Python库 - NumPy

## NumPy基础

```python
import numpy as np

# 创建数组
expression_array = np.array([5.2, 3.8, 7.1, 4.5, 6.3])
quality_scores = np.array([[30, 35, 40], [25, 30, 35], [35, 40, 45]])

# 基本统计
mean_expression = np.mean(expression_array)
std_expression = np.std(expression_array)
max_quality = np.max(quality_scores)

print(f"平均表达水平: {mean_expression:.2f}")
print(f"表达水平标准差: {std_expression:.2f}")
print(f"最高质量分数: {max_quality}")

# 数组操作
normalized = (expression_array - np.mean(expression_array)) / np.std(expression_array)
high_expression = expression_array[expression_array > 5.0]

print(f"标准化后的表达水平: {normalized}")
print(f"高表达基因的表达水平: {high_expression}")
```

---

<!-- _class: content -->
# 生物信息学Python库 - Pandas

## Pandas基础

```python
import pandas as pd

# 创建DataFrame
gene_data = pd.DataFrame({
    'gene_name': ['BRCA1', 'BRCA2', 'TP53', 'EGFR', 'MYC'],
    'chromosome': [17, 13, 17, 7, 8],
    'expression': [5.2, 3.8, 7.1, 4.5, 6.3],
    'mutation_status': ['mutated', 'wild_type', 'mutated', 'wild_type', 'mutated']
})

# 基本操作
print(gene_data.head())
print(gene_data.describe())
print(gene_data['expression'].mean())

# 数据筛选
high_expression_genes = gene_data[gene_data['expression'] > 5.0]
mutated_genes = gene_data[gene_data['mutation_status'] == 'mutated']

# 分组统计
mutation_stats = gene_data.groupby('mutation_status')['expression'].agg(['mean', 'std', 'count'])
print(mutation_stats)
```

---

<!-- _class: code -->
# 生物信息学Python库 - BioPython

```python
# 需要先安装：pip install biopython
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC

# 序列对象
dna_seq = Seq("ATCGATCGATCG")
print(f"DNA序列: {dna_seq}")
print(f"反向互补: {dna_seq.reverse_complement()}")
print(f"转录为RNA: {dna_seq.transcribe()}")
print(f"翻译为蛋白质: {dna_seq.translate()}")

# GC含量计算
gc_content = GC(dna_seq)
print(f"GC含量: {gc_content:.2f}%")

# 读取FASTA文件
def analyze_fasta_file(filename):
    """分析FASTA文件中的序列"""
    sequences = []
    for record in SeqIO.parse(filename, "fasta"):
        seq_info = {
            'id': record.id,
            'description': record.description,
            'length': len(record.seq),
            'gc_content': GC(record.seq)
        }
        sequences.append(seq_info)
    return sequences

# 使用示例（需要有FASTA文件）
# sequences = analyze_fasta_file("sequences.fasta")
# for seq in sequences:
#     print(f"序列ID: {seq['id']}, 长度: {seq['length']}, GC含量: {seq['gc_content']:.2f}%")
```

---

<!-- _class: content -->
# 数据可视化 - Matplotlib

## 基础绘图

```python
import matplotlib.pyplot as plt
import numpy as np

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False    # 用来正常显示负号

# 基因表达水平柱状图
genes = ['BRCA1', 'BRCA2', 'TP53', 'EGFR', 'MYC']
expression = [5.2, 3.8, 7.1, 4.5, 6.3]

plt.figure(figsize=(10, 6))
plt.bar(genes, expression, color=['red', 'blue', 'green', 'orange', 'purple'])
plt.title('基因表达水平')
plt.xlabel('基因名称')
plt.ylabel('表达水平 (log2 FPKM)')
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

# 质量分数分布直方图
quality_scores = np.random.normal(30, 5, 1000)  # 模拟质量分数
plt.figure(figsize=(8, 6))
plt.hist(quality_scores, bins=30, alpha=0.7, color='skyblue')
plt.title('测序质量分数分布')
plt.xlabel('质量分数')
plt.ylabel('频次')
plt.show()
```

---

<!-- _class: content -->
# 实际应用案例

## 序列分析工具

```python
class SequenceAnalyzer:
    """DNA序列分析工具类"""
    
    def __init__(self, sequence):
        self.sequence = sequence.upper()
    
    def gc_content(self):
        """计算GC含量"""
        gc_count = self.sequence.count('G') + self.sequence.count('C')
        return gc_count / len(self.sequence) if len(self.sequence) > 0 else 0
    
    def find_orfs(self, min_length=100):
        """寻找开放阅读框（ORF）"""
        start_codon = "ATG"
        stop_codons = ["TAA", "TAG", "TGA"]
        orfs = []
        
        for frame in range(3):  # 三个阅读框
            for i in range(frame, len(self.sequence) - 2, 3):
                codon = self.sequence[i:i+3]
                if codon == start_codon:
                    # 寻找终止密码子
                    for j in range(i + 3, len(self.sequence) - 2, 3):
                        stop_codon = self.sequence[j:j+3]
                        if stop_codon in stop_codons:
                            orf_length = j - i + 3
                            if orf_length >= min_length:
                                orfs.append({
                                    'start': i,
                                    'end': j + 3,
                                    'length': orf_length,
                                    'frame': frame + 1,
                                    'sequence': self.sequence[i:j+3]
                                })
                            break
        return orfs

# 使用示例
sequence = "ATGAAATTTCCCGGGAAATAG"
analyzer = SequenceAnalyzer(sequence)
print(f"GC含量: {analyzer.gc_content():.2%}")
orfs = analyzer.find_orfs(min_length=10)
for orf in orfs:
    print(f"ORF: 位置{orf['start']}-{orf['end']}, 长度{orf['length']}, 阅读框{orf['frame']}")
```

---

<!-- _class: content -->
# 调试和测试

## 调试技巧

```python
# 使用print调试
def debug_function(data):
    print(f"输入数据: {data}")  # 调试信息
    result = process_data(data)
    print(f"处理结果: {result}")  # 调试信息
    return result

# 使用断言
def calculate_gc_content(sequence):
    assert isinstance(sequence, str), "序列必须是字符串"
    assert len(sequence) > 0, "序列不能为空"
    assert all(base in "ATCG" for base in sequence), "序列只能包含ATCG"
    
    gc_count = sequence.count('G') + sequence.count('C')
    return gc_count / len(sequence)

# 简单的单元测试
def test_gc_content():
    """测试GC含量计算函数"""
    # 测试用例1：已知结果
    assert calculate_gc_content("ATCG") == 0.5
    # 测试用例2：全是GC
    assert calculate_gc_content("GCGC") == 1.0
    # 测试用例3：没有GC
    assert calculate_gc_content("ATAT") == 0.0
    print("所有测试通过！")

test_gc_content()
```

---

<!-- _class: content -->
# 代码规范和最佳实践

## PEP 8 编码规范

```python
# 好的命名规范
gene_expression_data = {}  # 变量名用下划线
MAX_SEQUENCE_LENGTH = 1000  # 常量用大写

def calculate_gc_content(dna_sequence):  # 函数名用下划线
    """
    计算DNA序列的GC含量
    
    Args:
        dna_sequence (str): DNA序列字符串
        
    Returns:
        float: GC含量百分比
    """
    # 函数体
    pass

class SequenceProcessor:  # 类名用驼峰命名
    """序列处理器类"""
    
    def __init__(self, sequence):
        self.sequence = sequence
    
    def process(self):
        """处理序列"""
        pass

# 代码组织
import os           # 标准库
import sys

import numpy as np  # 第三方库
import pandas as pd

from .utils import helper_function  # 本地模块
```

---

<!-- _class: summary -->
# 本次课程总结

## 主要内容回顾
- Python基础语法和数据类型
- 数据结构：列表、字典、元组、集合
- 控制流程：条件语句和循环
- 函数定义和文件操作
- 生物信息学Python库介绍

---

## 重要概念
- **Pythonic编程**：简洁、可读的代码风格
- **数据结构选择**：根据需求选择合适的数据结构
- **异常处理**：编写健壮的程序
- **库的使用**：站在巨人的肩膀上

## 下次课程预告
- R语言基础
- 统计分析和数据可视化
- R在生物信息学中的应用

**作业/练习：**
- 完成Python基础编程练习
- 编写简单的序列分析脚本

---

<!-- _class: end -->
# 谢谢大家！

**有问题请联系：**
- 邮箱：wangys@hunau.edu.cn
- 办公室：16教420室

**推荐学习资源：**
- Python官方教程：https://docs.python.org/3/tutorial/
- BioPython教程：https://biopython.org/wiki/Documentation
- 《Python生物信息学编程》