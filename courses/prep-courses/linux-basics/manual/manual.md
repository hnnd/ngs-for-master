# Linux系统基础实践操作手册

## 课程信息
- **课程名称**：高通量测序数据分析 - Linux系统基础
- **主讲教师**：王运生
- **联系邮箱**：wangys@hunau.edu.cn
- **办公室**：16教420室
- **上课地点**：105机房
- **课程时间**：4学时（理论1学时 + 实践3学时）

## 实验目标

### 主要目标
- 掌握Linux系统的基本操作和命令行使用
- 熟练使用文件和目录管理命令
- 学会使用文本处理工具处理生物信息学数据
- 理解Linux系统在生物信息学中的应用

### 预期成果
- 能够熟练使用命令行进行文件操作
- 掌握基本的文本处理和数据筛选技能
- 能够编写简单的Bash脚本自动化任务
- 理解Linux文件系统和权限管理

## 环境要求

### 软件环境
| 软件名称 | 版本要求 | 安装方式 | 说明 |
|---------|---------|---------|------|
| Ubuntu Linux | 18.04+ | 系统安装 | 操作系统 |
| Bash | 4.0+ | 系统自带 | 命令行解释器 |
| nano/vim | 任意版本 | `sudo apt install nano vim` | 文本编辑器 |
| tree | 任意版本 | `sudo apt install tree` | 目录树显示 |
| htop | 任意版本 | `sudo apt install htop` | 系统监控工具 |

### 硬件要求
- **内存**：至少 2 GB RAM
- **存储空间**：至少 5 GB 可用空间
- **CPU**：任意现代处理器
- **网络**：稳定的互联网连接（用于下载练习数据）

### 数据准备
| 数据文件 | 大小 | 下载链接/位置 | 说明 |
|---------|------|-------------|------|
| sample_sequences.fasta | ~1MB | 课程提供 | 示例DNA序列文件 |
| gene_expression.csv | ~500KB | 课程提供 | 基因表达数据 |
| quality_scores.txt | ~200KB | 课程提供 | 测序质量分数 |

## 操作步骤

### 步骤1：环境设置和准备工作

#### 1.1 创建工作目录
```bash
# 创建本次实验的工作目录
mkdir -p ~/linux-basics-lab
cd ~/linux-basics-lab

# 创建子目录结构
mkdir -p {data,scripts,results,logs,backup}

# 验证目录结构
tree
```

#### 1.2 检查系统环境
```bash
# 检查操作系统版本
cat /etc/os-release

# 检查当前用户
whoami
id

# 检查当前目录
pwd

# 检查磁盘空间
df -h

# 检查内存使用
free -h
```

**预期输出：**
```
NAME="Ubuntu"
VERSION="20.04.3 LTS (Focal Fossa)"
...
```

#### 1.3 下载和准备练习数据
```bash
# 进入数据目录
cd ~/linux-basics-lab/data

# 创建示例DNA序列文件
cat > sample_sequences.fasta << 'EOF'
>sequence1
ATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCG
>sequence2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
>sequence3
TTAATTAATTAATTAATTAATTAATTAATTAA
TTAATTAATTAATTAATTAATTAATTAATTAA
EOF

# 创建基因表达数据文件
cat > gene_expression.csv << 'EOF'
gene_id,sample1,sample2,sample3,sample4
GENE001,12.5,15.2,8.7,11.3
GENE002,45.1,42.8,48.3,44.7
GENE003,2.1,1.8,2.5,2.0
GENE004,78.9,82.1,75.6,79.4
GENE005,0.5,0.3,0.8,0.6
EOF

# 创建质量分数文件
cat > quality_scores.txt << 'EOF'
Read1: IIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
Read2: HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
Read3: GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
Read4: FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
Read5: EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
EOF

# 验证文件创建
ls -la
wc -l *.{fasta,csv,txt}
```

**检查点：** 确认所有数据文件已正确创建并位于 `data/` 目录中。

---

### 步骤2：基本文件操作练习

#### 2.1 文件查看和信息获取

**操作说明：**
学习使用各种命令查看文件内容和获取文件信息，这是处理生物信息学数据的基础技能。

**执行命令：**
```bash
# 回到工作目录
cd ~/linux-basics-lab

# 查看文件内容的不同方式
cat data/sample_sequences.fasta
echo "--- 分隔线 ---"
head -n 5 data/gene_expression.csv
echo "--- 分隔线 ---"
tail -n 3 data/quality_scores.txt

# 使用less分页查看（按q退出）
less data/sample_sequences.fasta

# 获取文件信息
ls -lh data/
file data/*
```

**参数解释：**
- `cat`：显示文件全部内容
- `head -n 5`：显示前5行
- `tail -n 3`：显示后3行
- `less`：分页查看，支持上下滚动
- `ls -lh`：详细列表，人类可读的文件大小
- `file`：显示文件类型

**预期输出：**
```
>sequence1
ATCGATCGATCGATCGATCGATCGATCGATCG
...
-rw-rw-r-- 1 user user 245 Dec 20 10:30 sample_sequences.fasta
```

#### 2.2 文件复制、移动和重命名

**执行命令：**
```bash
# 复制文件
cp data/sample_sequences.fasta backup/
cp data/gene_expression.csv backup/sequences_backup.csv

# 移动文件到results目录
cp data/quality_scores.txt results/
mv results/quality_scores.txt results/quality_backup.txt

# 验证操作结果
ls -la backup/
ls -la results/

# 创建符号链接
ln -s ../data/sample_sequences.fasta scripts/sequences_link.fasta
ls -la scripts/
```

**结果验证：**
```bash
# 验证文件是否正确复制
diff data/sample_sequences.fasta backup/sample_sequences.fasta
echo "如果没有输出，说明文件完全相同"

# 检查符号链接
readlink scripts/sequences_link.fasta
```

**检查点：** 确认文件已正确复制到backup目录，符号链接创建成功。

---

### 步骤3：文本处理和数据筛选

#### 3.1 使用grep搜索和筛选

**操作说明：**
grep是生物信息学中最常用的文本搜索工具，用于从大量数据中筛选特定信息。

**执行命令：**
```bash
# 搜索FASTA文件中的序列头
grep ">" data/sample_sequences.fasta

# 搜索包含特定基因的行
grep "GENE00[1-3]" data/gene_expression.csv

# 忽略大小写搜索
grep -i "gene" data/gene_expression.csv

# 显示行号
grep -n ">" data/sample_sequences.fasta

# 反向搜索（不包含指定模式的行）
grep -v ">" data/sample_sequences.fasta

# 统计匹配行数
grep -c ">" data/sample_sequences.fasta
```

**参数解释：**
- `grep "pattern"`：搜索包含模式的行
- `-i`：忽略大小写
- `-n`：显示行号
- `-v`：反向搜索
- `-c`：统计匹配行数

**预期输出：**
```
>sequence1
>sequence2
>sequence3
```

#### 3.2 使用cut和awk处理结构化数据

**执行命令：**
```bash
# 使用cut提取CSV文件的特定列
cut -d',' -f1,2 data/gene_expression.csv

# 提取基因ID列
cut -d',' -f1 data/gene_expression.csv | tail -n +2

# 使用awk处理数据
awk -F',' '{print $1, $2}' data/gene_expression.csv

# 计算第二列的平均值（跳过标题行）
awk -F',' 'NR>1 {sum+=$2; count++} END {print "平均值:", sum/count}' data/gene_expression.csv

# 筛选表达量大于40的基因
awk -F',' 'NR>1 && $2>40 {print $1, $2}' data/gene_expression.csv
```

**结果验证：**
```bash
# 验证提取的列数
cut -d',' -f1,2 data/gene_expression.csv | head -n 3
```

#### 3.3 排序和去重操作

**执行命令：**
```bash
# 创建测试数据
echo -e "apple\nbanana\napple\ncherry\nbanana" > results/fruits.txt

# 排序
sort results/fruits.txt

# 排序并去重
sort results/fruits.txt | uniq

# 统计重复次数
sort results/fruits.txt | uniq -c

# 对数值进行排序
cut -d',' -f2 data/gene_expression.csv | tail -n +2 | sort -n

# 找出表达量最高的基因
awk -F',' 'NR>1 {print $2, $1}' data/gene_expression.csv | sort -nr | head -n 1
```

**检查点：** 能够正确提取、排序和统计数据。

---

### 步骤4：管道和重定向操作

#### 4.1 输出重定向

**操作说明：**
重定向是将命令输出保存到文件的重要技术，在生物信息学分析中经常用于保存结果。

**执行命令：**
```bash
# 将结果重定向到文件
grep ">" data/sample_sequences.fasta > results/sequence_headers.txt

# 追加内容到文件
echo "分析完成时间: $(date)" >> results/sequence_headers.txt

# 同时保存输出和错误信息
ls data/ nonexistent_dir > results/ls_output.txt 2> results/ls_errors.txt

# 合并输出和错误
ls data/ nonexistent_dir > results/ls_combined.txt 2>&1

# 查看结果
cat results/sequence_headers.txt
cat results/ls_combined.txt
```

#### 4.2 管道操作

**执行命令：**
```bash
# 统计FASTA文件中的序列数量
grep ">" data/sample_sequences.fasta | wc -l

# 找出表达量最高的前3个基因
awk -F',' 'NR>1 {print $2, $1}' data/gene_expression.csv | sort -nr | head -n 3

# 复杂的管道操作：统计每个碱基的出现次数
grep -v ">" data/sample_sequences.fasta | tr -d '\n' | fold -w1 | sort | uniq -c | sort -nr

# 将结果保存到文件
grep -v ">" data/sample_sequences.fasta | tr -d '\n' | fold -w1 | sort | uniq -c | sort -nr > results/base_count.txt
```

**参数解释：**
- `tr -d '\n'`：删除换行符
- `fold -w1`：每行一个字符
- `uniq -c`：统计重复次数
- `sort -nr`：数值逆序排序

**检查点：** 成功使用管道组合多个命令，并将结果保存到文件。

---

### 步骤5：脚本编写和自动化

#### 5.1 创建第一个Bash脚本

**操作说明：**
编写脚本可以自动化重复的任务，这在处理大量生物信息学数据时非常重要。

**执行命令：**
```bash
# 创建脚本文件
cat > scripts/analyze_sequences.sh << 'EOF'
#!/bin/bash
# 序列分析脚本
# 作者：学生姓名
# 日期：$(date +%Y-%m-%d)

echo "开始分析序列文件..."

# 检查输入文件是否存在
if [ ! -f "data/sample_sequences.fasta" ]; then
    echo "错误：找不到输入文件"
    exit 1
fi

# 统计序列数量
seq_count=$(grep -c ">" data/sample_sequences.fasta)
echo "序列总数: $seq_count"

# 统计总碱基数
total_bases=$(grep -v ">" data/sample_sequences.fasta | tr -d '\n' | wc -c)
echo "总碱基数: $total_bases"

# 计算平均序列长度
avg_length=$((total_bases / seq_count))
echo "平均序列长度: $avg_length"

# 统计各碱基含量
echo "碱基组成统计:"
grep -v ">" data/sample_sequences.fasta | tr -d '\n' | fold -w1 | sort | uniq -c | sort -nr

echo "分析完成！"
EOF

# 给脚本添加执行权限
chmod +x scripts/analyze_sequences.sh

# 运行脚本
./scripts/analyze_sequences.sh
```

#### 5.2 创建数据处理脚本

**执行命令：**
```bash
# 创建基因表达数据处理脚本
cat > scripts/process_expression.sh << 'EOF'
#!/bin/bash
# 基因表达数据处理脚本

echo "处理基因表达数据..."

# 创建输出目录
mkdir -p results/expression_analysis

# 提取基因ID
cut -d',' -f1 data/gene_expression.csv | tail -n +2 > results/expression_analysis/gene_ids.txt

# 计算每个基因在所有样本中的平均表达量
awk -F',' 'NR>1 {
    avg = ($2 + $3 + $4 + $5) / 4
    print $1 "," avg
}' data/gene_expression.csv > results/expression_analysis/average_expression.csv

# 添加标题行
echo "gene_id,average_expression" > results/expression_analysis/gene_averages.csv
cat results/expression_analysis/average_expression.csv >> results/expression_analysis/gene_averages.csv

# 找出高表达基因（平均表达量>30）
echo "高表达基因（平均表达量>30）:"
awk -F',' '$2>30 {print $1, $2}' results/expression_analysis/average_expression.csv

# 统计信息
total_genes=$(wc -l < results/expression_analysis/gene_ids.txt)
high_expr_genes=$(awk -F',' '$2>30' results/expression_analysis/average_expression.csv | wc -l)

echo "总基因数: $total_genes"
echo "高表达基因数: $high_expr_genes"
echo "高表达基因比例: $(echo "scale=2; $high_expr_genes/$total_genes*100" | bc)%"

echo "数据处理完成！结果保存在 results/expression_analysis/ 目录中"
EOF

# 给脚本添加执行权限
chmod +x scripts/process_expression.sh

# 安装bc计算器（如果没有安装）
sudo apt install -y bc

# 运行脚本
./scripts/process_expression.sh
```

#### 5.3 创建系统监控脚本

**执行命令：**
```bash
# 创建系统监控脚本
cat > scripts/system_monitor.sh << 'EOF'
#!/bin/bash
# 系统监控脚本

echo "=== 系统监控报告 ==="
echo "生成时间: $(date)"
echo

# CPU信息
echo "=== CPU信息 ==="
lscpu | grep "Model name"
echo "CPU使用率:"
top -bn1 | grep "Cpu(s)" | awk '{print $2}' | cut -d'%' -f1

# 内存信息
echo
echo "=== 内存使用情况 ==="
free -h

# 磁盘使用情况
echo
echo "=== 磁盘使用情况 ==="
df -h | grep -E "^/dev/"

# 当前进程
echo
echo "=== 占用CPU最多的前5个进程 ==="
ps aux --sort=-%cpu | head -n 6

# 网络连接
echo
echo "=== 网络连接统计 ==="
netstat -tuln | wc -l
echo "总连接数: $(netstat -tuln | wc -l)"

echo
echo "=== 监控完成 ==="
EOF

# 给脚本添加执行权限
chmod +x scripts/system_monitor.sh

# 运行脚本
./scripts/system_monitor.sh > results/system_report.txt

# 查看报告
cat results/system_report.txt
```

**检查点：** 成功创建并运行了三个不同功能的脚本。

---

### 步骤6：进程管理和任务调度

#### 6.1 进程管理练习

**操作说明：**
学习管理长时间运行的生物信息学分析任务，这在处理大数据集时非常重要。

**执行命令：**
```bash
# 创建一个模拟长时间运行的任务
cat > scripts/long_task.sh << 'EOF'
#!/bin/bash
# 模拟长时间运行的分析任务

echo "开始长时间分析任务..."
for i in {1..30}; do
    echo "处理进度: $i/30"
    sleep 2
done
echo "任务完成！"
EOF

chmod +x scripts/long_task.sh

# 后台运行任务
./scripts/long_task.sh &

# 查看后台任务
jobs

# 查看进程
ps aux | grep long_task

# 将任务调到前台（可选）
# fg %1

# 如果需要终止任务
# kill %1
```

#### 6.2 使用nohup运行持久任务

**执行命令：**
```bash
# 使用nohup运行任务，即使终端关闭也会继续运行
nohup ./scripts/long_task.sh > results/long_task.log 2>&1 &

# 查看任务状态
jobs
ps aux | grep long_task

# 查看日志
tail -f results/long_task.log
# 按Ctrl+C退出tail命令
```

**检查点：** 理解如何管理后台进程和长时间运行的任务。

---

### 步骤7：文件权限和安全

#### 7.1 权限管理练习

**执行命令：**
```bash
# 查看当前文件权限
ls -la scripts/

# 修改脚本权限
chmod 755 scripts/*.sh

# 创建只读数据文件
cp data/sample_sequences.fasta results/readonly_sequences.fasta
chmod 444 results/readonly_sequences.fasta

# 尝试修改只读文件（应该失败）
echo "test" >> results/readonly_sequences.fasta

# 查看权限变化
ls -la results/readonly_sequences.fasta

# 恢复写权限
chmod 644 results/readonly_sequences.fasta

# 创建组共享目录
mkdir -p results/shared_data
chmod 775 results/shared_data
```

#### 7.2 文件所有权管理

**执行命令：**
```bash
# 查看文件所有者
ls -la data/

# 查看当前用户和组
id

# 创建备份并设置权限
cp -r data/ backup/data_backup
chmod -R 755 backup/data_backup

# 验证权限设置
ls -la backup/data_backup/
```

**检查点：** 理解Linux文件权限系统和安全管理。

---

## 预期结果

### 主要输出文件
1. **分析脚本**：`scripts/analyze_sequences.sh`
   - 内容：序列分析自动化脚本
   - 用途：统计FASTA文件中的序列信息

2. **处理结果**：`results/expression_analysis/`
   - 内容：基因表达数据分析结果
   - 用途：展示数据处理和统计分析能力

3. **系统报告**：`results/system_report.txt`
   - 内容：系统资源使用情况
   - 用途：系统监控和性能分析

### 关键结果指标
- 序列统计：应该正确统计出3个序列
- 碱基组成：A、T、C、G的数量统计准确
- 基因表达：正确计算平均表达量和筛选高表达基因

### 成功标准
- [ ] 所有脚本执行无错误
- [ ] 生成了预期的输出文件
- [ ] 数据处理结果准确
- [ ] 能够熟练使用命令行操作

## 故障排除

### 常见问题1：权限被拒绝
**症状：** `Permission denied` 错误
**原因：** 脚本文件没有执行权限
**解决方案：**
```bash
# 添加执行权限
chmod +x script_name.sh
```

### 常见问题2：命令未找到
**症状：** `command not found` 错误
**原因：** 软件包未安装或不在PATH中
**解决方案：**
```bash
# 安装缺失的软件包
sudo apt update
sudo apt install package_name

# 检查PATH环境变量
echo $PATH
```

### 常见问题3：文件不存在
**症状：** `No such file or directory` 错误
**原因：** 文件路径错误或文件未创建
**解决方案：**
```bash
# 检查当前目录
pwd
# 检查文件是否存在
ls -la filename
# 使用绝对路径
/full/path/to/file
```

### 获取帮助
如果遇到其他问题：
1. 查看命令帮助：`command --help` 或 `man command`
2. 检查错误日志：`cat logs/error.log`
3. 联系助教或老师：wangys@hunau.edu.cn

## 扩展练习

### 练习1：高级文本处理
**目标：** 掌握更复杂的文本处理技术
**任务：** 
1. 从FASTA文件中提取所有序列（不包括标题行）
2. 计算每个序列的GC含量
3. 生成GC含量统计报告

**提示：** 使用awk和正则表达式

### 练习2：批量文件处理
**目标：** 学会批量处理多个文件
**任务：**
1. 创建10个模拟的基因表达文件
2. 编写脚本批量处理所有文件
3. 生成汇总统计报告

**提示：** 使用for循环和通配符

### 练习3：日志分析
**目标：** 分析系统日志文件
**任务：**
1. 分析系统日志文件（/var/log/syslog）
2. 统计不同类型的日志条目
3. 找出最近的错误信息

**提示：** 使用grep、awk和时间处理

### 思考问题
1. 为什么生物信息学分析通常在Linux系统上进行？
2. 如何优化大文件的处理效率？
3. 在处理敏感的生物数据时，应该注意哪些安全问题？

## 参考资料

### 相关文献
1. Linux Command Line and Shell Scripting Bible. Richard Blum & Christine Bresnahan.
2. Bioinformatics Data Skills. Vince Buffalo. O'Reilly Media, 2015.

### 在线资源
- Linux命令大全：https://www.runoob.com/linux/linux-command-manual.html
- Bash脚本教程：https://wangdoc.com/bash/
- 生物信息学Linux教程：https://datacarpentry.org/shell-genomics/

### 软件文档
- Bash手册：https://www.gnu.org/software/bash/manual/
- GNU核心工具：https://www.gnu.org/software/coreutils/manual/

## 附录

### 附录A：常用命令速查表
```bash
# 文件操作
ls -la          # 详细列表
cp -r src dst   # 递归复制
mv src dst      # 移动/重命名
rm -rf dir      # 强制删除目录

# 文本处理
grep pattern file    # 搜索文本
sed 's/old/new/g'   # 替换文本
awk '{print $1}'    # 提取列
sort | uniq         # 排序去重

# 系统信息
ps aux         # 进程列表
top            # 实时监控
df -h          # 磁盘使用
free -h        # 内存使用
```

### 附录B：脚本模板
```bash
#!/bin/bash
# 脚本描述
# 作者：姓名
# 日期：日期

set -e  # 遇到错误立即退出
set -u  # 使用未定义变量时报错

# 参数检查
if [ $# -lt 1 ]; then
    echo "用法: $0 <参数>"
    exit 1
fi

# 主要功能
main() {
    echo "开始执行..."
    # 具体操作
    echo "执行完成"
}

main "$@"
```

### 附录C：FASTA格式说明
FASTA格式是生物信息学中最常用的序列格式：
- 标题行以 `>` 开头，包含序列标识符和描述
- 序列行包含实际的DNA、RNA或蛋白质序列
- 每行通常不超过80个字符

---

**实验完成时间：** 预计 3 小时  
**难度等级：** 初级  
**最后更新：** 2025年12月20日