#!/bin/bash
# Linux基础课程练习脚本
# 作者：王运生
# 日期：2025-12-20
# 用法：bash practice_exercises.sh

set -e  # 遇到错误立即退出

echo "=== Linux基础课程练习脚本 ==="
echo "开始时间: $(date)"
echo

# 检查工作目录
WORK_DIR="$HOME/linux-basics-lab"
if [ ! -d "$WORK_DIR" ]; then
    echo "错误：工作目录不存在，请先运行 setup.sh"
    exit 1
fi

cd "$WORK_DIR"

# 练习1：文件操作基础
echo "=== 练习1：文件操作基础 ==="
echo "1.1 创建练习目录..."
mkdir -p exercises/file_operations
cd exercises/file_operations

echo "1.2 创建测试文件..."
echo "这是测试文件1" > test1.txt
echo -e "apple\nbanana\ncherry\napple\ndate" > fruits.txt
echo -e "1\n5\n3\n9\n2\n7" > numbers.txt

echo "1.3 文件复制和移动练习..."
cp test1.txt test1_backup.txt
cp fruits.txt fruits_copy.txt
mkdir backup
mv test1_backup.txt backup/

echo "1.4 文件权限练习..."
chmod 755 test1.txt
chmod 644 fruits.txt
ls -la *.txt

echo "练习1完成！"
echo

# 练习2：文本处理
echo "=== 练习2：文本处理 ==="
cd "$WORK_DIR"
mkdir -p exercises/text_processing
cd exercises/text_processing

echo "2.1 创建文本处理测试数据..."
cat > sample_data.txt << 'EOF'
ID001,John,25,Engineer,50000
ID002,Alice,30,Manager,75000
ID003,Bob,28,Developer,60000
ID004,Carol,35,Designer,55000
ID005,David,32,Analyst,65000
EOF

echo "2.2 使用grep搜索..."
echo "搜索包含'Manager'的行："
grep "Manager" sample_data.txt

echo "2.3 使用cut提取字段..."
echo "提取姓名和职位："
cut -d',' -f2,4 sample_data.txt

echo "2.4 使用awk处理数据..."
echo "计算平均薪资："
awk -F',' '{sum+=$5; count++} END {print "平均薪资:", sum/count}' sample_data.txt

echo "2.5 排序练习..."
echo "按薪资排序："
sort -t',' -k5 -n sample_data.txt

echo "练习2完成！"
echo

# 练习3：管道和重定向
echo "=== 练习3：管道和重定向 ==="
cd "$WORK_DIR"
mkdir -p exercises/pipes_redirection
cd exercises/pipes_redirection

echo "3.1 创建日志文件..."
cat > access.log << 'EOF'
192.168.1.1 - - [20/Dec/2025:10:00:01] "GET /index.html HTTP/1.1" 200 1234
192.168.1.2 - - [20/Dec/2025:10:00:02] "POST /login HTTP/1.1" 200 567
192.168.1.1 - - [20/Dec/2025:10:00:03] "GET /about.html HTTP/1.1" 404 0
192.168.1.3 - - [20/Dec/2025:10:00:04] "GET /contact.html HTTP/1.1" 200 890
192.168.1.2 - - [20/Dec/2025:10:00:05] "GET /products.html HTTP/1.1" 200 2345
EOF

echo "3.2 管道操作练习..."
echo "统计不同状态码的数量："
awk '{print $9}' access.log | sort | uniq -c

echo "3.3 重定向练习..."
echo "将结果保存到文件："
awk '{print $9}' access.log | sort | uniq -c > status_codes.txt
cat status_codes.txt

echo "3.4 复杂管道操作..."
echo "找出访问最多的IP地址："
awk '{print $1}' access.log | sort | uniq -c | sort -nr | head -n 1

echo "练习3完成！"
echo

# 练习4：系统信息收集
echo "=== 练习4：系统信息收集 ==="
cd "$WORK_DIR"
mkdir -p exercises/system_info
cd exercises/system_info

echo "4.1 收集系统信息..."
cat > system_report.txt << EOF
=== 系统信息报告 ===
生成时间: $(date)
用户名: $(whoami)
主机名: $(hostname)
操作系统: $(cat /etc/os-release | grep PRETTY_NAME | cut -d'"' -f2)
内核版本: $(uname -r)
运行时间: $(uptime)

=== 硬件信息 ===
CPU信息: $(lscpu | grep "Model name" | cut -d':' -f2 | xargs)
内存信息:
$(free -h)

=== 磁盘使用情况 ===
$(df -h | grep -E "^/dev/")

=== 网络接口 ===
$(ip addr show | grep -E "^[0-9]+:" | cut -d':' -f2 | xargs)

=== 当前进程数 ===
总进程数: $(ps aux | wc -l)
EOF

echo "4.2 显示系统报告..."
cat system_report.txt

echo "练习4完成！"
echo

# 练习5：脚本编程
echo "=== 练习5：脚本编程 ==="
cd "$WORK_DIR"
mkdir -p exercises/scripting
cd exercises/scripting

echo "5.1 创建文件统计脚本..."
cat > file_stats.sh << 'EOF'
#!/bin/bash
# 文件统计脚本

if [ $# -eq 0 ]; then
    echo "用法: $0 <目录路径>"
    exit 1
fi

DIR="$1"
if [ ! -d "$DIR" ]; then
    echo "错误: $DIR 不是一个有效目录"
    exit 1
fi

echo "=== 目录统计: $DIR ==="
echo "文件总数: $(find "$DIR" -type f | wc -l)"
echo "目录总数: $(find "$DIR" -type d | wc -l)"
echo "总大小: $(du -sh "$DIR" | cut -f1)"

echo
echo "=== 文件类型统计 ==="
find "$DIR" -type f -name "*.*" | sed 's/.*\.//' | sort | uniq -c | sort -nr | head -5

echo
echo "=== 最大的5个文件 ==="
find "$DIR" -type f -exec ls -lh {} \; | sort -k5 -hr | head -5 | awk '{print $5, $9}'
EOF

chmod +x file_stats.sh

echo "5.2 测试文件统计脚本..."
./file_stats.sh "$WORK_DIR/data"

echo "5.3 创建备份脚本..."
cat > backup_script.sh << 'EOF'
#!/bin/bash
# 简单备份脚本

SOURCE_DIR="$1"
BACKUP_DIR="$2"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

if [ $# -ne 2 ]; then
    echo "用法: $0 <源目录> <备份目录>"
    exit 1
fi

if [ ! -d "$SOURCE_DIR" ]; then
    echo "错误: 源目录 $SOURCE_DIR 不存在"
    exit 1
fi

mkdir -p "$BACKUP_DIR"
BACKUP_NAME="backup_${TIMESTAMP}.tar.gz"

echo "开始备份 $SOURCE_DIR 到 $BACKUP_DIR/$BACKUP_NAME"
tar -czf "$BACKUP_DIR/$BACKUP_NAME" -C "$(dirname "$SOURCE_DIR")" "$(basename "$SOURCE_DIR")"

if [ $? -eq 0 ]; then
    echo "备份成功完成: $BACKUP_DIR/$BACKUP_NAME"
    echo "备份大小: $(ls -lh "$BACKUP_DIR/$BACKUP_NAME" | awk '{print $5}')"
else
    echo "备份失败"
    exit 1
fi
EOF

chmod +x backup_script.sh

echo "5.4 测试备份脚本..."
mkdir -p backups
./backup_script.sh "$WORK_DIR/data" "$WORK_DIR/exercises/scripting/backups"

echo "练习5完成！"
echo

# 练习6：高级文本处理
echo "=== 练习6：高级文本处理 ==="
cd "$WORK_DIR"
mkdir -p exercises/advanced_text
cd exercises/advanced_text

echo "6.1 处理FASTA文件..."
cp "$WORK_DIR/data/sample_sequences.fasta" .

echo "统计序列信息："
echo "序列总数: $(grep -c ">" sample_sequences.fasta)"
echo "总碱基数: $(grep -v ">" sample_sequences.fasta | tr -d '\n' | wc -c)"

echo "6.2 计算GC含量..."
cat > calculate_gc.sh << 'EOF'
#!/bin/bash
# 计算FASTA文件中每个序列的GC含量

FASTA_FILE="$1"
if [ ! -f "$FASTA_FILE" ]; then
    echo "错误: FASTA文件不存在"
    exit 1
fi

echo "序列ID,长度,GC含量(%)"
awk '
/^>/ {
    if (seq_id != "") {
        gc_count = gsub(/[GCgc]/, "", sequence)
        total_length = length(sequence)
        gc_percent = (gc_count / total_length) * 100
        printf "%s,%d,%.2f\n", seq_id, total_length, gc_percent
    }
    seq_id = substr($0, 2)
    sequence = ""
}
!/^>/ {
    sequence = sequence $0
}
END {
    if (seq_id != "") {
        gc_count = gsub(/[GCgc]/, "", sequence)
        total_length = length(sequence)
        gc_percent = (gc_count / total_length) * 100
        printf "%s,%d,%.2f\n", seq_id, total_length, gc_percent
    }
}' "$FASTA_FILE"
EOF

chmod +x calculate_gc.sh
./calculate_gc.sh sample_sequences.fasta > gc_content.csv

echo "6.3 GC含量统计结果："
cat gc_content.csv

echo "练习6完成！"
echo

# 生成练习总结报告
echo "=== 生成练习总结报告 ==="
cd "$WORK_DIR"
cat > exercises/practice_summary.txt << EOF
=== Linux基础课程练习总结 ===
完成时间: $(date)
学生: $(whoami)

练习完成情况:
✓ 练习1: 文件操作基础
✓ 练习2: 文本处理
✓ 练习3: 管道和重定向
✓ 练习4: 系统信息收集
✓ 练习5: 脚本编程
✓ 练习6: 高级文本处理

创建的文件和目录:
$(find exercises/ -type f | wc -l) 个文件
$(find exercises/ -type d | wc -l) 个目录

总用时: 约2-3小时
难度评估: 初级到中级

建议后续学习:
1. 深入学习正则表达式
2. 学习更多的文本处理工具（sed, awk高级用法）
3. 学习Shell脚本高级编程
4. 学习系统管理和网络操作

EOF

echo "练习总结报告已生成: exercises/practice_summary.txt"
echo
echo "=== 所有练习完成！==="
echo "请查看各个练习目录中的结果文件"
echo "练习目录结构："
tree exercises/ -L 2
echo
echo "完成时间: $(date)"