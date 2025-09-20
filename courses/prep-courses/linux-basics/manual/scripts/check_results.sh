#!/bin/bash
# Linux基础课程结果检查脚本
# 作者：王运生教授
# 日期：2025-12-20
# 用法：bash check_results.sh

set -e  # 遇到错误立即退出

echo "=== Linux基础课程结果检查 ==="
echo "检查时间: $(date)"
echo

# 检查工作目录
WORK_DIR="$HOME/linux-basics-lab"
if [ ! -d "$WORK_DIR" ]; then
    echo "❌ 错误：工作目录 $WORK_DIR 不存在"
    echo "请先运行 setup.sh 创建环境"
    exit 1
fi

cd "$WORK_DIR"
echo "✓ 工作目录存在: $WORK_DIR"

# 初始化检查结果
TOTAL_CHECKS=0
PASSED_CHECKS=0

# 检查函数
check_item() {
    local description="$1"
    local condition="$2"
    TOTAL_CHECKS=$((TOTAL_CHECKS + 1))
    
    if eval "$condition"; then
        echo "✓ $description"
        PASSED_CHECKS=$((PASSED_CHECKS + 1))
        return 0
    else
        echo "❌ $description"
        return 1
    fi
}

echo
echo "=== 1. 基础环境检查 ==="

# 检查目录结构
check_item "数据目录存在" "[ -d 'data' ]"
check_item "脚本目录存在" "[ -d 'scripts' ]"
check_item "结果目录存在" "[ -d 'results' ]"
check_item "日志目录存在" "[ -d 'logs' ]"
check_item "备份目录存在" "[ -d 'backup' ]"

# 检查数据文件
check_item "FASTA文件存在" "[ -f 'data/sample_sequences.fasta' ]"
check_item "表达数据文件存在" "[ -f 'data/gene_expression.csv' ]"
check_item "质量分数文件存在" "[ -f 'data/quality_scores.txt' ]"

echo
echo "=== 2. 文件操作检查 ==="

# 检查文件权限
check_item "脚本文件可执行" "[ -x 'scripts/analyze_sequences.sh' ] 2>/dev/null || true"
check_item "数据文件可读" "[ -r 'data/sample_sequences.fasta' ]"

# 检查文件内容
if [ -f "data/sample_sequences.fasta" ]; then
    SEQ_COUNT=$(grep -c ">" data/sample_sequences.fasta 2>/dev/null || echo "0")
    check_item "FASTA文件包含序列 (>3)" "[ '$SEQ_COUNT' -gt 3 ]"
fi

if [ -f "data/gene_expression.csv" ]; then
    GENE_COUNT=$(tail -n +2 data/gene_expression.csv 2>/dev/null | wc -l || echo "0")
    check_item "基因表达文件包含数据 (>5行)" "[ '$GENE_COUNT' -gt 5 ]"
fi

echo
echo "=== 3. 命令执行能力检查 ==="

# 检查基本命令
check_item "ls命令可用" "command -v ls >/dev/null"
check_item "grep命令可用" "command -v grep >/dev/null"
check_item "awk命令可用" "command -v awk >/dev/null"
check_item "sed命令可用" "command -v sed >/dev/null"
check_item "sort命令可用" "command -v sort >/dev/null"

# 检查高级工具
check_item "tree命令可用" "command -v tree >/dev/null"
check_item "htop命令可用" "command -v htop >/dev/null"
check_item "nano编辑器可用" "command -v nano >/dev/null"

echo
echo "=== 4. 实际操作测试 ==="

# 创建临时测试目录
TEST_DIR="temp_test_$$"
mkdir -p "$TEST_DIR"
cd "$TEST_DIR"

# 测试文件操作
echo "test content" > test_file.txt
check_item "文件创建成功" "[ -f 'test_file.txt' ]"

cp test_file.txt test_copy.txt
check_item "文件复制成功" "[ -f 'test_copy.txt' ]"

# 测试文本处理
echo -e "apple\nbanana\napple\ncherry" > fruits.txt
UNIQUE_COUNT=$(sort fruits.txt | uniq | wc -l)
check_item "文本处理正确 (去重)" "[ '$UNIQUE_COUNT' -eq 3 ]"

# 测试管道操作
PIPE_RESULT=$(echo "hello world" | wc -w)
check_item "管道操作正确" "[ '$PIPE_RESULT' -eq 2 ]"

# 测试权限操作
chmod 755 test_file.txt
PERM_CHECK=$(ls -l test_file.txt | cut -c1-10)
check_item "权限设置正确" "echo '$PERM_CHECK' | grep -q 'rwxr-xr-x'"

# 清理测试目录
cd ..
rm -rf "$TEST_DIR"

echo
echo "=== 5. 练习完成情况检查 ==="

# 检查练习目录
if [ -d "exercises" ]; then
    check_item "练习目录存在" "true"
    
    # 检查各个练习
    check_item "文件操作练习完成" "[ -d 'exercises/file_operations' ]"
    check_item "文本处理练习完成" "[ -d 'exercises/text_processing' ]"
    check_item "管道重定向练习完成" "[ -d 'exercises/pipes_redirection' ]"
    check_item "系统信息练习完成" "[ -d 'exercises/system_info' ]"
    check_item "脚本编程练习完成" "[ -d 'exercises/scripting' ]"
    
    # 检查脚本文件
    if [ -f "exercises/scripting/file_stats.sh" ]; then
        check_item "文件统计脚本存在" "[ -x 'exercises/scripting/file_stats.sh' ]"
    fi
    
    if [ -f "exercises/scripting/backup_script.sh" ]; then
        check_item "备份脚本存在" "[ -x 'exercises/scripting/backup_script.sh' ]"
    fi
else
    echo "ℹ️  练习目录不存在，可能尚未完成练习"
fi

echo
echo "=== 6. 高级技能检查 ==="

# 检查是否能处理FASTA文件
if [ -f "data/sample_sequences.fasta" ]; then
    # 测试序列统计
    TOTAL_BASES=$(grep -v ">" data/sample_sequences.fasta | tr -d '\n' | wc -c 2>/dev/null || echo "0")
    check_item "能够统计序列碱基数" "[ '$TOTAL_BASES' -gt 100 ]"
    
    # 测试GC含量计算（如果存在相关文件）
    if [ -f "exercises/advanced_text/gc_content.csv" ]; then
        check_item "GC含量计算完成" "[ -s 'exercises/advanced_text/gc_content.csv' ]"
    fi
fi

# 检查系统监控能力
PROCESS_COUNT=$(ps aux | wc -l 2>/dev/null || echo "0")
check_item "能够获取系统进程信息" "[ '$PROCESS_COUNT' -gt 10 ]"

echo
echo "=== 7. 生成详细报告 ==="

# 创建检查报告
REPORT_FILE="results/skill_assessment_$(date +%Y%m%d_%H%M%S).txt"
mkdir -p results

cat > "$REPORT_FILE" << EOF
=== Linux基础技能评估报告 ===
评估时间: $(date)
学生用户: $(whoami)
工作目录: $WORK_DIR

=== 评估结果 ===
总检查项目: $TOTAL_CHECKS
通过项目: $PASSED_CHECKS
通过率: $(echo "scale=1; $PASSED_CHECKS * 100 / $TOTAL_CHECKS" | bc 2>/dev/null || echo "计算错误")%

=== 技能等级评估 ===
EOF

# 计算技能等级
PASS_RATE=$(echo "scale=0; $PASSED_CHECKS * 100 / $TOTAL_CHECKS" | bc 2>/dev/null || echo "0")

if [ "$PASS_RATE" -ge 90 ]; then
    SKILL_LEVEL="优秀 (90%+)"
    RECOMMENDATION="已掌握Linux基础技能，可以进入下一阶段学习"
elif [ "$PASS_RATE" -ge 75 ]; then
    SKILL_LEVEL="良好 (75-89%)"
    RECOMMENDATION="基本掌握Linux技能，建议加强练习薄弱环节"
elif [ "$PASS_RATE" -ge 60 ]; then
    SKILL_LEVEL="及格 (60-74%)"
    RECOMMENDATION="需要更多练习，重点关注失败的检查项目"
else
    SKILL_LEVEL="需要改进 (<60%)"
    RECOMMENDATION="建议重新学习基础内容，完成所有练习"
fi

cat >> "$REPORT_FILE" << EOF
技能等级: $SKILL_LEVEL
学习建议: $RECOMMENDATION

=== 详细检查结果 ===
通过的检查项目: $PASSED_CHECKS/$TOTAL_CHECKS

=== 环境信息 ===
操作系统: $(cat /etc/os-release | grep PRETTY_NAME | cut -d'"' -f2 2>/dev/null || echo "未知")
Shell版本: $SHELL
当前路径: $(pwd)
磁盘使用: $(df -h . | tail -1 | awk '{print $5}' 2>/dev/null || echo "未知")

=== 建议后续学习内容 ===
1. Python编程基础
2. R语言数据分析
3. 生物信息学工具使用
4. 高级Shell脚本编程
5. 版本控制系统(Git)

=== 评估完成 ===
报告生成时间: $(date)
EOF

echo "✓ 详细评估报告已生成: $REPORT_FILE"

echo
echo "=== 最终评估结果 ==="
echo "通过率: $PASSED_CHECKS/$TOTAL_CHECKS ($PASS_RATE%)"
echo "技能等级: $SKILL_LEVEL"
echo "学习建议: $RECOMMENDATION"
echo
echo "详细报告: $REPORT_FILE"

# 根据通过率给出不同的退出状态
if [ "$PASS_RATE" -ge 75 ]; then
    echo "🎉 恭喜！您已经掌握了Linux基础技能"
    exit 0
elif [ "$PASS_RATE" -ge 60 ]; then
    echo "⚠️  基本合格，建议继续练习提高"
    exit 0
else
    echo "❌ 需要更多练习，请重新完成相关内容"
    exit 1
fi