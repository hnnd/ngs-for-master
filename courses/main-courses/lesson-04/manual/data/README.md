# Lesson-04 数据文件说明

## 数据准备

本课程的数据文件通过自动化脚本准备，请运行：

```bash
cd ~/ngs-analysis/lesson-04
bash scripts/download_data.sh
```

## 数据文件位置

**实际数据存储位置:** `~/ngs-analysis/lesson-04/data/`

所有实验数据将自动下载/生成到上述目录，包括：

1. **chr22_fragment.fa** - 参考基因组 (人类22号染色体5Mb片段)
2. **chr22_fragment.fa.fai** - 参考基因组索引
3. **sample.bam** - 已比对的测序数据 (~10万reads)
4. **sample.bam.bai** - BAM文件索引
5. **dbsnp.vcf** - 已知变异数据库 (~1000个位点)
6. **hapmap.vcf** - 高质量SNP集合 (~500个位点)

## 目录结构

- `sample_data/` - 保留用于未来扩展
- `expected_results/` - 保留用于存放标准结果

## 数据详细说明

详细的数据来源、格式和使用说明请参阅：
- **[DATA_SOURCES.md](../DATA_SOURCES.md)** - 完整的数据来源文档
- **[manual.md](../manual.md)** - 实验操作手册

## 数据特征

- **参考基因组大小:** ~5 MB (GRCh38 chr22:10M-15M)
- **测序数据量:** ~100,000 paired reads (~200,000 total)
- **覆盖深度:** ~30X
- **文件总大小:** ~150-200 MB

## 快速验证

验证数据是否准备就绪：

```bash
cd ~/ngs-analysis/lesson-04/data

# 检查所有必需文件
for file in chr22_fragment.fa chr22_fragment.fa.fai sample.bam sample.bam.bai dbsnp.vcf hapmap.vcf; do
    if [ -f "$file" ]; then
        echo "✓ $file"
    else
        echo "✗ $file 缺失"
    fi
done

# 查看数据统计
echo -e "\n数据统计:"
echo "参考基因组: $(grep -v '^>' chr22_fragment.fa | wc -c) bp"
echo "BAM文件: $(samtools view -c sample.bam 2>/dev/null) reads"
echo "dbSNP变异: $(grep -v '^#' dbsnp.vcf | wc -l) 个"
echo "HapMap变异: $(grep -v '^#' hapmap.vcf | wc -l) 个"
```

## 使用说明

1. **首次使用:** 运行 `bash scripts/download_data.sh` 准备数据
2. **验证数据:** 使用上面的快速验证命令检查
3. **开始实验:** 按照 [manual.md](../manual.md) 进行操作
4. **结果保存:** 所有分析结果保存在 `~/ngs-analysis/lesson-04/results/`

## 注意事项

⚠️ **重要提示:**
- 这是教学用数据，规模远小于真实全基因组数据
- dbSNP和HapMap是模拟生成的，重叠率会低于真实数据
- 不要用于实际研究或临床诊断
- 请勿修改 `~/ngs-analysis/lesson-04/data/` 目录下的原始文件
- 分析结果自动保存在 `results/` 目录
- 定期清理临时文件释放空间

## 故障排除

**数据缺失:**
```bash
# 重新运行数据准备脚本
cd ~/ngs-analysis/lesson-04
bash scripts/download_data.sh
```

**空间不足:**
```bash
# 检查可用空间
df -h ~

# 清理旧的结果（谨慎操作）
rm -rf ~/ngs-analysis/lesson-04/results/*
```

**文件损坏:**
```bash
# 删除特定文件并重新生成
cd ~/ngs-analysis/lesson-04/data
rm chr22_fragment.fa*
cd ..
bash scripts/download_data.sh
```

## 获取帮助

如遇问题，请联系：
- **教师:** 王运生
- **邮箱:** wangys@hunau.edu.cn
- **办公室:** 16教420室

---

**最后更新:** 2025年1月
