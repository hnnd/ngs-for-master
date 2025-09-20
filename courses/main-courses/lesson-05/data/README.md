# RNA-seq示例数据说明

## 数据概述

本目录包含用于第5次课转录组测序数据分析的示例数据。这些数据是经过处理的小规模测试数据，用于教学演示。

## 数据来源

- **物种**: 人类 (Homo sapiens)
- **测序平台**: Illumina HiSeq
- **测序类型**: 双端测序 (Paired-end)
- **读长**: 150bp
- **链特异性**: 是 (stranded)

## 文件结构

```
data/
├── README.md                    # 本文件
├── sample1_R1.fastq.gz         # 处理组样本1，正向reads
├── sample1_R2.fastq.gz         # 处理组样本1，反向reads
├── sample2_R1.fastq.gz         # 处理组样本2，正向reads
├── sample2_R2.fastq.gz         # 处理组样本2，反向reads
├── sample3_R1.fastq.gz         # 处理组样本3，正向reads
├── sample3_R2.fastq.gz         # 处理组样本3，反向reads
├── ctrl1_R1.fastq.gz           # 对照组样本1，正向reads
├── ctrl1_R2.fastq.gz           # 对照组样本1，反向reads
├── ctrl2_R1.fastq.gz           # 对照组样本2，正向reads
├── ctrl2_R2.fastq.gz           # 对照组样本2，反向reads
├── ctrl3_R1.fastq.gz           # 对照组样本3，正向reads
├── ctrl3_R2.fastq.gz           # 对照组样本3，反向reads
├── sample_info.txt             # 样本信息表
└── clean/                      # 清洗后数据目录（分析过程中生成）
```

## 样本信息

| 样本名称 | 分组 | 生物学重复 | 描述 |
|---------|------|-----------|------|
| sample1 | treatment | 1 | 处理组重复1 |
| sample2 | treatment | 2 | 处理组重复2 |
| sample3 | treatment | 3 | 处理组重复3 |
| ctrl1   | control   | 1 | 对照组重复1 |
| ctrl2   | control   | 2 | 对照组重复2 |
| ctrl3   | control   | 3 | 对照组重复3 |

## 实验设计

- **对照组**: 未处理的细胞
- **处理组**: 经过特定处理的细胞
- **生物学重复**: 每组3个独立的生物学重复
- **技术重复**: 无（每个生物学样本只测序一次）

## 数据特征

- **文件大小**: 每个FASTQ文件约100-200MB（压缩后）
- **reads数量**: 每个样本约100万条reads
- **质量**: 平均质量值 > Q30
- **GC含量**: 约45-50%
- **重复序列**: 低水平

## 下载说明

由于文件较大，实际的FASTQ文件需要从课程服务器下载：

```bash
# 创建数据目录
mkdir -p ~/ngs-analysis/lesson-05/data

# 从课程服务器下载数据
scp user@server:/shared/course_data/rnaseq/*.fastq.gz ~/ngs-analysis/lesson-05/data/

# 或者使用wget下载（如果有公开链接）
# wget -P ~/ngs-analysis/lesson-05/data/ http://course-server.edu/data/rnaseq/sample1_R1.fastq.gz
```

## 质量控制预期结果

运行FastQC后，预期看到：

- **Per base sequence quality**: 大部分位置质量值 > 28
- **Per sequence quality scores**: 峰值在质量值35左右
- **Per base sequence content**: 前几个位置可能有偏差（正常现象）
- **Per sequence GC content**: 正态分布，峰值约47%
- **Sequence Length Distribution**: 150bp
- **Sequence Duplication Levels**: 低重复水平
- **Overrepresented sequences**: 可能有少量接头序列

## 比对预期结果

使用HISAT2比对到人类基因组GRCh38后，预期：

- **总体比对率**: 85-95%
- **唯一比对**: 80-90%
- **多重比对**: 5-10%
- **未比对**: 5-15%

## 基因定量预期结果

使用featureCounts定量后，预期：

- **成功分配的reads**: 70-85%
- **检测到的基因数量**: 15,000-20,000个
- **高表达基因**: 管家基因（GAPDH, ACTB等）
- **低表达基因**: 组织特异性基因

## 差异表达预期结果

DESeq2分析后，预期：

- **显著差异基因**: 1,000-3,000个（padj < 0.05, |FC| > 2）
- **上调基因**: 约占差异基因的40-60%
- **下调基因**: 约占差异基因的40-60%
- **最大倍数变化**: 5-20倍

## 功能富集预期结果

GO和KEGG分析后，预期富集的通路：

- **细胞周期相关**: 如果处理影响细胞增殖
- **代谢通路**: 如果处理影响细胞代谢
- **信号转导**: 如果处理激活特定信号通路
- **免疫反应**: 如果处理引起免疫应答

## 注意事项

1. **文件完整性**: 下载后请检查文件大小和MD5值
2. **存储空间**: 确保有足够的磁盘空间（至少20GB）
3. **内存要求**: 建议至少8GB RAM
4. **处理时间**: 完整分析预计需要2-4小时

## 故障排除

### 文件下载问题
```bash
# 检查文件完整性
md5sum *.fastq.gz
# 重新下载损坏的文件
```

### 文件权限问题
```bash
# 修改文件权限
chmod 644 *.fastq.gz
```

### 磁盘空间不足
```bash
# 检查磁盘空间
df -h
# 清理临时文件
rm -rf /tmp/*
```

## 联系方式

如有问题，请联系：
- **教师**: 王运生
- **邮箱**: wangys@hunau.edu.cn
- **办公室**: 16教420室

## 更新日志

- **2024-XX-XX**: 初始版本创建
- **2024-XX-XX**: 添加质量控制说明
- **2024-XX-XX**: 更新预期结果描述