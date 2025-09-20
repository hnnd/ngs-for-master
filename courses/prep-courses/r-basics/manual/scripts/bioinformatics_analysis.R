#!/usr/bin/env Rscript
# R语言基础课程 - 生物信息学分析脚本
# 作者：王运生
# 日期：2025-01-20
# 用途：Biostrings和GenomicRanges包的生物信息学应用

# 加载必需的包
suppressPackageStartupMessages({
  library(Biostrings)
  library(GenomicRanges)
  library(IRanges)
  library(dplyr)
  library(ggplot2)
})

# 清理环境
rm(list = ls())

cat("生物信息学R包应用练习\n")
cat("======================\n\n")

# 练习1: Biostrings包 - DNA序列分析
cat("练习1: Biostrings包 - DNA序列分析\n")
cat("----------------------------------\n")

# 创建DNA序列
dna_sequences <- DNAStringSet(c(
  BRCA1 = "ATGAAATTTCCCGGGAAATAGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
  TP53 = "ATGCCCGGGAAATTTAAATTTCCCGGGAAATAGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC",
  EGFR = "ATGAAATTTCCCGGGAAATAGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
  MYC = "ATGCCCAAATTTGGGCCCAAATAGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC",
  KRAS = "ATGAAACCCGGGTTTATATAGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
))

cat("DNA序列集合:\n")
print(dna_sequences)
cat("\n")

# 序列基本信息
cat("序列基本信息:\n")
cat("序列数量:", length(dna_sequences), "\n")
cat("序列名称:", names(dna_sequences), "\n")
cat("序列长度:\n")
seq_lengths <- width(dna_sequences)
for (i in 1:length(seq_lengths)) {
  cat(" ", names(dna_sequences)[i], ":", seq_lengths[i], "bp\n")
}
cat("\n")

# 序列操作
cat("序列操作示例:\n")

# 反向互补序列
cat("反向互补序列:\n")
rev_comp_sequences <- reverseComplement(dna_sequences)
for (i in 1:2) {  # 只显示前两个
  cat(" ", names(dna_sequences)[i], ":\n")
  cat("   原序列:", as.character(subseq(dna_sequences[i], 1, 50)), "...\n")
  cat("   反向互补:", as.character(subseq(rev_comp_sequences[i], 1, 50)), "...\n")
}
cat("\n")

# 翻译为蛋白质序列
cat("翻译为蛋白质序列:\n")
protein_sequences <- translate(dna_sequences)
for (i in 1:length(protein_sequences)) {
  cat(" ", names(protein_sequences)[i], ":", 
      as.character(subseq(protein_sequences[i], 1, 20)), "...\n")
}
cat("\n")

# 碱基频率分析
cat("碱基频率分析:\n")
base_freq <- alphabetFrequency(dna_sequences)
cat("各序列的碱基组成:\n")
print(base_freq[, 1:4])  # 只显示A, C, G, T
cat("\n")

# GC含量计算
cat("GC含量分析:\n")
gc_content <- letterFrequency(dna_sequences, "GC", as.prob = TRUE)
gc_df <- data.frame(
  Gene = names(dna_sequences),
  GC_Content = round(as.numeric(gc_content), 3),
  stringsAsFactors = FALSE
)
print(gc_df)
cat("\n")

# AT含量计算
at_content <- letterFrequency(dna_sequences, "AT", as.prob = TRUE)
gc_df$AT_Content <- round(as.numeric(at_content), 3)
gc_df$GC_AT_Ratio <- round(gc_df$GC_Content / gc_df$AT_Content, 3)

cat("GC/AT含量比较:\n")
print(gc_df)
cat("\n")

# 二核苷酸频率
cat("二核苷酸频率分析:\n")
dinuc_freq <- dinucleotideFrequency(dna_sequences)
cat("CpG二核苷酸频率:\n")
cpg_freq <- dinuc_freq[, "CG"]
for (i in 1:length(cpg_freq)) {
  cat(" ", names(dna_sequences)[i], ":", cpg_freq[i], "\n")
}
cat("\n")

# 模式匹配
cat("模式匹配分析:\n")

# 查找起始密码子ATG
cat("起始密码子ATG位置:\n")
atg_matches <- vmatchPattern("ATG", dna_sequences)
for (i in 1:length(atg_matches)) {
  positions <- start(atg_matches[[i]])
  cat(" ", names(dna_sequences)[i], ":", 
      ifelse(length(positions) > 0, paste(positions, collapse = ", "), "未找到"), "\n")
}
cat("\n")

# 计算特定模式的数量
patterns <- c("CG", "ATG", "TAA", "TAG", "TGA")
cat("特定模式计数:\n")
pattern_counts <- matrix(0, nrow = length(dna_sequences), ncol = length(patterns))
rownames(pattern_counts) <- names(dna_sequences)
colnames(pattern_counts) <- patterns

for (i in 1:length(patterns)) {
  pattern_counts[, i] <- vcountPattern(patterns[i], dna_sequences)
}

print(pattern_counts)
cat("\n")

# 练习2: GenomicRanges包 - 基因组区间操作
cat("练习2: GenomicRanges包 - 基因组区间操作\n")
cat("---------------------------------------\n")

# 创建基因信息数据
gene_info <- data.frame(
  gene_name = c("BRCA1", "BRCA2", "TP53", "EGFR", "MYC"),
  chromosome = c(17, 13, 17, 7, 8),
  start_position = c(43044295, 32315086, 7668402, 55019017, 127735434),
  end_position = c(43125483, 32400266, 7687538, 55207337, 127742951),
  strand = c("+", "+", "+", "+", "+"),
  expression = c(5.2, 3.8, 7.1, 4.5, 6.3),
  mutation_status = c("mutated", "wild_type", "mutated", "wild_type", "mutated"),
  stringsAsFactors = FALSE
)

# 创建GenomicRanges对象
gene_ranges <- GRanges(
  seqnames = paste0("chr", gene_info$chromosome),
  ranges = IRanges(
    start = gene_info$start_position,
    end = gene_info$end_position
  ),
  strand = gene_info$strand,
  gene_name = gene_info$gene_name,
  expression = gene_info$expression,
  mutation_status = gene_info$mutation_status
)

cat("基因组区间对象:\n")
print(gene_ranges)
cat("\n")

# 区间基本信息
cat("区间基本信息:\n")
cat("区间数量:", length(gene_ranges), "\n")
cat("染色体:", as.character(unique(seqnames(gene_ranges))), "\n")
cat("链方向:", as.character(unique(strand(gene_ranges))), "\n")

cat("基因长度 (bp):\n")
gene_widths <- width(gene_ranges)
for (i in 1:length(gene_widths)) {
  cat(" ", gene_ranges$gene_name[i], ":", gene_widths[i], "\n")
}
cat("\n")

# 区间操作
cat("区间操作示例:\n")

# 获取启动子区域 (上游2kb)
cat("启动子区域 (上游2kb):\n")
promoter_regions <- flank(gene_ranges, width = 2000, start = TRUE)
promoter_regions$region_type <- "promoter"
for (i in 1:length(promoter_regions)) {
  cat(" ", promoter_regions$gene_name[i], ":", 
      as.character(seqnames(promoter_regions[i])), ":",
      start(promoter_regions[i]), "-", end(promoter_regions[i]), "\n")
}
cat("\n")

# 获取基因下游区域 (下游1kb)
cat("下游区域 (下游1kb):\n")
downstream_regions <- flank(gene_ranges, width = 1000, start = FALSE)
downstream_regions$region_type <- "downstream"
for (i in 1:2) {  # 只显示前两个
  cat(" ", downstream_regions$gene_name[i], ":", 
      as.character(seqnames(downstream_regions[i])), ":",
      start(downstream_regions[i]), "-", end(downstream_regions[i]), "\n")
}
cat("\n")

# 调整基因区间大小
cat("调整基因区间大小:\n")
resized_genes <- resize(gene_ranges, width = 10000, fix = "center")
cat("调整为10kb (以中心为准):\n")
for (i in 1:3) {  # 只显示前三个
  cat(" ", resized_genes$gene_name[i], ":", 
      start(resized_genes[i]), "-", end(resized_genes[i]), 
      " (长度:", width(resized_genes[i]), ")\n")
}
cat("\n")

# 平移区间
cat("平移区间:\n")
shifted_genes <- shift(gene_ranges, 5000)
cat("向右平移5kb:\n")
for (i in 1:3) {
  original_start <- start(gene_ranges[i])
  shifted_start <- start(shifted_genes[i])
  cat(" ", shifted_genes$gene_name[i], ":", 
      original_start, "->", shifted_start, 
      " (差值:", shifted_start - original_start, ")\n")
}
cat("\n")

# 查找重叠
cat("重叠分析:\n")
overlaps <- findOverlaps(gene_ranges, gene_ranges)
cat("基因间重叠情况:\n")
overlap_df <- data.frame(
  Query = gene_ranges$gene_name[queryHits(overlaps)],
  Subject = gene_ranges$gene_name[subjectHits(overlaps)],
  stringsAsFactors = FALSE
)
# 移除自身重叠
overlap_df <- overlap_df[overlap_df$Query != overlap_df$Subject, ]
if (nrow(overlap_df) > 0) {
  print(overlap_df)
} else {
  cat("没有发现基因间重叠\n")
}
cat("\n")

# 按染色体分组
cat("按染色体分组:\n")
by_chr <- split(gene_ranges, seqnames(gene_ranges))
chr_summary <- data.frame(
  Chromosome = names(by_chr),
  Gene_Count = lengths(by_chr),
  stringsAsFactors = FALSE
)
print(chr_summary)
cat("\n")

# 练习3: 综合生物信息学分析
cat("练习3: 综合生物信息学分析\n")
cat("---------------------------\n")

# 整合序列和基因组位置信息
cat("整合分析:\n")

# 创建综合分析函数
comprehensive_gene_analysis <- function(sequences, ranges) {
  # 确保序列和区间的名称匹配
  common_genes <- intersect(names(sequences), ranges$gene_name)
  
  if (length(common_genes) == 0) {
    stop("序列和基因组区间没有共同的基因名称")
  }
  
  # 筛选匹配的数据
  matched_sequences <- sequences[common_genes]
  matched_ranges <- ranges[ranges$gene_name %in% common_genes]
  
  # 创建结果数据框
  results <- data.frame(
    gene_name = common_genes,
    chromosome = as.character(seqnames(matched_ranges)),
    start_pos = start(matched_ranges),
    end_pos = end(matched_ranges),
    genomic_length = width(matched_ranges),
    sequence_length = width(matched_sequences),
    stringsAsFactors = FALSE
  )
  
  # 添加序列分析结果
  results$gc_content <- round(as.numeric(letterFrequency(matched_sequences, "GC", as.prob = TRUE)), 3)
  results$at_content <- round(as.numeric(letterFrequency(matched_sequences, "AT", as.prob = TRUE)), 3)
  
  # 添加模式计数
  results$atg_count <- vcountPattern("ATG", matched_sequences)
  results$cg_count <- vcountPattern("CG", matched_sequences)
  results$stop_codon_count <- vcountPattern("TAA", matched_sequences) + 
                             vcountPattern("TAG", matched_sequences) + 
                             vcountPattern("TGA", matched_sequences)
  
  # 添加表达信息
  results$expression <- matched_ranges$expression
  results$mutation_status <- matched_ranges$mutation_status
  
  return(results)
}

# 执行综合分析
comprehensive_results <- comprehensive_gene_analysis(dna_sequences, gene_ranges)

cat("综合基因分析结果:\n")
print(comprehensive_results)
cat("\n")

# 生成分析报告
cat("分析报告:\n")
cat("=========\n")

cat("1. 基本统计:\n")
cat("   - 分析基因数量:", nrow(comprehensive_results), "\n")
cat("   - 涉及染色体:", length(unique(comprehensive_results$chromosome)), "个\n")
cat("   - 平均基因组长度:", round(mean(comprehensive_results$genomic_length), 0), "bp\n")
cat("   - 平均序列长度:", round(mean(comprehensive_results$sequence_length), 0), "bp\n")
cat("   - 平均GC含量:", round(mean(comprehensive_results$gc_content), 3), "\n")
cat("   - 平均表达水平:", round(mean(comprehensive_results$expression), 2), "\n\n")

cat("2. GC含量分析:\n")
high_gc_genes <- comprehensive_results$gene_name[comprehensive_results$gc_content > 0.5]
low_gc_genes <- comprehensive_results$gene_name[comprehensive_results$gc_content <= 0.5]
cat("   - 高GC含量基因 (>50%):", paste(high_gc_genes, collapse = ", "), "\n")
cat("   - 低GC含量基因 (≤50%):", paste(low_gc_genes, collapse = ", "), "\n\n")

cat("3. 表达水平分析:\n")
high_expr_genes <- comprehensive_results$gene_name[comprehensive_results$expression > 5.0]
low_expr_genes <- comprehensive_results$gene_name[comprehensive_results$expression <= 5.0]
cat("   - 高表达基因 (>5.0):", paste(high_expr_genes, collapse = ", "), "\n")
cat("   - 低表达基因 (≤5.0):", paste(low_expr_genes, collapse = ", "), "\n\n")

cat("4. 突变状态分析:\n")
mutation_summary <- table(comprehensive_results$mutation_status)
for (status in names(mutation_summary)) {
  genes_with_status <- comprehensive_results$gene_name[comprehensive_results$mutation_status == status]
  cat("   -", status, ":", mutation_summary[status], "个基因 -", 
      paste(genes_with_status, collapse = ", "), "\n")
}
cat("\n")

# 相关性分析
cat("5. 相关性分析:\n")
cor_gc_expr <- cor(comprehensive_results$gc_content, comprehensive_results$expression)
cor_length_expr <- cor(comprehensive_results$genomic_length, comprehensive_results$expression)
cor_gc_length <- cor(comprehensive_results$gc_content, comprehensive_results$genomic_length)

cat("   - GC含量与表达水平相关性:", round(cor_gc_expr, 3), "\n")
cat("   - 基因长度与表达水平相关性:", round(cor_length_expr, 3), "\n")
cat("   - GC含量与基因长度相关性:", round(cor_gc_length, 3), "\n\n")

# 创建可视化
cat("6. 创建可视化图表:\n")

# GC含量 vs 表达水平散点图
p1 <- ggplot(comprehensive_results, aes(x = gc_content, y = expression)) +
  geom_point(aes(color = mutation_status, size = genomic_length/1000), alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_color_manual(values = c("mutated" = "red", "wild_type" = "blue"),
                     labels = c("突变", "野生型")) +
  labs(title = "GC含量与表达水平的关系",
       x = "GC含量", y = "表达水平",
       color = "突变状态", size = "基因长度 (kb)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("gc_vs_expression.png", p1, width = 10, height = 6, dpi = 300)
print(p1)

# 基因特征热图
heatmap_data <- comprehensive_results %>%
  select(gene_name, gc_content, expression, genomic_length) %>%
  mutate(genomic_length = genomic_length / 1000) %>%  # 转换为kb
  column_to_rownames("gene_name") %>%
  scale() %>%  # 标准化
  as.data.frame() %>%
  rownames_to_column("gene_name") %>%
  reshape2::melt(id.vars = "gene_name")

p2 <- ggplot(heatmap_data, aes(x = variable, y = gene_name, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       name = "标准化值") +
  labs(title = "基因特征热图",
       x = "特征", y = "基因") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_x_discrete(labels = c("gc_content" = "GC含量", 
                              "expression" = "表达水平",
                              "genomic_length" = "基因长度"))

ggsave("gene_features_heatmap.png", p2, width = 8, height = 6, dpi = 300)
print(p2)

# 保存分析结果
write.csv(comprehensive_results, "comprehensive_gene_analysis.csv", row.names = FALSE)
cat("   - 综合分析结果已保存为 'comprehensive_gene_analysis.csv'\n")
cat("   - GC含量vs表达水平图已保存为 'gc_vs_expression.png'\n")
cat("   - 基因特征热图已保存为 'gene_features_heatmap.png'\n\n")

cat("生物信息学分析练习完成！\n")
cat("==========================\n")