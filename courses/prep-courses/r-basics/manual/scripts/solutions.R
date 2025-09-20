#!/usr/bin/env Rscript
# R语言基础课程 - 练习答案脚本
# 作者：王运生
# 日期：2025-01-20
# 用途：课程练习题的参考答案

# 加载必需的包
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(Biostrings)
  library(GenomicRanges)
})

# 清理环境
rm(list = ls())

cat("R语言基础课程练习答案\n")
cat("======================\n\n")

# 练习题1: 基础数据类型和向量操作
cat("练习题1: 基础数据类型和向量操作\n")
cat("-------------------------------\n")

# 题目：创建一个包含5个基因名称的向量，计算每个基因名称的字符数，
# 并找出字符数最多的基因名称

# 答案：
gene_names <- c("BRCA1", "TP53", "EGFR", "MYC", "KRAS")
char_counts <- nchar(gene_names)
names(char_counts) <- gene_names

cat("基因名称及其字符数:\n")
print(char_counts)

max_char_gene <- names(char_counts)[which.max(char_counts)]
cat("字符数最多的基因:", max_char_gene, "(", max(char_counts), "个字符)\n\n")

# 练习题2: 数据框操作和筛选
cat("练习题2: 数据框操作和筛选\n")
cat("-------------------------\n")

# 题目：创建一个基因表达数据框，包含基因名、染色体、表达水平和突变状态，
# 筛选出17号染色体上表达水平大于5.0的基因

# 答案：
gene_data <- data.frame(
  gene_name = c("BRCA1", "BRCA2", "TP53", "EGFR", "MYC", "PTEN", "RB1"),
  chromosome = c(17, 13, 17, 7, 8, 10, 13),
  expression = c(5.2, 3.8, 7.1, 4.5, 6.3, 4.2, 5.8),
  mutation_status = c("mutated", "wild_type", "mutated", "wild_type", 
                      "mutated", "wild_type", "mutated"),
  stringsAsFactors = FALSE
)

cat("原始基因数据:\n")
print(gene_data)

# 筛选条件
filtered_genes <- gene_data %>%
  filter(chromosome == 17 & expression > 5.0)

cat("\n17号染色体上表达水平>5.0的基因:\n")
print(filtered_genes)
cat("\n")

# 练习题3: 字符串处理和序列分析
cat("练习题3: 字符串处理和序列分析\n")
cat("-----------------------------\n")

# 题目：编写函数计算DNA序列的GC含量，并分析一组序列的GC含量分布

# 答案：
calculate_gc_content <- function(sequence) {
  # 转换为大写
  seq_upper <- toupper(sequence)
  
  # 计算总长度
  total_length <- nchar(seq_upper)
  
  if (total_length == 0) {
    return(0)
  }
  
  # 计算G和C的数量
  g_count <- lengths(regmatches(seq_upper, gregexpr("G", seq_upper)))
  c_count <- lengths(regmatches(seq_upper, gregexpr("C", seq_upper)))
  
  # 计算GC含量
  gc_content <- (g_count + c_count) / total_length
  
  return(round(gc_content, 3))
}

# 测试序列
test_sequences <- c(
  "ATCGATCGATCG",
  "GCGCGCGCGCGC", 
  "ATATATATATATAT",
  "CGATCGATCGATCG",
  "TTTTTTTTTTTTTT"
)

cat("序列GC含量分析:\n")
for (i in 1:length(test_sequences)) {
  seq <- test_sequences[i]
  gc <- calculate_gc_content(seq)
  cat("序列", i, ":", seq, "- GC含量:", gc, "\n")
}

# 统计分析
gc_values <- sapply(test_sequences, calculate_gc_content)
cat("\nGC含量统计:\n")
cat("平均值:", round(mean(gc_values), 3), "\n")
cat("标准差:", round(sd(gc_values), 3), "\n")
cat("范围:", range(gc_values), "\n\n")

# 练习题4: 统计分析和假设检验
cat("练习题4: 统计分析和假设检验\n")
cat("---------------------------\n")

# 题目：比较突变基因和野生型基因的表达水平是否有显著差异

# 答案：
mutated_expr <- gene_data$expression[gene_data$mutation_status == "mutated"]
wildtype_expr <- gene_data$expression[gene_data$mutation_status == "wild_type"]

cat("突变基因表达水平:", mutated_expr, "\n")
cat("野生型基因表达水平:", wildtype_expr, "\n\n")

# 描述性统计
cat("描述性统计:\n")
cat("突变基因 - 平均值:", round(mean(mutated_expr), 2), 
    ", 标准差:", round(sd(mutated_expr), 2), "\n")
cat("野生型基因 - 平均值:", round(mean(wildtype_expr), 2), 
    ", 标准差:", round(sd(wildtype_expr), 2), "\n\n")

# t检验
t_test_result <- t.test(mutated_expr, wildtype_expr)
cat("双样本t检验结果:\n")
cat("t统计量:", round(t_test_result$statistic, 4), "\n")
cat("p值:", round(t_test_result$p.value, 4), "\n")
cat("95%置信区间: [", round(t_test_result$conf.int[1], 3), 
    ", ", round(t_test_result$conf.int[2], 3), "]\n")

# 结论
if (t_test_result$p.value < 0.05) {
  cat("结论: 拒绝原假设，突变基因和野生型基因的表达水平有显著差异 (p < 0.05)\n")
} else {
  cat("结论: 不能拒绝原假设，突变基因和野生型基因的表达水平无显著差异 (p ≥ 0.05)\n")
}
cat("\n")

# 练习题5: 数据可视化
cat("练习题5: 数据可视化\n")
cat("-------------------\n")

# 题目：创建基因表达水平的箱线图，按突变状态分组

# 答案：
p1 <- ggplot(gene_data, aes(x = mutation_status, y = expression)) +
  geom_boxplot(aes(fill = mutation_status), alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2) +
  scale_fill_manual(values = c("mutated" = "red", "wild_type" = "blue"),
                    labels = c("突变", "野生型")) +
  scale_x_discrete(labels = c("突变", "野生型")) +
  labs(title = "基因表达水平按突变状态分组",
       x = "突变状态", y = "表达水平",
       fill = "突变状态") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  ) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
               fill = "yellow", color = "black")

ggsave("expression_by_mutation_boxplot.png", p1, width = 8, height = 6, dpi = 300)
print(p1)

cat("箱线图已保存为 'expression_by_mutation_boxplot.png'\n\n")

# 练习题6: dplyr数据处理
cat("练习题6: dplyr数据处理\n")
cat("---------------------\n")

# 题目：使用dplyr计算每个染色体上基因的平均表达水平，并按平均表达水平排序

# 答案：
chr_summary <- gene_data %>%
  group_by(chromosome) %>%
  summarise(
    gene_count = n(),
    mean_expression = round(mean(expression), 2),
    sd_expression = round(sd(expression), 2),
    min_expression = min(expression),
    max_expression = max(expression),
    .groups = 'drop'
  ) %>%
  arrange(desc(mean_expression))

cat("按染色体分组的表达统计 (按平均表达水平降序排列):\n")
print(chr_summary)
cat("\n")

# 练习题7: 函数编写
cat("练习题7: 函数编写\n")
cat("-----------------\n")

# 题目：编写一个函数，输入基因表达数据框，返回表达水平分类的汇总统计

# 答案：
classify_and_summarize_expression <- function(data, high_cutoff = 6.0, low_cutoff = 3.0) {
  # 添加表达分类
  classified_data <- data %>%
    mutate(
      expression_category = case_when(
        expression > high_cutoff ~ "高表达",
        expression > low_cutoff ~ "中等表达",
        TRUE ~ "低表达"
      )
    )
  
  # 汇总统计
  summary_stats <- classified_data %>%
    group_by(expression_category) %>%
    summarise(
      count = n(),
      mean_expression = round(mean(expression), 2),
      genes = paste(gene_name, collapse = ", "),
      .groups = 'drop'
    ) %>%
    arrange(desc(mean_expression))
  
  # 返回结果列表
  return(list(
    classified_data = classified_data,
    summary = summary_stats
  ))
}

# 测试函数
result <- classify_and_summarize_expression(gene_data)

cat("表达水平分类结果:\n")
print(result$classified_data)

cat("\n表达分类汇总:\n")
print(result$summary)
cat("\n")

# 练习题8: Biostrings应用
cat("练习题8: Biostrings应用\n")
cat("-----------------------\n")

# 题目：使用Biostrings分析DNA序列，找出所有可能的开放阅读框(ORF)

# 答案：
find_orfs <- function(dna_sequence, min_length = 30) {
  # 转换为DNAString对象
  if (is.character(dna_sequence)) {
    dna_seq <- DNAString(dna_sequence)
  } else {
    dna_seq <- dna_sequence
  }
  
  # 起始和终止密码子
  start_codon <- "ATG"
  stop_codons <- c("TAA", "TAG", "TGA")
  
  orfs <- list()
  
  # 检查三个阅读框
  for (frame in 0:2) {
    seq_frame <- subseq(dna_seq, start = frame + 1)
    
    # 查找起始密码子
    start_positions <- start(matchPattern(start_codon, seq_frame))
    
    for (start_pos in start_positions) {
      # 从起始位置开始查找终止密码子
      subseq_from_start <- subseq(seq_frame, start = start_pos + 3)
      
      for (stop_codon in stop_codons) {
        stop_matches <- matchPattern(stop_codon, subseq_from_start)
        
        if (length(stop_matches) > 0) {
          # 取第一个终止密码子
          first_stop <- start(stop_matches)[1]
          orf_length <- first_stop + 2  # 包含终止密码子
          
          if (orf_length >= min_length) {
            orf_seq <- subseq(seq_frame, start = start_pos, width = orf_length + 3)
            
            orfs <- append(orfs, list(list(
              frame = frame + 1,
              start = start_pos + frame,
              end = start_pos + frame + orf_length + 2,
              length = orf_length + 3,
              sequence = as.character(orf_seq)
            )))
          }
          break  # 找到第一个终止密码子就停止
        }
      }
    }
  }
  
  return(orfs)
}

# 测试序列
test_dna <- "ATGAAATTTCCCGGGAAATAGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGTAG"

cat("测试DNA序列:", test_dna, "\n")
cat("序列长度:", nchar(test_dna), "bp\n\n")

# 查找ORF
orfs <- find_orfs(test_dna, min_length = 30)

cat("找到的开放阅读框:\n")
if (length(orfs) > 0) {
  for (i in 1:length(orfs)) {
    orf <- orfs[[i]]
    cat("ORF", i, ":\n")
    cat("  阅读框:", orf$frame, "\n")
    cat("  位置:", orf$start, "-", orf$end, "\n")
    cat("  长度:", orf$length, "bp\n")
    cat("  序列:", substr(orf$sequence, 1, 50), "...\n\n")
  }
} else {
  cat("未找到符合条件的ORF\n")
}

# 练习题9: GenomicRanges应用
cat("练习题9: GenomicRanges应用\n")
cat("-------------------------\n")

# 题目：创建基因组区间对象，计算基因密度（每Mb基因数）

# 答案：
# 创建GenomicRanges对象
gene_ranges <- GRanges(
  seqnames = paste0("chr", gene_data$chromosome),
  ranges = IRanges(
    start = c(43044295, 32315086, 7668402, 55019017, 127735434, 89623194, 48877297),
    end = c(43125483, 32400266, 7687538, 55207337, 127742951, 89728532, 48981047)
  ),
  strand = "+",
  gene_name = gene_data$gene_name,
  expression = gene_data$expression
)

cat("基因组区间对象:\n")
print(gene_ranges)
cat("\n")

# 计算每个染色体的基因密度
chr_density <- gene_ranges %>%
  as.data.frame() %>%
  group_by(seqnames) %>%
  summarise(
    gene_count = n(),
    total_span = max(end) - min(start) + 1,
    span_mb = round(total_span / 1e6, 2),
    density_per_mb = round(gene_count / (total_span / 1e6), 2),
    .groups = 'drop'
  ) %>%
  arrange(desc(density_per_mb))

cat("染色体基因密度分析:\n")
print(chr_density)
cat("\n")

# 练习题10: 综合分析项目
cat("练习题10: 综合分析项目\n")
cat("----------------------\n")

# 题目：整合所有分析，创建一个基因特征的综合报告

# 答案：
create_gene_report <- function(gene_data, sequences = NULL) {
  report <- list()
  
  # 基本统计
  report$basic_stats <- list(
    total_genes = nrow(gene_data),
    unique_chromosomes = length(unique(gene_data$chromosome)),
    mean_expression = round(mean(gene_data$expression), 2),
    expression_range = range(gene_data$expression)
  )
  
  # 突变状态分析
  mutation_table <- table(gene_data$mutation_status)
  report$mutation_analysis <- list(
    mutation_counts = mutation_table,
    mutation_percentages = round(prop.table(mutation_table) * 100, 1)
  )
  
  # 表达水平分析
  high_expr <- gene_data$gene_name[gene_data$expression > 5.0]
  low_expr <- gene_data$gene_name[gene_data$expression <= 5.0]
  
  report$expression_analysis <- list(
    high_expression_genes = high_expr,
    low_expression_genes = low_expr,
    high_expr_count = length(high_expr),
    low_expr_count = length(low_expr)
  )
  
  # 染色体分布
  chr_dist <- table(gene_data$chromosome)
  report$chromosome_distribution <- chr_dist
  
  # 序列分析（如果提供）
  if (!is.null(sequences)) {
    gc_contents <- sapply(sequences, calculate_gc_content)
    report$sequence_analysis <- list(
      mean_gc_content = round(mean(gc_contents), 3),
      gc_range = range(gc_contents),
      high_gc_genes = names(gc_contents)[gc_contents > 0.5]
    )
  }
  
  return(report)
}

# 创建测试序列
test_sequences <- c(
  BRCA1 = "ATGAAATTTCCCGGGAAATAGATCGATCGATCG",
  BRCA2 = "ATGCCCGGGAAATTTAAATTTCCCGGGAAATAG", 
  TP53 = "ATGCCCGGGAAATTTAAATTTCCCGGGAAATAG",
  EGFR = "ATGAAATTTCCCGGGAAATAGATCGATCGATCG",
  MYC = "ATGCCCAAATTTGGGCCCAAATAGATCGATCG",
  PTEN = "ATGAAACCCGGGTTTATATAGATCGATCGATC",
  RB1 = "ATGCCCGGGAAATTTAAATTTCCCGGGAAATAG"
)

# 生成报告
comprehensive_report <- create_gene_report(gene_data, test_sequences)

cat("基因特征综合报告:\n")
cat("==================\n\n")

cat("1. 基本统计信息:\n")
cat("   总基因数:", comprehensive_report$basic_stats$total_genes, "\n")
cat("   涉及染色体数:", comprehensive_report$basic_stats$unique_chromosomes, "\n")
cat("   平均表达水平:", comprehensive_report$basic_stats$mean_expression, "\n")
cat("   表达水平范围:", comprehensive_report$basic_stats$expression_range, "\n\n")

cat("2. 突变状态分析:\n")
for (status in names(comprehensive_report$mutation_analysis$mutation_counts)) {
  count <- comprehensive_report$mutation_analysis$mutation_counts[status]
  percent <- comprehensive_report$mutation_analysis$mutation_percentages[status]
  cat("  ", status, ":", count, "个 (", percent, "%)\n")
}
cat("\n")

cat("3. 表达水平分析:\n")
cat("   高表达基因 (>5.0):", comprehensive_report$expression_analysis$high_expr_count, "个\n")
cat("   ", paste(comprehensive_report$expression_analysis$high_expression_genes, collapse = ", "), "\n")
cat("   低表达基因 (≤5.0):", comprehensive_report$expression_analysis$low_expr_count, "个\n")
cat("   ", paste(comprehensive_report$expression_analysis$low_expression_genes, collapse = ", "), "\n\n")

cat("4. 染色体分布:\n")
for (chr in names(comprehensive_report$chromosome_distribution)) {
  count <- comprehensive_report$chromosome_distribution[chr]
  cat("   染色体", chr, ":", count, "个基因\n")
}
cat("\n")

cat("5. 序列分析:\n")
cat("   平均GC含量:", comprehensive_report$sequence_analysis$mean_gc_content, "\n")
cat("   GC含量范围:", comprehensive_report$sequence_analysis$gc_range, "\n")
cat("   高GC含量基因:", paste(comprehensive_report$sequence_analysis$high_gc_genes, collapse = ", "), "\n\n")

# 保存报告
report_text <- capture.output({
  cat("基因特征综合报告\n")
  cat("==================\n\n")
  cat("生成时间:", Sys.time(), "\n")
  cat("分析基因数:", nrow(gene_data), "\n\n")
  print(comprehensive_report)
})

writeLines(report_text, "gene_analysis_report.txt")
cat("详细报告已保存为 'gene_analysis_report.txt'\n\n")

cat("所有练习题答案完成！\n")
cat("======================\n")