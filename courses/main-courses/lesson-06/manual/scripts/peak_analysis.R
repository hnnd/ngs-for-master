# ChIP-seq Peak注释和功能分析脚本
# 作者：王运生
# 日期：2025-01-21
# 用法：Rscript peak_analysis.R [peak_file] [output_dir]

# 加载必要的包
suppressMessages({
    library(ChIPseeker)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(ggplot2)
    library(dplyr)
    library(VennDiagram)
})

# 参数设置
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    peak_file <- "results/peaks/H3K4me3_peaks.narrowPeak"
    output_dir <- "results/annotation"
} else {
    peak_file <- args[1]
    output_dir <- args[2]
}

# 创建输出目录
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

cat("ChIP-seq Peak注释和功能分析\n")
cat("Peak文件:", peak_file, "\n")
cat("输出目录:", output_dir, "\n")

# 检查输入文件
if (!file.exists(peak_file)) {
    stop("错误：Peak文件不存在 - ", peak_file)
}

# 1. 读取Peak文件
cat("\n=== 步骤1：读取Peak数据 ===\n")
peaks <- readPeakFile(peak_file)
cat("读取到", length(peaks), "个peaks\n")

# Peak基本统计
peak_widths <- width(peaks)
cat("Peak宽度统计：\n")
cat("  最小值:", min(peak_widths), "bp\n")
cat("  最大值:", max(peak_widths), "bp\n")
cat("  中位数:", median(peak_widths), "bp\n")
cat("  平均值:", round(mean(peak_widths)), "bp\n")

# 2. Peak基因组注释
cat("\n=== 步骤2：Peak基因组注释 ===\n")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# 执行注释
peakAnno <- annotatePeak(peaks, 
                        tssRegion = c(-2000, 2000),
                        TxDb = txdb,
                        annoDb = "org.Hs.eg.db")

cat("Peak注释完成\n")

# 保存注释结果
anno_df <- as.data.frame(peakAnno)
write.csv(anno_df, file.path(output_dir, "peak_annotation.csv"), row.names = FALSE)
cat("注释结果已保存到:", file.path(output_dir, "peak_annotation.csv"), "\n")

# 3. 基因组分布可视化
cat("\n=== 步骤3：基因组分布可视化 ===\n")

# 饼图
pdf("figures/peak_annotation_pie.pdf", width = 8, height = 6)
plotAnnoPie(peakAnno)
dev.off()

# 柱状图
pdf("figures/peak_annotation_bar.pdf", width = 10, height = 6)
plotAnnoBar(peakAnno)
dev.off()

# TSS距离分布
pdf("figures/peak_distance_to_TSS.pdf", width = 10, height = 6)
plotDistToTSS(peakAnno, title = "Distribution of ChIP peaks relative to TSS")
dev.off()

cat("基因组分布图已保存到 figures/ 目录\n")

# 4. 注释统计
cat("\n=== 步骤4：注释统计 ===\n")
anno_stats <- table(anno_df$annotation)
anno_percent <- round(anno_stats / sum(anno_stats) * 100, 2)

cat("Peak基因组分布：\n")
for (i in 1:length(anno_stats)) {
    cat(sprintf("  %s: %d (%.2f%%)\n", 
                names(anno_stats)[i], 
                anno_stats[i], 
                anno_percent[i]))
}

# 保存统计结果
stats_df <- data.frame(
    Annotation = names(anno_stats),
    Count = as.numeric(anno_stats),
    Percentage = as.numeric(anno_percent)
)
write.csv(stats_df, file.path(output_dir, "annotation_statistics.csv"), row.names = FALSE)

# 5. 获取Peak相关基因
cat("\n=== 步骤5：提取Peak相关基因 ===\n")
genes <- anno_df$geneId
genes <- genes[!is.na(genes)]
genes <- unique(genes)
cat("获得", length(genes), "个唯一基因\n")

# 保存基因列表
gene_symbols <- anno_df$SYMBOL[!is.na(anno_df$geneId)]
gene_symbols <- unique(gene_symbols[!is.na(gene_symbols)])

gene_list <- data.frame(
    GeneID = genes,
    Symbol = mapIds(org.Hs.eg.db, genes, "SYMBOL", "ENTREZID")
)
write.csv(gene_list, file.path(output_dir, "peak_associated_genes.csv"), row.names = FALSE)

# 6. GO富集分析
cat("\n=== 步骤6：GO富集分析 ===\n")
tryCatch({
    # 生物过程
    ego_bp <- enrichGO(gene = genes,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2,
                       readable = TRUE)
    
    if (nrow(as.data.frame(ego_bp)) > 0) {
        write.csv(as.data.frame(ego_bp), 
                  file.path(output_dir, "GO_BP_enrichment.csv"), 
                  row.names = FALSE)
        
        # 可视化
        pdf("figures/GO_BP_dotplot.pdf", width = 12, height = 8)
        print(dotplot(ego_bp, showCategory = 20) + 
              ggtitle("GO Biological Process Enrichment"))
        dev.off()
        
        pdf("figures/GO_BP_barplot.pdf", width = 12, height = 8)
        print(barplot(ego_bp, showCategory = 20) + 
              ggtitle("GO Biological Process Enrichment"))
        dev.off()
        
        cat("GO生物过程富集分析完成，发现", nrow(as.data.frame(ego_bp)), "个显著条目\n")
    } else {
        cat("GO生物过程富集分析无显著结果\n")
    }
    
    # 分子功能
    ego_mf <- enrichGO(gene = genes,
                       OrgDb = org.Hs.eg.db,
                       ont = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2,
                       readable = TRUE)
    
    if (nrow(as.data.frame(ego_mf)) > 0) {
        write.csv(as.data.frame(ego_mf), 
                  file.path(output_dir, "GO_MF_enrichment.csv"), 
                  row.names = FALSE)
        
        pdf("figures/GO_MF_dotplot.pdf", width = 12, height = 8)
        print(dotplot(ego_mf, showCategory = 20) + 
              ggtitle("GO Molecular Function Enrichment"))
        dev.off()
        
        cat("GO分子功能富集分析完成，发现", nrow(as.data.frame(ego_mf)), "个显著条目\n")
    } else {
        cat("GO分子功能富集分析无显著结果\n")
    }
    
}, error = function(e) {
    cat("GO富集分析出错:", e$message, "\n")
})

# 7. KEGG通路分析
cat("\n=== 步骤7：KEGG通路分析 ===\n")
tryCatch({
    kk <- enrichKEGG(gene = genes,
                     organism = 'hsa',
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2)
    
    if (nrow(as.data.frame(kk)) > 0) {
        write.csv(as.data.frame(kk), 
                  file.path(output_dir, "KEGG_enrichment.csv"), 
                  row.names = FALSE)
        
        pdf("figures/KEGG_dotplot.pdf", width = 12, height = 8)
        print(dotplot(kk, showCategory = 20) + 
              ggtitle("KEGG Pathway Enrichment"))
        dev.off()
        
        pdf("figures/KEGG_barplot.pdf", width = 12, height = 8)
        print(barplot(kk, showCategory = 20) + 
              ggtitle("KEGG Pathway Enrichment"))
        dev.off()
        
        cat("KEGG通路分析完成，发现", nrow(as.data.frame(kk)), "个显著通路\n")
    } else {
        cat("KEGG通路分析无显著结果\n")
    }
    
}, error = function(e) {
    cat("KEGG通路分析出错:", e$message, "\n")
})

# 8. Peak热图可视化
cat("\n=== 步骤8：Peak热图可视化 ===\n")
tryCatch({
    # 按染色体分组的Peak分布
    peak_chr <- table(seqnames(peaks))
    peak_chr_df <- data.frame(
        Chromosome = names(peak_chr),
        Count = as.numeric(peak_chr)
    )
    peak_chr_df <- peak_chr_df[order(peak_chr_df$Count, decreasing = TRUE), ]
    
    # 染色体分布图
    pdf("figures/peak_chromosome_distribution.pdf", width = 12, height = 6)
    p <- ggplot(peak_chr_df[1:22, ], aes(x = reorder(Chromosome, Count), y = Count)) +
        geom_bar(stat = "identity", fill = "steelblue") +
        coord_flip() +
        labs(title = "Peak Distribution Across Chromosomes",
             x = "Chromosome", y = "Number of Peaks") +
        theme_minimal()
    print(p)
    dev.off()
    
    cat("Peak染色体分布图已保存\n")
    
}, error = function(e) {
    cat("Peak可视化出错:", e$message, "\n")
})

# 9. 生成分析总结
cat("\n=== 步骤9：生成分析总结 ===\n")
summary_text <- paste0(
    "ChIP-seq Peak注释和功能分析总结\n",
    "=====================================\n\n",
    "分析时间: ", Sys.time(), "\n",
    "Peak文件: ", peak_file, "\n",
    "总Peak数: ", length(peaks), "\n",
    "相关基因数: ", length(genes), "\n\n",
    "Peak基因组分布:\n",
    paste(sprintf("  %s: %d (%.2f%%)", 
                  names(anno_stats), 
                  anno_stats, 
                  anno_percent), collapse = "\n"), "\n\n",
    "输出文件:\n",
    "  - Peak注释: ", file.path(output_dir, "peak_annotation.csv"), "\n",
    "  - 基因列表: ", file.path(output_dir, "peak_associated_genes.csv"), "\n",
    "  - GO富集: ", file.path(output_dir, "GO_*_enrichment.csv"), "\n",
    "  - KEGG富集: ", file.path(output_dir, "KEGG_enrichment.csv"), "\n",
    "  - 可视化图: figures/目录下的PDF文件\n\n",
    "分析建议:\n",
    "1. 检查Peak基因组分布是否符合预期\n",
    "2. 关注GO和KEGG富集的生物学意义\n",
    "3. 结合其他组学数据进行整合分析\n"
)

writeLines(summary_text, file.path(output_dir, "analysis_summary.txt"))
cat(summary_text)

cat("\n=== ChIP-seq Peak注释和功能分析完成 ===\n")
cat("所有结果已保存到:", output_dir, "和 figures/ 目录\n")