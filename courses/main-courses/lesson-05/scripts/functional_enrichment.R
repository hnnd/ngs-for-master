#!/usr/bin/env Rscript

# 功能富集分析脚本
# 课程：高通量测序数据分析 - 第5次课
# 作者：王运生
# 邮箱：wangys@hunau.edu.cn

# 加载必要的包
suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(DOSE)
    library(enrichplot)
    library(ggplot2)
    library(dplyr)
    library(stringr)
})

# 设置参数
WORK_DIR <- "~/ngs-analysis/lesson-05"
PADJ_CUTOFF <- 0.05
FC_CUTOFF <- 1.0

# 日志函数
log_info <- function(msg) {
    cat("[INFO]", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", msg, "\n")
}

log_error <- function(msg) {
    cat("[ERROR]", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", msg, "\n")
}

# 设置工作目录
setwd(WORK_DIR)

# 创建输出目录
dir.create("results/functional", showWarnings = FALSE, recursive = TRUE)
dir.create("results/plots/functional", showWarnings = FALSE, recursive = TRUE)

log_info("开始功能富集分析")

# 1. 读取差异表达结果
log_info("读取差异表达分析结果...")

sig_file <- "results/deseq2/significant_genes.csv"
if (!file.exists(sig_file)) {
    log_error(paste("显著差异基因文件不存在:", sig_file))
    quit(status = 1)
}

sig_genes <- read.csv(sig_file, stringsAsFactors = FALSE)
log_info(paste("读取到", nrow(sig_genes), "个显著差异基因"))

if (nrow(sig_genes) == 0) {
    log_error("没有显著差异基因，无法进行功能富集分析")
    quit(status = 1)
}

# 2. 基因ID转换
log_info("进行基因ID转换...")

# 提取Ensembl基因ID（去除版本号）
ensembl_ids <- str_replace(sig_genes$gene_id, "\\.\\d+$", "")

# 转换为Gene Symbol和ENTREZID
gene_mapping <- bitr(ensembl_ids, 
                     fromType = "ENSEMBL", 
                     toType = c("SYMBOL", "ENTREZID"), 
                     OrgDb = org.Hs.eg.db)

log_info(paste("成功转换", nrow(gene_mapping), "个基因ID"))

# 合并转换结果
sig_genes$ensembl_clean <- str_replace(sig_genes$gene_id, "\\.\\d+$", "")
sig_genes_mapped <- merge(sig_genes, gene_mapping, 
                         by.x = "ensembl_clean", by.y = "ENSEMBL", 
                         all.x = TRUE)

# 分离上调和下调基因
up_genes <- sig_genes_mapped[sig_genes_mapped$log2FoldChange > 0 & !is.na(sig_genes_mapped$SYMBOL), ]
down_genes <- sig_genes_mapped[sig_genes_mapped$log2FoldChange < 0 & !is.na(sig_genes_mapped$SYMBOL), ]

log_info(paste("上调基因:", nrow(up_genes)))
log_info(paste("下调基因:", nrow(down_genes)))

# 3. GO富集分析
log_info("进行GO富集分析...")

# GO富集分析 - 上调基因
if (nrow(up_genes) > 0) {
    log_info("分析上调基因的GO富集...")
    
    ego_up_bp <- enrichGO(gene = up_genes$SYMBOL,
                          OrgDb = org.Hs.eg.db,
                          ont = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.01,
                          qvalueCutoff = 0.05,
                          readable = TRUE)
    
    ego_up_mf <- enrichGO(gene = up_genes$SYMBOL,
                          OrgDb = org.Hs.eg.db,
                          ont = "MF",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.01,
                          qvalueCutoff = 0.05,
                          readable = TRUE)
    
    ego_up_cc <- enrichGO(gene = up_genes$SYMBOL,
                          OrgDb = org.Hs.eg.db,
                          ont = "CC",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.01,
                          qvalueCutoff = 0.05,
                          readable = TRUE)
    
    # 保存结果
    if (nrow(ego_up_bp) > 0) {
        write.csv(as.data.frame(ego_up_bp), "results/functional/GO_BP_upregulated.csv", row.names = FALSE)
        log_info(paste("上调基因BP富集:", nrow(ego_up_bp), "个terms"))
    }
    
    if (nrow(ego_up_mf) > 0) {
        write.csv(as.data.frame(ego_up_mf), "results/functional/GO_MF_upregulated.csv", row.names = FALSE)
        log_info(paste("上调基因MF富集:", nrow(ego_up_mf), "个terms"))
    }
    
    if (nrow(ego_up_cc) > 0) {
        write.csv(as.data.frame(ego_up_cc), "results/functional/GO_CC_upregulated.csv", row.names = FALSE)
        log_info(paste("上调基因CC富集:", nrow(ego_up_cc), "个terms"))
    }
}

# GO富集分析 - 下调基因
if (nrow(down_genes) > 0) {
    log_info("分析下调基因的GO富集...")
    
    ego_down_bp <- enrichGO(gene = down_genes$SYMBOL,
                            OrgDb = org.Hs.eg.db,
                            ont = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.01,
                            qvalueCutoff = 0.05,
                            readable = TRUE)
    
    ego_down_mf <- enrichGO(gene = down_genes$SYMBOL,
                            OrgDb = org.Hs.eg.db,
                            ont = "MF",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.01,
                            qvalueCutoff = 0.05,
                            readable = TRUE)
    
    ego_down_cc <- enrichGO(gene = down_genes$SYMBOL,
                            OrgDb = org.Hs.eg.db,
                            ont = "CC",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.01,
                            qvalueCutoff = 0.05,
                            readable = TRUE)
    
    # 保存结果
    if (nrow(ego_down_bp) > 0) {
        write.csv(as.data.frame(ego_down_bp), "results/functional/GO_BP_downregulated.csv", row.names = FALSE)
        log_info(paste("下调基因BP富集:", nrow(ego_down_bp), "个terms"))
    }
    
    if (nrow(ego_down_mf) > 0) {
        write.csv(as.data.frame(ego_down_mf), "results/functional/GO_MF_downregulated.csv", row.names = FALSE)
        log_info(paste("下调基因MF富集:", nrow(ego_down_mf), "个terms"))
    }
    
    if (nrow(ego_down_cc) > 0) {
        write.csv(as.data.frame(ego_down_cc), "results/functional/GO_CC_downregulated.csv", row.names = FALSE)
        log_info(paste("下调基因CC富集:", nrow(ego_down_cc), "个terms"))
    }
}

# 4. KEGG通路分析
log_info("进行KEGG通路分析...")

# KEGG分析 - 上调基因
if (nrow(up_genes) > 0 && sum(!is.na(up_genes$ENTREZID)) > 0) {
    log_info("分析上调基因的KEGG通路...")
    
    kk_up <- enrichKEGG(gene = up_genes$ENTREZID[!is.na(up_genes$ENTREZID)],
                        organism = 'hsa',
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH")
    
    if (nrow(kk_up) > 0) {
        write.csv(as.data.frame(kk_up), "results/functional/KEGG_upregulated.csv", row.names = FALSE)
        log_info(paste("上调基因KEGG通路:", nrow(kk_up), "个通路"))
    }
}

# KEGG分析 - 下调基因
if (nrow(down_genes) > 0 && sum(!is.na(down_genes$ENTREZID)) > 0) {
    log_info("分析下调基因的KEGG通路...")
    
    kk_down <- enrichKEGG(gene = down_genes$ENTREZID[!is.na(down_genes$ENTREZID)],
                          organism = 'hsa',
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH")
    
    if (nrow(kk_down) > 0) {
        write.csv(as.data.frame(kk_down), "results/functional/KEGG_downregulated.csv", row.names = FALSE)
        log_info(paste("下调基因KEGG通路:", nrow(kk_down), "个通路"))
    }
}

# 5. 生成可视化图表
log_info("生成功能富集可视化图表...")

# GO富集可视化 - 上调基因
if (exists("ego_up_bp") && nrow(ego_up_bp) > 0) {
    # 点图
    p1 <- dotplot(ego_up_bp, showCategory = 20) + 
        ggtitle("GO BP Enrichment - Upregulated Genes")
    ggsave("results/plots/functional/GO_BP_upregulated_dotplot.pdf", p1, width = 12, height = 8)
    
    # 条形图
    p2 <- barplot(ego_up_bp, showCategory = 15) + 
        ggtitle("GO BP Enrichment - Upregulated Genes")
    ggsave("results/plots/functional/GO_BP_upregulated_barplot.pdf", p2, width = 10, height = 8)
    
    # 网络图（如果terms数量适中）
    if (nrow(ego_up_bp) >= 5 && nrow(ego_up_bp) <= 50) {
        p3 <- emapplot(pairwise_termsim(ego_up_bp), showCategory = 20) + 
            ggtitle("GO BP Network - Upregulated Genes")
        ggsave("results/plots/functional/GO_BP_upregulated_network.pdf", p3, width = 12, height = 10)
    }
}

# GO富集可视化 - 下调基因
if (exists("ego_down_bp") && nrow(ego_down_bp) > 0) {
    # 点图
    p4 <- dotplot(ego_down_bp, showCategory = 20) + 
        ggtitle("GO BP Enrichment - Downregulated Genes")
    ggsave("results/plots/functional/GO_BP_downregulated_dotplot.pdf", p4, width = 12, height = 8)
    
    # 条形图
    p5 <- barplot(ego_down_bp, showCategory = 15) + 
        ggtitle("GO BP Enrichment - Downregulated Genes")
    ggsave("results/plots/functional/GO_BP_downregulated_barplot.pdf", p5, width = 10, height = 8)
    
    # 网络图
    if (nrow(ego_down_bp) >= 5 && nrow(ego_down_bp) <= 50) {
        p6 <- emapplot(pairwise_termsim(ego_down_bp), showCategory = 20) + 
            ggtitle("GO BP Network - Downregulated Genes")
        ggsave("results/plots/functional/GO_BP_downregulated_network.pdf", p6, width = 12, height = 10)
    }
}

# KEGG通路可视化
if (exists("kk_up") && nrow(kk_up) > 0) {
    p7 <- dotplot(kk_up, showCategory = 15) + 
        ggtitle("KEGG Pathway - Upregulated Genes")
    ggsave("results/plots/functional/KEGG_upregulated_dotplot.pdf", p7, width = 12, height = 8)
}

if (exists("kk_down") && nrow(kk_down) > 0) {
    p8 <- dotplot(kk_down, showCategory = 15) + 
        ggtitle("KEGG Pathway - Downregulated Genes")
    ggsave("results/plots/functional/KEGG_downregulated_dotplot.pdf", p8, width = 12, height = 8)
}

# 6. 比较分析
log_info("进行比较分析...")

# 如果上调和下调基因都有富集结果，进行比较
if (exists("ego_up_bp") && exists("ego_down_bp") && 
    nrow(ego_up_bp) > 0 && nrow(ego_down_bp) > 0) {
    
    # 创建基因列表用于比较
    gene_list <- list(
        "Upregulated" = up_genes$SYMBOL[!is.na(up_genes$SYMBOL)],
        "Downregulated" = down_genes$SYMBOL[!is.na(down_genes$SYMBOL)]
    )
    
    # 比较GO富集
    compare_go <- compareCluster(geneClusters = gene_list,
                                fun = "enrichGO",
                                OrgDb = org.Hs.eg.db,
                                ont = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.01,
                                qvalueCutoff = 0.05)
    
    if (nrow(compare_go) > 0) {
        # 保存比较结果
        write.csv(as.data.frame(compare_go), "results/functional/GO_comparison.csv", row.names = FALSE)
        
        # 生成比较图
        p9 <- dotplot(compare_go, showCategory = 10) + 
            ggtitle("GO BP Enrichment Comparison")
        ggsave("results/plots/functional/GO_comparison_dotplot.pdf", p9, width = 14, height = 10)
    }
}

# 7. 生成功能分析摘要
log_info("生成功能分析摘要...")

summary_file <- "results/functional/functional_analysis_summary.txt"
cat("功能富集分析摘要\n", file = summary_file)
cat("================\n\n", file = summary_file, append = TRUE)
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", file = summary_file, append = TRUE)

cat("基因数量统计:\n", file = summary_file, append = TRUE)
cat("- 总差异基因:", nrow(sig_genes), "\n", file = summary_file, append = TRUE)
cat("- 成功转换ID的基因:", nrow(gene_mapping), "\n", file = summary_file, append = TRUE)
cat("- 上调基因:", nrow(up_genes), "\n", file = summary_file, append = TRUE)
cat("- 下调基因:", nrow(down_genes), "\n", file = summary_file, append = TRUE)

cat("\nGO富集结果:\n", file = summary_file, append = TRUE)
if (exists("ego_up_bp")) {
    cat("- 上调基因BP terms:", nrow(ego_up_bp), "\n", file = summary_file, append = TRUE)
}
if (exists("ego_up_mf")) {
    cat("- 上调基因MF terms:", nrow(ego_up_mf), "\n", file = summary_file, append = TRUE)
}
if (exists("ego_up_cc")) {
    cat("- 上调基因CC terms:", nrow(ego_up_cc), "\n", file = summary_file, append = TRUE)
}
if (exists("ego_down_bp")) {
    cat("- 下调基因BP terms:", nrow(ego_down_bp), "\n", file = summary_file, append = TRUE)
}
if (exists("ego_down_mf")) {
    cat("- 下调基因MF terms:", nrow(ego_down_mf), "\n", file = summary_file, append = TRUE)
}
if (exists("ego_down_cc")) {
    cat("- 下调基因CC terms:", nrow(ego_down_cc), "\n", file = summary_file, append = TRUE)
}

cat("\nKEGG通路结果:\n", file = summary_file, append = TRUE)
if (exists("kk_up")) {
    cat("- 上调基因通路:", nrow(kk_up), "\n", file = summary_file, append = TRUE)
}
if (exists("kk_down")) {
    cat("- 下调基因通路:", nrow(kk_down), "\n", file = summary_file, append = TRUE)
}

# 8. 保存R对象
log_info("保存R对象...")
objects_to_save <- c("sig_genes_mapped", "up_genes", "down_genes", "gene_mapping")

# 添加存在的富集结果对象
if (exists("ego_up_bp")) objects_to_save <- c(objects_to_save, "ego_up_bp")
if (exists("ego_up_mf")) objects_to_save <- c(objects_to_save, "ego_up_mf")
if (exists("ego_up_cc")) objects_to_save <- c(objects_to_save, "ego_up_cc")
if (exists("ego_down_bp")) objects_to_save <- c(objects_to_save, "ego_down_bp")
if (exists("ego_down_mf")) objects_to_save <- c(objects_to_save, "ego_down_mf")
if (exists("ego_down_cc")) objects_to_save <- c(objects_to_save, "ego_down_cc")
if (exists("kk_up")) objects_to_save <- c(objects_to_save, "kk_up")
if (exists("kk_down")) objects_to_save <- c(objects_to_save, "kk_down")
if (exists("compare_go")) objects_to_save <- c(objects_to_save, "compare_go")

save(list = objects_to_save, file = "results/functional/functional_objects.RData")

log_info("功能富集分析完成！")
log_info("结果文件保存在 results/functional/ 目录")
log_info("图表文件保存在 results/plots/functional/ 目录")

# 列出生成的文件
log_info("生成的主要文件:")
functional_files <- list.files("results/functional", full.names = TRUE)
plot_files <- list.files("results/plots/functional", full.names = TRUE)

for (file in c(functional_files, plot_files)) {
    log_info(paste("✓", file))
}