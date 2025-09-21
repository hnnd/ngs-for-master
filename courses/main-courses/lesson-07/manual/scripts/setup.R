# 单细胞RNA测序数据分析环境设置脚本
# 课程：高通量测序数据分析 - 第7次课
# 作者：王运生
# 日期：2025-01-21
# 用法：source("scripts/setup.R")

# 设置工作环境
cat("=== 单细胞RNA测序分析环境设置 ===\n")

# 检查并创建必要的目录
required_dirs <- c("data", "scripts", "results", "plots")
for (dir in required_dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    cat("创建目录:", dir, "\n")
  } else {
    cat("目录已存在:", dir, "\n")
  }
}

# 检查R版本
r_version <- R.version.string
cat("R版本:", r_version, "\n")

if (as.numeric(R.version$major) < 4) {
  warning("建议使用R 4.0.0或更高版本")
}

# 需要安装的包列表
required_packages <- c(
  "Seurat",
  "dplyr", 
  "ggplot2",
  "patchwork",
  "Matrix",
  "scales",
  "cowplot",
  "RColorBrewer",
  "viridis"
)

# 检查并安装缺失的包
cat("\n=== 检查R包依赖 ===\n")
missing_packages <- c()

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    missing_packages <- c(missing_packages, pkg)
    cat("缺失包:", pkg, "\n")
  } else {
    cat("已安装:", pkg, "\n")
  }
}

# 安装缺失的包
if (length(missing_packages) > 0) {
  cat("\n正在安装缺失的包...\n")
  
  # 首先尝试从CRAN安装
  cran_packages <- missing_packages[!missing_packages %in% c("Seurat")]
  if (length(cran_packages) > 0) {
    install.packages(cran_packages, dependencies = TRUE)
  }
  
  # 安装Seurat（可能需要特殊处理）
  if ("Seurat" %in% missing_packages) {
    cat("安装Seurat包...\n")
    install.packages("Seurat", dependencies = TRUE)
  }
  
  cat("包安装完成\n")
} else {
  cat("所有必需的包都已安装\n")
}

# 加载包并检查版本
cat("\n=== 加载R包并检查版本 ===\n")
package_versions <- list()

for (pkg in required_packages) {
  tryCatch({
    library(pkg, character.only = TRUE, quietly = TRUE)
    version <- as.character(packageVersion(pkg))
    package_versions[[pkg]] <- version
    cat(sprintf("%-12s: %s\n", pkg, version))
  }, error = function(e) {
    cat("加载", pkg, "失败:", e$message, "\n")
  })
}

# 检查关键包的版本要求
version_requirements <- list(
  "Seurat" = "4.0.0",
  "dplyr" = "1.0.0", 
  "ggplot2" = "3.3.0"
)

cat("\n=== 版本兼容性检查 ===\n")
for (pkg in names(version_requirements)) {
  if (pkg %in% names(package_versions)) {
    current_version <- package_versions[[pkg]]
    required_version <- version_requirements[[pkg]]
    
    if (compareVersion(current_version, required_version) >= 0) {
      cat("✓", pkg, "版本符合要求\n")
    } else {
      cat("⚠", pkg, "版本过低，建议升级到", required_version, "或更高版本\n")
    }
  }
}

# 设置全局选项
options(stringsAsFactors = FALSE)
options(dplyr.summarise.inform = FALSE)

# 设置ggplot2主题
theme_set(theme_bw() + 
          theme(panel.grid.minor = element_blank(),
                strip.background = element_blank(),
                text = element_text(size = 12)))

# 设置随机种子以确保结果可重现
set.seed(42)

# 检查系统内存
memory_info <- function() {
  if (.Platform$OS.type == "windows") {
    mem_limit <- memory.limit()
    cat("Windows内存限制:", mem_limit, "MB\n")
    if (mem_limit < 8000) {
      cat("⚠ 建议增加内存限制: memory.limit(16000)\n")
    }
  } else {
    # Unix/Linux系统
    cat("Unix/Linux系统 - 请确保有足够的内存（建议8GB+）\n")
  }
}

cat("\n=== 系统内存检查 ===\n")
memory_info()

# 创建分析配置
analysis_config <- list(
  project_name = "pbmc3k_analysis",
  min_cells = 3,
  min_features = 200,
  max_features = 2500,
  max_mt_percent = 20,
  n_variable_features = 2000,
  n_pcs = 15,
  clustering_resolution = 0.5,
  random_seed = 42
)

# 保存配置
saveRDS(analysis_config, "results/analysis_config.rds")
cat("\n分析配置已保存到 results/analysis_config.rds\n")

# 创建辅助函数
source_if_exists <- function(file) {
  if (file.exists(file)) {
    source(file)
    cat("已加载:", file, "\n")
  }
}

# 打印设置完成信息
cat("\n=== 环境设置完成 ===\n")
cat("工作目录:", getwd(), "\n")
cat("随机种子:", analysis_config$random_seed, "\n")
cat("准备开始单细胞数据分析...\n")

# 清理临时变量
rm(required_packages, missing_packages, cran_packages, 
   version_requirements, package_versions, required_dirs)

cat("\n环境设置脚本执行完成！\n")
cat("现在可以开始数据分析了。\n")