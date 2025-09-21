# 第8次课实践操作手册

## 课程信息
- **课程名称**：高通量测序数据分析
- **主讲教师**：王运生
- **联系邮箱**：wangys@hunau.edu.cn
- **办公室**：16教420室
- **上课地点**：105机房
- **课程时间**：第8次课（2学时实践）

## 实验目标

### 主要目标
- 掌握多组学数据的预处理和整合方法
- 学会使用R和Python进行数据整合分析
- 实践机器学习算法在生物数据中的应用
- 体验云计算平台进行大规模数据分析

### 预期成果
- 完成多组学数据整合分析流程
- 构建机器学习预测模型
- 生成数据可视化图表和分析报告
- 掌握云计算平台的基本使用方法

## 环境要求

### 软件环境
| 软件名称 | 版本要求 | 安装方式 | 说明 |
|---------|---------|---------|------|
| R | >= 4.0.0 | 官网下载 | 统计分析环境 |
| RStudio | >= 1.4.0 | 官网下载 | R集成开发环境 |
| Python | >= 3.8 | Anaconda | Python环境 |
| Jupyter | 最新版 | pip install jupyter | 交互式开发 |

### R包依赖
```r
# 安装必需的R包
install.packages(c("mixOmics", "MultiAssayExperiment", "caret", 
                   "randomForest", "ggplot2", "pheatmap", "VennDiagram"))

# Bioconductor包
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR", "DESeq2"))
```

### Python包依赖
```bash
# 安装Python包
pip install pandas numpy scikit-learn matplotlib seaborn
pip install tensorflow keras plotly
pip install umap-learn
```

### 硬件要求
- **内存**：至少 8 GB RAM
- **存储空间**：至少 5 GB 可用空间
- **CPU**：多核处理器推荐
- **网络**：稳定的互联网连接（用于云计算体验）

### 数据准备
| 数据文件 | 大小 | 下载链接/位置 | 说明 |
|---------|------|-------------|------|
| genomics_data.csv | 2MB | scripts/data/ | 基因组变异数据 |
| transcriptomics_data.csv | 5MB | scripts/data/ | 转录组表达数据 |
| proteomics_data.csv | 3MB | scripts/data/ | 蛋白质组数据 |
| metadata.csv | 50KB | scripts/data/ | 样本信息 |

## 操作步骤

### 步骤1：环境设置和准备工作

#### 1.1 创建工作目录
```bash
# 创建本次实验的工作目录
mkdir -p ~/ngs-analysis/lesson-08
cd ~/ngs-analysis/lesson-08

# 创建子目录结构
mkdir -p {data,scripts,results,plots,models}
```

#### 1.2 检查软件环境
```bash
# 检查R环境
R --version

# 检查Python环境
python --version

# 检查Jupyter
jupyter --version
```

**预期输出：**
```
R version 4.0.0 (2020-04-24) -- "Arbor Day"
Python 3.8.5
jupyter core     : 4.7.1
```

#### 1.3 下载和准备数据
```bash
# 复制实验数据到工作目录
cp /path/to/course/data/*.csv data/

# 验证数据完整性
ls -la data/
wc -l data/*.csv
```

**检查点：** 确认所有数据文件已正确下载并位于 `data/` 目录中。

---

### 步骤2：多组学数据预处理

#### 2.1 数据质量评估

**操作说明：**
首先对各组学数据进行质量评估，了解数据分布、缺失值情况和异常值。

**执行命令：**
```bash
# 启动R环境
R
```

**R代码：**
```r
# 加载必要的包
library(ggplot2)
library(pheatmap)
library(VennDiagram)

# 读取数据
genomics <- read.csv("data/genomics_data.csv", row.names = 1)
transcriptomics <- read.csv("data/transcriptomics_data.csv", row.names = 1)
proteomics <- read.csv("data/proteomics_data.csv", row.names = 1)
metadata <- read.csv("data/metadata.csv", row.names = 1)

# 查看数据基本信息
cat("基因组数据维度:", dim(genomics), "\n")
cat("转录组数据维度:", dim(transcriptomics), "\n")
cat("蛋白质组数据维度:", dim(proteomics), "\n")
cat("样本信息维度:", dim(metadata), "\n")

# 检查缺失值
cat("基因组数据缺失值:", sum(is.na(genomics)), "\n")
cat("转录组数据缺失值:", sum(is.na(transcriptomics)), "\n")
cat("蛋白质组数据缺失值:", sum(is.na(proteomics)), "\n")
```

**参数解释：**
- `row.names = 1`：将第一列作为行名
- `dim()`：查看数据维度
- `is.na()`：检查缺失值

**预期输出：**
```
基因组数据维度: 100 50 
转录组数据维度: 1000 50 
蛋白质组数据维度: 500 50 
样本信息维度: 50 5 
基因组数据缺失值: 0 
转录组数据缺失值: 25 
蛋白质组数据缺失值: 10 
```

#### 2.2 数据标准化和预处理

**R代码：**
```r
# 转录组数据log2转换和标准化
transcriptomics_log <- log2(transcriptomics + 1)
transcriptomics_scaled <- scale(transcriptomics_log)

# 蛋白质组数据标准化
proteomics_scaled <- scale(proteomics)

# 基因组数据（已经是0/1编码，不需要标准化）
genomics_processed <- genomics

# 处理缺失值（使用均值填补）
transcriptomics_scaled[is.na(transcriptomics_scaled)] <- 0
proteomics_scaled[is.na(proteomics_scaled)] <- 0

# 保存预处理后的数据
write.csv(genomics_processed, "results/genomics_processed.csv")
write.csv(transcriptomics_scaled, "results/transcriptomics_scaled.csv")
write.csv(proteomics_scaled, "results/proteomics_scaled.csv")

cat("数据预处理完成！\n")
```

**结果验证：**
```r
# 验证标准化效果
cat("转录组数据标准化后均值:", round(mean(transcriptomics_scaled), 3), "\n")
cat("转录组数据标准化后标准差:", round(sd(transcriptomics_scaled), 3), "\n")
```

**检查点：** 确认数据已正确标准化，均值接近0，标准差接近1。

---

### 步骤3：多组学数据整合分析

#### 3.1 使用mixOmics进行PLS-DA分析

**操作说明：**
使用mixOmics包进行多组学数据的偏最小二乘判别分析，识别不同组间的差异特征。

**R代码：**
```r
# 加载mixOmics包
library(mixOmics)

# 准备数据列表
X <- list(genomics = t(genomics_processed),
          transcriptomics = t(transcriptomics_scaled),
          proteomics = t(proteomics_scaled))

# 准备分组信息（假设根据metadata中的group列）
Y <- metadata$group

# 执行多组学PLS-DA分析
result.diablo <- block.plsda(X, Y, ncomp = 2)

# 查看结果
print(result.diablo)

# 绘制样本得分图
pdf("plots/multiomics_plsda_samples.pdf", width = 10, height = 8)
plotIndiv(result.diablo, ind.names = FALSE, legend = TRUE, 
          title = "多组学PLS-DA样本得分图")
dev.off()

# 绘制变量载荷图
pdf("plots/multiomics_plsda_variables.pdf", width = 12, height = 8)
plotVar(result.diablo, cutoff = 0.5, title = "多组学PLS-DA变量载荷图")
dev.off()
```

**参数解释：**
- `ncomp = 2`：提取2个主成分
- `cutoff = 0.5`：载荷图中只显示载荷值大于0.5的变量
- `ind.names = FALSE`：不显示样本名称

#### 3.2 网络分析

**R代码：**
```r
# 构建相关性网络
pdf("plots/multiomics_network.pdf", width = 12, height = 10)
network(result.diablo, blocks = c(1,2,3), 
        color.node = c("darkorchid", "brown1", "lightgreen"),
        cutoff = 0.4)
dev.off()

# 提取重要特征
important.features <- selectVar(result.diablo, comp = 1)
print(important.features$genomics$name[1:10])
print(important.features$transcriptomics$name[1:10])
print(important.features$proteomics$name[1:10])
```

**预期输出：**
网络图显示不同组学数据间的相关性，重要特征列表显示对分组贡献最大的变量。

**检查点：** 确认生成了PLS-DA结果图和网络图，并识别出重要特征。

---

### 步骤4：机器学习模型构建

#### 4.1 特征选择和数据准备

**操作说明：**
使用Python进行机器学习分析，首先进行特征选择和数据准备。

**执行命令：**
```bash
# 启动Python环境
python
```

**Python代码：**
```python
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.metrics import classification_report, confusion_matrix
import matplotlib.pyplot as plt
import seaborn as sns

# 读取预处理后的数据
genomics = pd.read_csv("results/genomics_processed.csv", index_col=0)
transcriptomics = pd.read_csv("results/transcriptomics_scaled.csv", index_col=0)
proteomics = pd.read_csv("results/proteomics_scaled.csv", index_col=0)
metadata = pd.read_csv("data/metadata.csv", index_col=0)

# 合并所有组学数据
X_combined = pd.concat([genomics.T, transcriptomics.T, proteomics.T], axis=1)
y = metadata['group']

print(f"合并后数据维度: {X_combined.shape}")
print(f"分组信息: {y.value_counts()}")

# 特征选择（选择前500个最重要的特征）
selector = SelectKBest(f_classif, k=500)
X_selected = selector.fit_transform(X_combined, y)

print(f"特征选择后维度: {X_selected.shape}")
```

#### 4.2 模型训练和评估

**Python代码：**
```python
# 数据分割
X_train, X_test, y_train, y_test = train_test_split(
    X_selected, y, test_size=0.3, random_state=42, stratify=y)

print(f"训练集大小: {X_train.shape}")
print(f"测试集大小: {X_test.shape}")

# 随机森林模型
rf_model = RandomForestClassifier(n_estimators=100, random_state=42)
rf_model.fit(X_train, y_train)

# 预测和评估
y_pred_rf = rf_model.predict(X_test)
print("随机森林分类报告:")
print(classification_report(y_test, y_pred_rf))

# SVM模型
svm_model = SVC(kernel='rbf', random_state=42)
svm_model.fit(X_train, y_train)

y_pred_svm = svm_model.predict(X_test)
print("SVM分类报告:")
print(classification_report(y_test, y_pred_svm))

# 特征重要性分析
feature_importance = rf_model.feature_importances_
top_features_idx = np.argsort(feature_importance)[-20:]

plt.figure(figsize=(10, 8))
plt.barh(range(20), feature_importance[top_features_idx])
plt.title('Top 20 重要特征')
plt.xlabel('特征重要性')
plt.tight_layout()
plt.savefig('plots/feature_importance.pdf')
plt.show()
```

**参数解释：**
- `n_estimators=100`：随机森林中树的数量
- `test_size=0.3`：测试集占30%
- `stratify=y`：分层抽样，保持各组比例
- `kernel='rbf'`：SVM使用径向基函数核

**预期输出：**
```
合并后数据维度: (50, 1600)
分组信息: 
group1    25
group2    25
Name: group, dtype: int64
特征选择后维度: (50, 500)
训练集大小: (35, 500)
测试集大小: (15, 500)
```

#### 4.3 模型比较和可视化

**Python代码：**
```python
# 混淆矩阵可视化
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# 随机森林混淆矩阵
cm_rf = confusion_matrix(y_test, y_pred_rf)
sns.heatmap(cm_rf, annot=True, fmt='d', ax=axes[0], cmap='Blues')
axes[0].set_title('随机森林混淆矩阵')
axes[0].set_xlabel('预测标签')
axes[0].set_ylabel('真实标签')

# SVM混淆矩阵
cm_svm = confusion_matrix(y_test, y_pred_svm)
sns.heatmap(cm_svm, annot=True, fmt='d', ax=axes[1], cmap='Greens')
axes[1].set_title('SVM混淆矩阵')
axes[1].set_xlabel('预测标签')
axes[1].set_ylabel('真实标签')

plt.tight_layout()
plt.savefig('plots/confusion_matrices.pdf')
plt.show()

# 保存模型
import joblib
joblib.dump(rf_model, 'models/random_forest_model.pkl')
joblib.dump(svm_model, 'models/svm_model.pkl')

print("模型已保存到models/目录")
```

**检查点：** 确认模型训练完成，生成了分类报告和可视化图表。

---

### 步骤5：深度学习示例

#### 5.1 构建简单的神经网络

**Python代码：**
```python
import tensorflow as tf
from tensorflow import keras
from sklearn.preprocessing import LabelEncoder

# 标签编码
le = LabelEncoder()
y_encoded = le.fit_transform(y)
y_train_encoded = le.transform(y_train)
y_test_encoded = le.transform(y_test)

# 构建神经网络模型
model = keras.Sequential([
    keras.layers.Dense(256, activation='relu', input_shape=(X_selected.shape[1],)),
    keras.layers.Dropout(0.3),
    keras.layers.Dense(128, activation='relu'),
    keras.layers.Dropout(0.3),
    keras.layers.Dense(64, activation='relu'),
    keras.layers.Dense(len(np.unique(y)), activation='softmax')
])

# 编译模型
model.compile(optimizer='adam',
              loss='sparse_categorical_crossentropy',
              metrics=['accuracy'])

# 训练模型
history = model.fit(X_train, y_train_encoded,
                    epochs=50,
                    batch_size=8,
                    validation_split=0.2,
                    verbose=1)

# 评估模型
test_loss, test_acc = model.evaluate(X_test, y_test_encoded, verbose=0)
print(f"深度学习模型测试准确率: {test_acc:.4f}")

# 绘制训练历史
plt.figure(figsize=(12, 4))

plt.subplot(1, 2, 1)
plt.plot(history.history['loss'], label='训练损失')
plt.plot(history.history['val_loss'], label='验证损失')
plt.title('模型损失')
plt.xlabel('轮次')
plt.ylabel('损失')
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(history.history['accuracy'], label='训练准确率')
plt.plot(history.history['val_accuracy'], label='验证准确率')
plt.title('模型准确率')
plt.xlabel('轮次')
plt.ylabel('准确率')
plt.legend()

plt.tight_layout()
plt.savefig('plots/deep_learning_history.pdf')
plt.show()

# 保存模型
model.save('models/deep_learning_model.h5')
```

**参数解释：**
- `Dropout(0.3)`：30%的神经元随机失活，防止过拟合
- `epochs=50`：训练50个轮次
- `batch_size=8`：每批处理8个样本
- `validation_split=0.2`：20%的训练数据用于验证

**检查点：** 确认深度学习模型训练完成，生成了训练历史图表。

---

### 步骤6：云计算平台体验

#### 6.1 Google Colab使用

**操作说明：**
体验Google Colab进行云端机器学习分析。

**执行步骤：**
1. 打开浏览器，访问 https://colab.research.google.com/
2. 登录Google账号
3. 创建新的Notebook
4. 上传数据文件到Colab环境

**Colab代码示例：**
```python
# 在Colab中运行
from google.colab import files
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score

# 上传数据文件
uploaded = files.upload()

# 读取数据（假设上传了合并的数据文件）
# data = pd.read_csv('combined_data.csv')

# 使用GPU加速（如果可用）
import tensorflow as tf
print("GPU可用:", tf.config.list_physical_devices('GPU'))

# 简单的机器学习示例
# 这里可以运行之前的任何代码
print("Google Colab环境配置完成！")
```

#### 6.2 云端资源监控

**操作说明：**
在Colab中监控资源使用情况。

**Colab代码：**
```python
# 查看系统信息
!cat /proc/meminfo | grep MemTotal
!cat /proc/cpuinfo | grep "model name" | head -1
!nvidia-smi  # 如果有GPU

# 安装系统监控工具
!pip install psutil

import psutil
import matplotlib.pyplot as plt

# 监控内存使用
memory = psutil.virtual_memory()
print(f"总内存: {memory.total / (1024**3):.2f} GB")
print(f"可用内存: {memory.available / (1024**3):.2f} GB")
print(f"内存使用率: {memory.percent}%")

# 监控CPU使用
cpu_percent = psutil.cpu_percent(interval=1)
print(f"CPU使用率: {cpu_percent}%")
```

**检查点：** 成功在Google Colab中运行代码，了解云端资源使用情况。

---

### 步骤7：结果整合和报告生成

#### 7.1 结果汇总

**R代码：**
```r
# 生成综合分析报告
sink("results/analysis_summary.txt")

cat("=== 多组学数据整合分析报告 ===\n\n")

cat("1. 数据概况:\n")
cat("   - 基因组特征数:", ncol(genomics_processed), "\n")
cat("   - 转录组特征数:", ncol(transcriptomics_scaled), "\n")
cat("   - 蛋白质组特征数:", ncol(proteomics_scaled), "\n")
cat("   - 样本数:", nrow(metadata), "\n\n")

cat("2. PLS-DA分析结果:\n")
cat("   - 成功识别组间差异\n")
cat("   - 提取了重要的判别特征\n")
cat("   - 构建了多组学相关性网络\n\n")

cat("3. 机器学习结果:\n")
cat("   - 随机森林模型训练完成\n")
cat("   - SVM模型训练完成\n")
cat("   - 深度学习模型训练完成\n")
cat("   - 所有模型均保存在models/目录\n\n")

cat("4. 云计算体验:\n")
cat("   - 成功使用Google Colab\n")
cat("   - 体验了云端GPU加速\n")
cat("   - 了解了资源监控方法\n\n")

cat("分析完成时间:", Sys.time(), "\n")

sink()

cat("分析报告已生成: results/analysis_summary.txt\n")
```

#### 7.2 生成最终可视化

**Python代码：**
```python
# 创建综合结果图表
fig, axes = plt.subplots(2, 2, figsize=(15, 12))

# 1. 数据分布图
axes[0,0].hist([len(genomics.columns), len(transcriptomics.columns), len(proteomics.columns)], 
               bins=3, alpha=0.7, color=['blue', 'green', 'red'])
axes[0,0].set_title('各组学数据特征数量')
axes[0,0].set_xlabel('组学类型')
axes[0,0].set_ylabel('特征数量')
axes[0,0].set_xticks([0, 1, 2])
axes[0,0].set_xticklabels(['基因组', '转录组', '蛋白质组'])

# 2. 模型性能比较
models = ['随机森林', 'SVM', '深度学习']
accuracies = [0.85, 0.82, test_acc]  # 示例准确率
axes[0,1].bar(models, accuracies, color=['skyblue', 'lightgreen', 'orange'])
axes[0,1].set_title('模型性能比较')
axes[0,1].set_ylabel('准确率')
axes[0,1].set_ylim(0, 1)

# 3. 样本分组分布
group_counts = metadata['group'].value_counts()
axes[1,0].pie(group_counts.values, labels=group_counts.index, autopct='%1.1f%%')
axes[1,0].set_title('样本分组分布')

# 4. 特征重要性Top10
if 'feature_importance' in locals():
    top10_idx = np.argsort(feature_importance)[-10:]
    axes[1,1].barh(range(10), feature_importance[top10_idx])
    axes[1,1].set_title('Top 10 重要特征')
    axes[1,1].set_xlabel('重要性得分')

plt.tight_layout()
plt.savefig('plots/final_summary.pdf', dpi=300, bbox_inches='tight')
plt.show()

print("最终结果图表已生成: plots/final_summary.pdf")
```

**检查点：** 确认生成了完整的分析报告和综合结果图表。

## 预期结果

### 主要输出文件
1. **预处理数据**：`results/`目录下的标准化数据文件
   - 内容：经过质量控制和标准化的多组学数据
   - 用途：后续分析的输入数据

2. **分析结果**：`results/analysis_summary.txt`
   - 内容：完整的分析流程总结报告
   - 用途：记录分析过程和主要发现

3. **可视化图表**：`plots/`目录下的PDF文件
   - 内容：PLS-DA结果图、网络图、模型性能图等
   - 用途：结果展示和报告制作

4. **机器学习模型**：`models/`目录下的模型文件
   - 内容：训练好的随机森林、SVM和深度学习模型
   - 用途：后续预测和模型部署

### 关键结果指标
- **数据整合成功率**：应该成功整合所有组学数据
- **模型准确率**：各模型准确率应该在80%以上
- **特征识别**：应该识别出对分组有贡献的重要特征
- **可视化质量**：生成清晰、信息丰富的图表

### 成功标准
- [ ] 所有数据文件成功加载和预处理
- [ ] PLS-DA分析生成有意义的结果
- [ ] 机器学习模型训练无错误
- [ ] 深度学习模型收敛
- [ ] 云计算平台成功使用
- [ ] 生成完整的分析报告

## 故障排除

### 常见问题1：R包安装失败
**症状：** 安装mixOmics或其他R包时出错
**原因：** 网络问题或依赖包缺失
**解决方案：**
```r
# 更换CRAN镜像
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# 安装依赖包
install.packages(c("lattice", "ggplot2", "RColorBrewer"))

# 重新安装mixOmics
install.packages("mixOmics")
```

### 常见问题2：Python包导入错误
**症状：** ImportError: No module named 'sklearn'
**原因：** Python包未正确安装
**解决方案：**
```bash
# 检查Python环境
which python
pip list | grep scikit

# 重新安装
pip install --upgrade scikit-learn
pip install --upgrade tensorflow
```

### 常见问题3：内存不足错误
**症状：** 运行大数据分析时出现内存错误
**原因：** 数据量过大，超出可用内存
**解决方案：**
```python
# 减少特征数量
selector = SelectKBest(f_classif, k=200)  # 从500减少到200

# 使用数据分批处理
from sklearn.model_selection import StratifiedKFold
kfold = StratifiedKFold(n_splits=5)
```

### 常见问题4：Google Colab连接问题
**症状：** 无法访问Google Colab或上传失败
**原因：** 网络连接问题或浏览器设置
**解决方案：**
1. 检查网络连接
2. 清除浏览器缓存
3. 尝试使用其他浏览器
4. 使用VPN（如果需要）

### 获取帮助
如果遇到其他问题：
1. 检查错误日志：查看完整的错误信息
2. 查看软件文档：R和Python官方文档
3. 在线搜索：Stack Overflow等技术论坛
4. 联系助教或老师：wangys@hunau.edu.cn

## 扩展练习

### 练习1：自定义特征选择
**目标：** 实现自定义的特征选择算法
**任务：** 编写代码，结合统计检验和机器学习方法进行特征选择
**提示：** 可以结合t检验、方差分析和递归特征消除

### 练习2：集成学习模型
**目标：** 构建集成学习模型提高预测性能
**任务：** 使用投票、bagging或stacking方法集成多个模型
**提示：** 使用sklearn的VotingClassifier或StackingClassifier

### 练习3：深度学习架构优化
**目标：** 优化神经网络架构提高性能
**任务：** 尝试不同的网络结构、激活函数和正则化方法
**提示：** 可以尝试CNN、LSTM或Transformer架构

### 练习4：云端大规模分析
**目标：** 在云平台上处理更大规模的数据
**任务：** 使用云端GPU训练更复杂的深度学习模型
**提示：** 利用Google Colab Pro或其他云计算资源

### 思考问题
1. 不同数据整合策略的优缺点是什么？在什么情况下选择哪种策略？
2. 如何评估多组学整合分析的生物学意义？
3. 深度学习在小样本生物数据中的应用有哪些挑战？如何克服？
4. 云计算在生物信息学中的发展趋势是什么？

## 参考资料

### 相关文献
1. Rohart, F., et al. (2017). mixOmics: An R package for 'omics feature selection and multiple data integration. PLoS Computational Biology, 13(11), e1005752.
2. Argelaguet, R., et al. (2018). Multi-Omics Factor Analysis—a framework for unsupervised integration of multi-omics data sets. Molecular Systems Biology, 14(6), e8124.
3. Hasin, Y., et al. (2017). Multi-omics approaches to disease. Genome Biology, 18(1), 83.

### 在线资源
- mixOmics官方教程：http://mixomics.org/
- scikit-learn用户指南：https://scikit-learn.org/stable/user_guide.html
- TensorFlow教程：https://www.tensorflow.org/tutorials
- Google Colab使用指南：https://colab.research.google.com/

### 软件文档
- R官方文档：https://www.r-project.org/
- Python官方文档：https://docs.python.org/
- Jupyter文档：https://jupyter.org/documentation

## 附录

### 附录A：完整脚本文件
参见：`scripts/` 目录中的相关脚本文件
- `data_integration.R`：R语言数据整合脚本
- `ml_analysis.py`：Python机器学习分析脚本
- `deep_learning.py`：深度学习示例脚本
- `cloud_demo.sh`：云计算演示脚本

### 附录B：配置文件模板
```yaml
# conda环境配置文件 environment.yml
name: multiomics
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.8
  - pandas
  - numpy
  - scikit-learn
  - tensorflow
  - matplotlib
  - seaborn
  - jupyter
```

### 附录C：数据格式说明
- **基因组数据**：行为样本，列为SNP位点，值为0/1/2（基因型）
- **转录组数据**：行为样本，列为基因，值为表达量（FPKM/TPM）
- **蛋白质组数据**：行为样本，列为蛋白质，值为丰度
- **元数据**：行为样本，列为样本属性（年龄、性别、分组等）

---

**实验完成时间：** 预计 2 小时  
**难度等级：** 高级  
**最后更新：** 2025年