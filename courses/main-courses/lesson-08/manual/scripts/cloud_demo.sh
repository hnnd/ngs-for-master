#!/bin/bash
# 云计算平台演示脚本
# 作者：王运生
# 日期：2025-01-21
# 用法：bash cloud_demo.sh

set -e  # 遇到错误立即退出

echo "=== 云计算平台演示 ==="

# 检查工作目录
if [ ! -d "~/ngs-analysis/lesson-08" ]; then
    echo "错误：请先运行setup.sh设置环境"
    exit 1
fi

cd ~/ngs-analysis/lesson-08

# 创建云计算演示目录
mkdir -p cloud_demo
cd cloud_demo

echo "1. 创建Google Colab演示文件..."

# 创建简化的演示脚本
cat > colab_demo.py << 'EOF'
"""
Google Colab 多组学数据分析演示
适用于第8次课云计算体验
"""

# 检查环境
print("=== Google Colab 环境检查 ===")
import sys
print(f"Python版本: {sys.version}")

# 检查GPU
try:
    import tensorflow as tf
    print(f"TensorFlow版本: {tf.__version__}")
    print(f"GPU可用: {tf.config.list_physical_devices('GPU')}")
except ImportError:
    print("TensorFlow未安装")

# 安装必要包
print("\n安装必要的包...")
import subprocess
import sys

packages = ['pandas', 'numpy', 'scikit-learn', 'matplotlib', 'seaborn']
for package in packages:
    try:
        __import__(package)
        print(f"✓ {package} 已安装")
    except ImportError:
        print(f"安装 {package}...")
        subprocess.check_call([sys.executable, '-m', 'pip', 'install', package])

# 生成示例数据
print("\n=== 生成示例多组学数据 ===")
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report
import matplotlib.pyplot as plt

# 设置随机种子
np.random.seed(42)

# 生成多组学数据
n_samples = 100
genomics_data = np.random.choice([0, 1, 2], size=(n_samples, 200))
transcriptomics_data = np.random.lognormal(5, 2, size=(n_samples, 500))
proteomics_data = np.random.normal(10, 3, size=(n_samples, 300))

# 合并数据
X = np.hstack([genomics_data, transcriptomics_data, proteomics_data])
y = np.array(['group1'] * 50 + ['group2'] * 50)

print(f"数据维度: {X.shape}")
print(f"分组分布: {np.unique(y, return_counts=True)}")

# 机器学习分析
print("\n=== 机器学习分析 ===")
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

# 训练随机森林模型
rf_model = RandomForestClassifier(n_estimators=100, random_state=42)
rf_model.fit(X_train, y_train)
y_pred = rf_model.predict(X_test)

accuracy = accuracy_score(y_test, y_pred)
print(f"随机森林准确率: {accuracy:.3f}")
print("\n分类报告:")
print(classification_report(y_test, y_pred))

# 可视化结果
print("\n=== 结果可视化 ===")
from sklearn.decomposition import PCA

# PCA降维可视化
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X)

plt.figure(figsize=(10, 6))
for i, group in enumerate(np.unique(y)):
    mask = y == group
    plt.scatter(X_pca[mask, 0], X_pca[mask, 1], label=group, alpha=0.7)

plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%})')
plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%})')
plt.title('多组学数据PCA可视化')
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()

print("\n=== 演示完成 ===")
print("这个演示展示了如何在云端进行多组学数据分析")
EOF

echo "✓ Google Colab演示脚本已创建"

echo ""
echo "2. 创建云平台使用指南..."

cat > cloud_usage_guide.md << 'EOF'
# 云计算平台使用指南

## Google Colab 使用步骤

### 1. 访问和设置
1. 打开浏览器，访问 https://colab.research.google.com/
2. 使用Google账号登录
3. 创建新的Notebook或上传现有文件

### 2. 环境配置
```python
# 检查Python环境
import sys
print(f"Python版本: {sys.version}")

# 检查GPU可用性
import tensorflow as tf
print("GPU可用:", tf.config.list_physical_devices('GPU'))

# 安装额外的包
!pip install pandas numpy scikit-learn matplotlib seaborn
```

### 3. 数据上传
```python
# 方法1: 直接上传文件
from google.colab import files
uploaded = files.upload()

# 方法2: 从Google Drive加载
from google.colab import drive
drive.mount('/content/drive')
```

### 4. 运行分析
- 运行单个代码块: Shift + Enter
- 运行所有代码: Runtime → Run all
- 重启运行时: Runtime → Restart runtime

### 5. 保存结果
```python
# 下载文件
from google.colab import files
files.download('result.csv')

# 保存到Google Drive
import shutil
shutil.copy('result.csv', '/content/drive/My Drive/')
```

## 其他云平台简介

### AWS (Amazon Web Services)
- **优势**: 服务最全面，生态系统成熟
- **适用**: 企业级应用，大规模部署
- **成本**: 相对较高，但功能强大
- **入门**: 12个月免费套餐

### 阿里云
- **优势**: 中国本土化，价格优势
- **适用**: 中国用户，成本敏感项目
- **成本**: 相对较低
- **入门**: 新用户优惠活动

### Azure (微软云)
- **优势**: 与Microsoft生态集成好
- **适用**: 企业用户，混合云
- **成本**: 中等水平
- **入门**: 12个月免费套餐

## 云计算优势

### 1. 成本效益
- 按需付费，无需前期投资
- 避免硬件维护成本
- 弹性扩缩容

### 2. 便利性
- 随时随地访问
- 无需本地配置
- 自动备份和恢复

### 3. 协作性
- 易于分享和协作
- 版本控制
- 团队协作

### 4. 性能
- 高性能计算资源
- GPU加速
- 大规模并行处理

## 注意事项

### 1. 数据安全
- 注意数据隐私保护
- 了解数据存储位置
- 使用加密传输

### 2. 成本控制
- 监控资源使用
- 设置预算告警
- 及时释放不用的资源

### 3. 网络依赖
- 需要稳定的网络连接
- 大文件传输耗时
- 考虑网络成本

## 实践建议

### 对于学习者
- 从免费服务开始 (Google Colab)
- 熟悉基本操作
- 逐步尝试付费服务

### 对于研究者
- 根据数据位置选择平台
- 考虑合规要求
- 评估长期成本

### 对于企业
- 多云策略
- 专业技术支持
- 安全和合规优先
EOF

echo "✓ 云平台使用指南已创建"

echo ""
echo "3. 创建资源监控脚本..."

cat > monitor_system.py << 'EOF'
#!/usr/bin/env python3
"""
简化的系统资源监控脚本
用于云计算环境演示
"""

import time
import psutil
from datetime import datetime

def check_system_resources():
    """检查系统资源使用情况"""
    print("=== 系统资源监控 ===")
    print(f"检查时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("-" * 40)
    
    # CPU信息
    cpu_count = psutil.cpu_count()
    cpu_percent = psutil.cpu_percent(interval=1)
    print(f"CPU核心数: {cpu_count}")
    print(f"CPU使用率: {cpu_percent:.1f}%")
    
    # 内存信息
    memory = psutil.virtual_memory()
    print(f"总内存: {memory.total / (1024**3):.2f} GB")
    print(f"可用内存: {memory.available / (1024**3):.2f} GB")
    print(f"内存使用率: {memory.percent:.1f}%")
    
    # 磁盘信息
    disk = psutil.disk_usage('/')
    print(f"磁盘总空间: {disk.total / (1024**3):.2f} GB")
    print(f"磁盘可用空间: {disk.free / (1024**3):.2f} GB")
    print(f"磁盘使用率: {disk.percent:.1f}%")
    
    # 网络信息
    net_io = psutil.net_io_counters()
    print(f"网络发送: {net_io.bytes_sent / (1024**2):.1f} MB")
    print(f"网络接收: {net_io.bytes_recv / (1024**2):.1f} MB")
    
    print("-" * 40)
    
    # 性能建议
    if cpu_percent > 80:
        print("⚠️  CPU使用率较高")
    if memory.percent > 80:
        print("⚠️  内存使用率较高")
    if disk.percent > 90:
        print("⚠️  磁盘空间不足")
    
    if cpu_percent < 50 and memory.percent < 50:
        print("✅ 系统资源充足")

def monitor_for_duration(duration=30):
    """持续监控指定时间"""
    print(f"开始监控 {duration} 秒...")
    start_time = time.time()
    
    while time.time() - start_time < duration:
        cpu = psutil.cpu_percent(interval=1)
        memory = psutil.virtual_memory().percent
        
        print(f"\r时间: {int(time.time() - start_time):2d}s | "
              f"CPU: {cpu:5.1f}% | 内存: {memory:5.1f}%", 
              end='', flush=True)
        
        time.sleep(2)
    
    print("\n监控完成!")

if __name__ == "__main__":
    try:
        check_system_resources()
        print("\n是否进行持续监控? (y/n): ", end='')
        choice = input().lower()
        
        if choice == 'y':
            monitor_for_duration(30)
    
    except KeyboardInterrupt:
        print("\n监控被中断")
    except Exception as e:
        print(f"监控出错: {e}")
EOF

chmod +x monitor_system.py

echo "✓ 系统监控脚本已创建"

echo ""
echo "4. 创建云计算对比表..."

cat > cloud_comparison.txt << 'EOF'
云计算平台对比表

平台特性对比:
=====================================
特性          Google Cloud  AWS    Azure   阿里云   华为云
-------------------------------------
免费额度      Colab免费     300$   200$    优惠     优惠
GPU支持       免费GPU       付费   付费    付费     付费
中文支持      有限          有限   有限    完整     完整
国内访问      需要VPN       需要   部分    直接     直接
学习资源      丰富          丰富   中等    中等     一般

生物信息学适用性:
=====================================
用途                推荐平台
-------------------------------------
学习和教学          Google Colab
小规模研究          阿里云/华为云
大规模分析          AWS/Azure
企业级应用          AWS/Azure/阿里云
成本敏感项目        阿里云/华为云
国际合作项目        AWS/Azure
政府/国企项目       华为云/阿里云

价格对比 (大概):
=====================================
资源类型        Google Cloud  AWS    阿里云
-------------------------------------
2核4GB/小时     $0.10        $0.10  $0.05
GPU/小时        $2.50        $3.00  $1.80
存储/GB/月      $0.02        $0.02  $0.015

选择建议:
=====================================
1. 学习阶段: Google Colab (免费GPU)
2. 研究项目: 阿里云 (成本低，中文支持)
3. 企业应用: AWS (功能全面)
4. 国际项目: AWS/Azure (全球覆盖)
5. 国产化需求: 华为云 (自主可控)
EOF

echo "✓ 云计算对比表已创建"

echo ""
echo "=== 云计算演示环境设置完成 ==="
echo ""
echo "创建的文件:"
echo "  1. colab_demo.py - Google Colab演示脚本"
echo "  2. cloud_usage_guide.md - 云平台使用指南"
echo "  3. monitor_system.py - 系统监控脚本"
echo "  4. cloud_comparison.txt - 云平台对比表"
echo ""
echo "使用方法:"
echo "  • 将 colab_demo.py 内容复制到 Google Colab"
echo "  • 运行 python monitor_system.py 监控系统"
echo "  • 查看指南了解各云平台特点"
echo ""
echo "云计算平台体验准备完成！"