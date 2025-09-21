#!/usr/bin/env python3
"""
机器学习分析脚本
作者：王运生
日期：2025-01-21
用法：python ml_analysis.py
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split, cross_val_score, GridSearchCV
from sklearn.feature_selection import SelectKBest, f_classif, RFE
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report, confusion_matrix, roc_curve, auc
from sklearn.decomposition import PCA
import joblib
import os
import warnings
warnings.filterwarnings('ignore')

def setup_environment():
    """设置工作环境"""
    print("=== 机器学习分析 ===")
    
    # 检查工作目录
    work_dir = os.path.expanduser("~/ngs-analysis/lesson-08")
    if not os.path.exists(work_dir):
        print("错误：请先运行setup.sh设置环境")
        return False
    
    os.chdir(work_dir)
    
    # 创建输出目录
    os.makedirs("results", exist_ok=True)
    os.makedirs("plots", exist_ok=True)
    os.makedirs("models", exist_ok=True)
    
    return True

def load_data():
    """加载数据"""
    print("1. 加载数据...")
    
    try:
        # 加载预处理后的数据
        genomics = pd.read_csv("results/genomics_processed.csv", index_col=0)
        transcriptomics = pd.read_csv("results/transcriptomics_scaled.csv", index_col=0)
        proteomics = pd.read_csv("results/proteomics_scaled.csv", index_col=0)
        metadata = pd.read_csv("data/metadata.csv", index_col=0)
        
        print(f"  基因组数据: {genomics.shape}")
        print(f"  转录组数据: {transcriptomics.shape}")
        print(f"  蛋白质组数据: {proteomics.shape}")
        print(f"  样本信息: {metadata.shape}")
        
        return genomics, transcriptomics, proteomics, metadata
    
    except FileNotFoundError as e:
        print(f"错误：数据文件未找到 - {e}")
        print("请先运行data_integration.R进行数据预处理")
        return None, None, None, None

def prepare_features(genomics, transcriptomics, proteomics, metadata):
    """准备特征和标签"""
    print("\n2. 准备特征和标签...")
    
    # 合并所有组学数据
    X_combined = pd.concat([genomics.T, transcriptomics.T, proteomics.T], axis=1)
    y = metadata['group']
    
    print(f"  合并后数据维度: {X_combined.shape}")
    print(f"  分组分布: {y.value_counts().to_dict()}")
    
    # 检查数据质量
    print(f"  缺失值数量: {X_combined.isnull().sum().sum()}")
    print(f"  无穷值数量: {np.isinf(X_combined.values).sum()}")
    
    # 处理异常值
    X_combined = X_combined.fillna(0)
    X_combined = X_combined.replace([np.inf, -np.inf], 0)
    
    return X_combined, y

def feature_selection(X, y, n_features=500):
    """特征选择"""
    print(f"\n3. 特征选择（选择前{n_features}个特征）...")
    
    # 方法1：单变量特征选择
    selector_univariate = SelectKBest(f_classif, k=n_features)
    X_selected_univariate = selector_univariate.fit_transform(X, y)
    
    # 方法2：递归特征消除（使用随机森林）
    rf_temp = RandomForestClassifier(n_estimators=50, random_state=42)
    selector_rfe = RFE(rf_temp, n_features_to_select=n_features)
    X_selected_rfe = selector_rfe.fit_transform(X, y)
    
    print(f"  单变量选择后维度: {X_selected_univariate.shape}")
    print(f"  RFE选择后维度: {X_selected_rfe.shape}")
    
    # 使用单变量选择的结果
    selected_features = selector_univariate.get_support(indices=True)
    feature_names = X.columns[selected_features]
    
    return X_selected_univariate, feature_names, selector_univariate

def split_data(X, y, test_size=0.3):
    """数据分割"""
    print(f"\n4. 数据分割（测试集比例: {test_size}）...")
    
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_size, random_state=42, stratify=y
    )
    
    print(f"  训练集大小: {X_train.shape}")
    print(f"  测试集大小: {X_test.shape}")
    print(f"  训练集分组: {pd.Series(y_train).value_counts().to_dict()}")
    print(f"  测试集分组: {pd.Series(y_test).value_counts().to_dict()}")
    
    return X_train, X_test, y_train, y_test

def train_models(X_train, X_test, y_train, y_test):
    """训练多个机器学习模型"""
    print("\n5. 训练机器学习模型...")
    
    models = {}
    results = {}
    
    # 随机森林
    print("  训练随机森林...")
    rf_model = RandomForestClassifier(n_estimators=100, random_state=42)
    rf_model.fit(X_train, y_train)
    y_pred_rf = rf_model.predict(X_test)
    models['RandomForest'] = rf_model
    results['RandomForest'] = {
        'predictions': y_pred_rf,
        'probabilities': rf_model.predict_proba(X_test)
    }
    
    # SVM
    print("  训练SVM...")
    svm_model = SVC(kernel='rbf', probability=True, random_state=42)
    svm_model.fit(X_train, y_train)
    y_pred_svm = svm_model.predict(X_test)
    models['SVM'] = svm_model
    results['SVM'] = {
        'predictions': y_pred_svm,
        'probabilities': svm_model.predict_proba(X_test)
    }
    
    # 逻辑回归
    print("  训练逻辑回归...")
    lr_model = LogisticRegression(random_state=42, max_iter=1000)
    lr_model.fit(X_train, y_train)
    y_pred_lr = lr_model.predict(X_test)
    models['LogisticRegression'] = lr_model
    results['LogisticRegression'] = {
        'predictions': y_pred_lr,
        'probabilities': lr_model.predict_proba(X_test)
    }
    
    return models, results

def evaluate_models(models, results, y_test):
    """评估模型性能"""
    print("\n6. 模型性能评估...")
    
    evaluation_results = {}
    
    for model_name, model in models.items():
        print(f"\n  {model_name} 分类报告:")
        y_pred = results[model_name]['predictions']
        
        # 分类报告
        report = classification_report(y_test, y_pred, output_dict=True)
        print(classification_report(y_test, y_pred))
        
        # 准确率
        accuracy = report['accuracy']
        evaluation_results[model_name] = {
            'accuracy': accuracy,
            'classification_report': report,
            'predictions': y_pred
        }
    
    return evaluation_results

def visualize_results(models, results, evaluation_results, y_test, feature_names):
    """可视化结果"""
    print("\n7. 生成可视化图表...")
    
    # 设置中文字体
    plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    
    # 1. 模型性能比较
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # 准确率比较
    model_names = list(evaluation_results.keys())
    accuracies = [evaluation_results[name]['accuracy'] for name in model_names]
    
    axes[0,0].bar(model_names, accuracies, color=['skyblue', 'lightgreen', 'orange'])
    axes[0,0].set_title('模型准确率比较')
    axes[0,0].set_ylabel('准确率')
    axes[0,0].set_ylim(0, 1)
    for i, acc in enumerate(accuracies):
        axes[0,0].text(i, acc + 0.01, f'{acc:.3f}', ha='center')
    
    # 混淆矩阵（随机森林）
    cm_rf = confusion_matrix(y_test, results['RandomForest']['predictions'])
    sns.heatmap(cm_rf, annot=True, fmt='d', ax=axes[0,1], cmap='Blues')
    axes[0,1].set_title('随机森林混淆矩阵')
    axes[0,1].set_xlabel('预测标签')
    axes[0,1].set_ylabel('真实标签')
    
    # 特征重要性（随机森林）
    rf_model = models['RandomForest']
    feature_importance = rf_model.feature_importances_
    top_features_idx = np.argsort(feature_importance)[-20:]
    
    axes[1,0].barh(range(20), feature_importance[top_features_idx])
    axes[1,0].set_title('Top 20 重要特征（随机森林）')
    axes[1,0].set_xlabel('特征重要性')
    
    # ROC曲线
    for model_name in model_names:
        y_proba = results[model_name]['probabilities'][:, 1]  # 假设是二分类
        fpr, tpr, _ = roc_curve(y_test == y_test.unique()[1], y_proba)
        roc_auc = auc(fpr, tpr)
        axes[1,1].plot(fpr, tpr, label=f'{model_name} (AUC = {roc_auc:.3f})')
    
    axes[1,1].plot([0, 1], [0, 1], 'k--')
    axes[1,1].set_xlim([0.0, 1.0])
    axes[1,1].set_ylim([0.0, 1.05])
    axes[1,1].set_xlabel('假正率')
    axes[1,1].set_ylabel('真正率')
    axes[1,1].set_title('ROC曲线')
    axes[1,1].legend()
    
    plt.tight_layout()
    plt.savefig('plots/ml_model_comparison.pdf', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 2. 详细的混淆矩阵
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    for i, model_name in enumerate(model_names):
        cm = confusion_matrix(y_test, results[model_name]['predictions'])
        sns.heatmap(cm, annot=True, fmt='d', ax=axes[i], cmap='Blues')
        axes[i].set_title(f'{model_name} 混淆矩阵')
        axes[i].set_xlabel('预测标签')
        axes[i].set_ylabel('真实标签')
    
    plt.tight_layout()
    plt.savefig('plots/confusion_matrices.pdf', dpi=300, bbox_inches='tight')
    plt.show()

def hyperparameter_tuning(X_train, y_train):
    """超参数调优"""
    print("\n8. 超参数调优...")
    
    # 随机森林超参数调优
    print("  随机森林超参数调优...")
    rf_param_grid = {
        'n_estimators': [50, 100, 200],
        'max_depth': [None, 10, 20],
        'min_samples_split': [2, 5, 10]
    }
    
    rf_grid = GridSearchCV(
        RandomForestClassifier(random_state=42),
        rf_param_grid,
        cv=5,
        scoring='accuracy',
        n_jobs=-1
    )
    rf_grid.fit(X_train, y_train)
    
    print(f"    最佳参数: {rf_grid.best_params_}")
    print(f"    最佳得分: {rf_grid.best_score_:.3f}")
    
    # SVM超参数调优
    print("  SVM超参数调优...")
    svm_param_grid = {
        'C': [0.1, 1, 10],
        'gamma': ['scale', 'auto', 0.1, 1]
    }
    
    svm_grid = GridSearchCV(
        SVC(random_state=42),
        svm_param_grid,
        cv=5,
        scoring='accuracy',
        n_jobs=-1
    )
    svm_grid.fit(X_train, y_train)
    
    print(f"    最佳参数: {svm_grid.best_params_}")
    print(f"    最佳得分: {svm_grid.best_score_:.3f}")
    
    return rf_grid.best_estimator_, svm_grid.best_estimator_

def cross_validation(models, X_train, y_train):
    """交叉验证"""
    print("\n9. 交叉验证...")
    
    cv_results = {}
    
    for model_name, model in models.items():
        scores = cross_val_score(model, X_train, y_train, cv=5, scoring='accuracy')
        cv_results[model_name] = {
            'mean': scores.mean(),
            'std': scores.std(),
            'scores': scores
        }
        print(f"  {model_name}: {scores.mean():.3f} (+/- {scores.std() * 2:.3f})")
    
    return cv_results

def save_models(models):
    """保存模型"""
    print("\n10. 保存模型...")
    
    for model_name, model in models.items():
        filename = f"models/{model_name.lower()}_model.pkl"
        joblib.dump(model, filename)
        print(f"  {model_name} 已保存到 {filename}")

def generate_report(evaluation_results, cv_results, feature_names):
    """生成分析报告"""
    print("\n11. 生成分析报告...")
    
    with open("results/ml_analysis_report.txt", "w", encoding='utf-8') as f:
        f.write("=== 机器学习分析报告 ===\n")
        f.write(f"分析时间: {pd.Timestamp.now()}\n\n")
        
        f.write("1. 数据概况:\n")
        f.write(f"   特征数量: {len(feature_names)}\n")
        f.write(f"   样本数量: 训练集 + 测试集\n\n")
        
        f.write("2. 模型性能:\n")
        for model_name, results in evaluation_results.items():
            f.write(f"   {model_name}:\n")
            f.write(f"     准确率: {results['accuracy']:.3f}\n")
            f.write(f"     交叉验证: {cv_results[model_name]['mean']:.3f} (+/- {cv_results[model_name]['std'] * 2:.3f})\n")
        
        f.write("\n3. 特征选择:\n")
        f.write(f"   选择的特征数量: {len(feature_names)}\n")
        f.write("   特征选择方法: 单变量统计检验\n\n")
        
        f.write("4. 输出文件:\n")
        f.write("   模型文件: models/\n")
        f.write("   可视化图表: plots/\n")
        f.write("   分析报告: results/ml_analysis_report.txt\n\n")
        
        f.write("分析完成！\n")
    
    print("  分析报告已保存到 results/ml_analysis_report.txt")

def main():
    """主函数"""
    # 设置环境
    if not setup_environment():
        return
    
    # 加载数据
    genomics, transcriptomics, proteomics, metadata = load_data()
    if genomics is None:
        return
    
    # 准备特征
    X_combined, y = prepare_features(genomics, transcriptomics, proteomics, metadata)
    
    # 特征选择
    X_selected, feature_names, selector = feature_selection(X_combined, y)
    
    # 数据分割
    X_train, X_test, y_train, y_test = split_data(X_selected, y)
    
    # 训练模型
    models, results = train_models(X_train, X_test, y_train, y_test)
    
    # 评估模型
    evaluation_results = evaluate_models(models, results, y_test)
    
    # 交叉验证
    cv_results = cross_validation(models, X_train, y_train)
    
    # 可视化结果
    visualize_results(models, results, evaluation_results, y_test, feature_names)
    
    # 超参数调优（可选）
    # best_rf, best_svm = hyperparameter_tuning(X_train, y_train)
    
    # 保存模型
    save_models(models)
    
    # 生成报告
    generate_report(evaluation_results, cv_results, feature_names)
    
    print("\n=== 机器学习分析完成 ===")
    print("输出文件:")
    print("  - 模型文件: models/")
    print("  - 可视化图表: plots/")
    print("  - 分析报告: results/ml_analysis_report.txt")

if __name__ == "__main__":
    main()