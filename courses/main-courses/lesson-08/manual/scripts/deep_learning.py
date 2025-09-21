#!/usr/bin/env python3
"""
深度学习分析脚本
作者：王运生
日期：2025-01-21
用法：python deep_learning.py
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import classification_report, confusion_matrix
import os
import warnings
warnings.filterwarnings('ignore')

# 设置随机种子
np.random.seed(42)
tf.random.set_seed(42)

def setup_environment():
    """设置工作环境"""
    print("=== 深度学习分析 ===")
    
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
    
    # 检查GPU
    print("GPU可用性检查:")
    gpus = tf.config.list_physical_devices('GPU')
    if gpus:
        print(f"  发现 {len(gpus)} 个GPU设备")
        for gpu in gpus:
            print(f"    {gpu}")
    else:
        print("  未发现GPU设备，将使用CPU")
    
    return True

def load_and_prepare_data():
    """加载和准备数据"""
    print("\n1. 加载和准备数据...")
    
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
        
        # 合并所有组学数据
        X_combined = pd.concat([genomics.T, transcriptomics.T, proteomics.T], axis=1)
        y = metadata['group']
        
        print(f"  合并后数据维度: {X_combined.shape}")
        print(f"  分组分布: {y.value_counts().to_dict()}")
        
        # 处理异常值
        X_combined = X_combined.fillna(0)
        X_combined = X_combined.replace([np.inf, -np.inf], 0)
        
        # 标签编码
        le = LabelEncoder()
        y_encoded = le.fit_transform(y)
        
        print(f"  标签编码: {dict(zip(le.classes_, range(len(le.classes_))))}")
        
        return X_combined.values, y_encoded, le, X_combined.columns
    
    except FileNotFoundError as e:
        print(f"错误：数据文件未找到 - {e}")
        print("请先运行data_integration.R进行数据预处理")
        return None, None, None, None

def create_autoencoder(input_dim, encoding_dim=64):
    """创建自编码器模型"""
    print(f"\n2. 创建自编码器模型（编码维度: {encoding_dim}）...")
    
    # 编码器
    input_layer = keras.Input(shape=(input_dim,))
    encoded = layers.Dense(512, activation='relu')(input_layer)
    encoded = layers.Dropout(0.3)(encoded)
    encoded = layers.Dense(256, activation='relu')(encoded)
    encoded = layers.Dropout(0.3)(encoded)
    encoded = layers.Dense(encoding_dim, activation='relu')(encoded)
    
    # 解码器
    decoded = layers.Dense(256, activation='relu')(encoded)
    decoded = layers.Dropout(0.3)(decoded)
    decoded = layers.Dense(512, activation='relu')(decoded)
    decoded = layers.Dropout(0.3)(decoded)
    decoded = layers.Dense(input_dim, activation='linear')(decoded)
    
    # 自编码器模型
    autoencoder = keras.Model(input_layer, decoded)
    encoder = keras.Model(input_layer, encoded)
    
    autoencoder.compile(optimizer='adam', loss='mse', metrics=['mae'])
    
    print("  自编码器架构:")
    autoencoder.summary()
    
    return autoencoder, encoder

def create_classifier(input_dim, num_classes):
    """创建分类器模型"""
    print(f"\n3. 创建分类器模型（类别数: {num_classes}）...")
    
    model = keras.Sequential([
        layers.Dense(512, activation='relu', input_shape=(input_dim,)),
        layers.BatchNormalization(),
        layers.Dropout(0.4),
        
        layers.Dense(256, activation='relu'),
        layers.BatchNormalization(),
        layers.Dropout(0.4),
        
        layers.Dense(128, activation='relu'),
        layers.BatchNormalization(),
        layers.Dropout(0.3),
        
        layers.Dense(64, activation='relu'),
        layers.Dropout(0.3),
        
        layers.Dense(num_classes, activation='softmax')
    ])
    
    model.compile(
        optimizer='adam',
        loss='sparse_categorical_crossentropy',
        metrics=['accuracy']
    )
    
    print("  分类器架构:")
    model.summary()
    
    return model

def train_autoencoder(autoencoder, X_train, X_val, epochs=100):
    """训练自编码器"""
    print(f"\n4. 训练自编码器（{epochs}轮）...")
    
    # 回调函数
    callbacks = [
        keras.callbacks.EarlyStopping(patience=10, restore_best_weights=True),
        keras.callbacks.ReduceLROnPlateau(factor=0.5, patience=5)
    ]
    
    # 训练
    history = autoencoder.fit(
        X_train, X_train,
        epochs=epochs,
        batch_size=16,
        validation_data=(X_val, X_val),
        callbacks=callbacks,
        verbose=1
    )
    
    return history

def train_classifier(model, X_train, X_val, y_train, y_val, epochs=100):
    """训练分类器"""
    print(f"\n5. 训练分类器（{epochs}轮）...")
    
    # 回调函数
    callbacks = [
        keras.callbacks.EarlyStopping(patience=15, restore_best_weights=True),
        keras.callbacks.ReduceLROnPlateau(factor=0.5, patience=7)
    ]
    
    # 训练
    history = model.fit(
        X_train, y_train,
        epochs=epochs,
        batch_size=16,
        validation_data=(X_val, y_val),
        callbacks=callbacks,
        verbose=1
    )
    
    return history

def visualize_training_history(history, title, save_path):
    """可视化训练历史"""
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    
    # 损失
    axes[0].plot(history.history['loss'], label='训练损失')
    axes[0].plot(history.history['val_loss'], label='验证损失')
    axes[0].set_title(f'{title} - 损失')
    axes[0].set_xlabel('轮次')
    axes[0].set_ylabel('损失')
    axes[0].legend()
    axes[0].grid(True)
    
    # 准确率（如果有）
    if 'accuracy' in history.history:
        axes[1].plot(history.history['accuracy'], label='训练准确率')
        axes[1].plot(history.history['val_accuracy'], label='验证准确率')
        axes[1].set_title(f'{title} - 准确率')
        axes[1].set_xlabel('轮次')
        axes[1].set_ylabel('准确率')
        axes[1].legend()
        axes[1].grid(True)
    else:
        # MAE（自编码器）
        axes[1].plot(history.history['mae'], label='训练MAE')
        axes[1].plot(history.history['val_mae'], label='验证MAE')
        axes[1].set_title(f'{title} - MAE')
        axes[1].set_xlabel('轮次')
        axes[1].set_ylabel('MAE')
        axes[1].legend()
        axes[1].grid(True)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()

def evaluate_classifier(model, X_test, y_test, le):
    """评估分类器"""
    print("\n6. 评估分类器...")
    
    # 预测
    y_pred_proba = model.predict(X_test)
    y_pred = np.argmax(y_pred_proba, axis=1)
    
    # 计算准确率
    test_loss, test_acc = model.evaluate(X_test, y_test, verbose=0)
    print(f"  测试准确率: {test_acc:.4f}")
    print(f"  测试损失: {test_loss:.4f}")
    
    # 分类报告
    print("\n  分类报告:")
    print(classification_report(y_test, y_pred, target_names=le.classes_))
    
    # 混淆矩阵
    cm = confusion_matrix(y_test, y_pred)
    
    plt.figure(figsize=(8, 6))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', 
                xticklabels=le.classes_, yticklabels=le.classes_)
    plt.title('深度学习模型混淆矩阵')
    plt.xlabel('预测标签')
    plt.ylabel('真实标签')
    plt.tight_layout()
    plt.savefig('plots/dl_confusion_matrix.pdf', dpi=300, bbox_inches='tight')
    plt.show()
    
    return test_acc, y_pred, y_pred_proba

def visualize_encoded_features(encoder, X_test, y_test, le):
    """可视化编码特征"""
    print("\n7. 可视化编码特征...")
    
    # 获取编码特征
    encoded_features = encoder.predict(X_test)
    
    # 如果编码维度大于2，使用PCA降维到2D
    if encoded_features.shape[1] > 2:
        from sklearn.decomposition import PCA
        pca = PCA(n_components=2)
        encoded_2d = pca.fit_transform(encoded_features)
        print(f"  PCA解释方差比: {pca.explained_variance_ratio_}")
    else:
        encoded_2d = encoded_features
    
    # 绘制2D散点图
    plt.figure(figsize=(10, 8))
    colors = ['red', 'blue', 'green', 'orange', 'purple']
    
    for i, class_name in enumerate(le.classes_):
        mask = y_test == i
        plt.scatter(encoded_2d[mask, 0], encoded_2d[mask, 1], 
                   c=colors[i % len(colors)], label=class_name, alpha=0.7)
    
    plt.xlabel('编码特征 1')
    plt.ylabel('编码特征 2')
    plt.title('自编码器编码特征可视化')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('plots/encoded_features_visualization.pdf', dpi=300, bbox_inches='tight')
    plt.show()

def create_variational_autoencoder(input_dim, latent_dim=32):
    """创建变分自编码器（VAE）"""
    print(f"\n8. 创建变分自编码器（潜在维度: {latent_dim}）...")
    
    # 编码器
    encoder_inputs = keras.Input(shape=(input_dim,))
    x = layers.Dense(512, activation='relu')(encoder_inputs)
    x = layers.Dense(256, activation='relu')(x)
    
    z_mean = layers.Dense(latent_dim)(x)
    z_log_var = layers.Dense(latent_dim)(x)
    
    # 采样函数
    def sampling(args):
        z_mean, z_log_var = args
        batch = tf.shape(z_mean)[0]
        dim = tf.shape(z_mean)[1]
        epsilon = tf.keras.backend.random_normal(shape=(batch, dim))
        return z_mean + tf.exp(0.5 * z_log_var) * epsilon
    
    z = layers.Lambda(sampling)([z_mean, z_log_var])
    
    # 解码器
    decoder_inputs = keras.Input(shape=(latent_dim,))
    x = layers.Dense(256, activation='relu')(decoder_inputs)
    x = layers.Dense(512, activation='relu')(x)
    decoder_outputs = layers.Dense(input_dim, activation='linear')(x)
    
    # 模型
    encoder = keras.Model(encoder_inputs, [z_mean, z_log_var, z])
    decoder = keras.Model(decoder_inputs, decoder_outputs)
    
    # VAE
    vae_outputs = decoder(z)
    vae = keras.Model(encoder_inputs, vae_outputs)
    
    # 损失函数
    reconstruction_loss = keras.losses.mse(encoder_inputs, vae_outputs)
    reconstruction_loss *= input_dim
    kl_loss = 1 + z_log_var - tf.square(z_mean) - tf.exp(z_log_var)
    kl_loss = tf.reduce_mean(kl_loss)
    kl_loss *= -0.5
    vae_loss = tf.reduce_mean(reconstruction_loss + kl_loss)
    
    vae.add_loss(vae_loss)
    vae.compile(optimizer='adam')
    
    return vae, encoder, decoder

def compare_models():
    """比较不同深度学习模型"""
    print("\n9. 模型性能比较...")
    
    # 这里可以加载之前训练的模型进行比较
    # 或者记录训练过程中的性能指标
    
    models_performance = {
        '标准神经网络': 0.85,  # 示例数据
        '自编码器+分类器': 0.87,
        '变分自编码器': 0.83
    }
    
    plt.figure(figsize=(10, 6))
    models = list(models_performance.keys())
    performances = list(models_performance.values())
    
    bars = plt.bar(models, performances, color=['skyblue', 'lightgreen', 'lightcoral'])
    plt.title('深度学习模型性能比较')
    plt.ylabel('准确率')
    plt.ylim(0, 1)
    
    # 添加数值标签
    for bar, perf in zip(bars, performances):
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01, 
                f'{perf:.3f}', ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig('plots/dl_models_comparison.pdf', dpi=300, bbox_inches='tight')
    plt.show()

def generate_report(test_acc, feature_names):
    """生成深度学习分析报告"""
    print("\n10. 生成分析报告...")
    
    with open("results/deep_learning_report.txt", "w", encoding='utf-8') as f:
        f.write("=== 深度学习分析报告 ===\n")
        f.write(f"分析时间: {pd.Timestamp.now()}\n\n")
        
        f.write("1. 模型架构:\n")
        f.write("   - 标准神经网络分类器\n")
        f.write("   - 自编码器特征提取\n")
        f.write("   - 变分自编码器（可选）\n\n")
        
        f.write("2. 数据概况:\n")
        f.write(f"   - 特征数量: {len(feature_names)}\n")
        f.write("   - 数据预处理: 标准化、缺失值处理\n\n")
        
        f.write("3. 模型性能:\n")
        f.write(f"   - 测试准确率: {test_acc:.4f}\n")
        f.write("   - 使用早停和学习率调度\n")
        f.write("   - 批量归一化和Dropout正则化\n\n")
        
        f.write("4. 技术特点:\n")
        f.write("   - 多层全连接网络\n")
        f.write("   - 自适应学习率\n")
        f.write("   - 特征编码和可视化\n\n")
        
        f.write("5. 输出文件:\n")
        f.write("   - 模型文件: models/\n")
        f.write("   - 可视化图表: plots/\n")
        f.write("   - 分析报告: results/deep_learning_report.txt\n\n")
        
        f.write("分析完成！\n")
    
    print("  深度学习分析报告已保存到 results/deep_learning_report.txt")

def main():
    """主函数"""
    # 设置环境
    if not setup_environment():
        return
    
    # 加载数据
    X, y, le, feature_names = load_and_prepare_data()
    if X is None:
        return
    
    # 数据分割
    X_train, X_temp, y_train, y_temp = train_test_split(
        X, y, test_size=0.4, random_state=42, stratify=y
    )
    X_val, X_test, y_val, y_test = train_test_split(
        X_temp, y_temp, test_size=0.5, random_state=42, stratify=y_temp
    )
    
    print(f"训练集: {X_train.shape}, 验证集: {X_val.shape}, 测试集: {X_test.shape}")
    
    # 创建和训练自编码器
    autoencoder, encoder = create_autoencoder(X.shape[1], encoding_dim=64)
    ae_history = train_autoencoder(autoencoder, X_train, X_val, epochs=50)
    
    # 可视化自编码器训练历史
    visualize_training_history(ae_history, "自编码器", "plots/autoencoder_history.pdf")
    
    # 使用编码器提取特征
    X_train_encoded = encoder.predict(X_train)
    X_val_encoded = encoder.predict(X_val)
    X_test_encoded = encoder.predict(X_test)
    
    # 创建和训练分类器（使用编码特征）
    classifier = create_classifier(X_train_encoded.shape[1], len(le.classes_))
    clf_history = train_classifier(classifier, X_train_encoded, X_val_encoded, 
                                 y_train, y_val, epochs=100)
    
    # 可视化分类器训练历史
    visualize_training_history(clf_history, "分类器", "plots/classifier_history.pdf")
    
    # 评估模型
    test_acc, y_pred, y_pred_proba = evaluate_classifier(classifier, X_test_encoded, y_test, le)
    
    # 可视化编码特征
    visualize_encoded_features(encoder, X_test, y_test, le)
    
    # 创建和训练VAE（可选）
    try:
        vae, vae_encoder, vae_decoder = create_variational_autoencoder(X.shape[1], latent_dim=32)
        vae_history = vae.fit(X_train, epochs=30, batch_size=16, 
                             validation_data=(X_val, None), verbose=1)
        visualize_training_history(vae_history, "变分自编码器", "plots/vae_history.pdf")
    except Exception as e:
        print(f"VAE训练出错: {e}")
    
    # 模型比较
    compare_models()
    
    # 保存模型
    autoencoder.save('models/autoencoder.h5')
    encoder.save('models/encoder.h5')
    classifier.save('models/classifier.h5')
    print("模型已保存到 models/ 目录")
    
    # 生成报告
    generate_report(test_acc, feature_names)
    
    print("\n=== 深度学习分析完成 ===")
    print("输出文件:")
    print("  - 模型文件: models/")
    print("  - 可视化图表: plots/")
    print("  - 分析报告: results/deep_learning_report.txt")

if __name__ == "__main__":
    main()