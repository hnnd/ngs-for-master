---
marp: true
theme: ngs-course
paginate: false
header: '高通量测序数据分析'
footer: '王运生 | 2025'
---

<!-- 
Linux系统基础课程幻灯片
课程名称：高通量测序数据分析 - 预备课程
主讲教师：王运生教授
联系邮箱：wangys@hunau.edu.cn
办公室：16教420室
上课地点：105机房
-->

<!-- _class: title -->
# Linux系统基础
## 高通量测序数据分析 - 预备课程

**主讲教师：** 王运生  
**联系邮箱：** wangys@hunau.edu.cn  
**办公室：** 16教420室  
**上课地点：** 105机房  

---

<!-- _class: toc -->
# 本次课程内容

1. Linux系统概述与架构
2. 命令行基础操作
3. 文件系统与目录管理
4. 文本处理工具
5. 进程管理与系统监控

---

<!-- _class: content -->
# Linux系统概述
**学习目标：**
- 理解Linux系统架构和特点
- 掌握基本的命令行操作
- 熟练使用文件管理命令
- 学会使用文本处理工具

---

<!-- _class: content -->
# Linux系统概述

<div class="columns">
<div class="column">

## 什么是Linux？

- **开源操作系统**：基于Unix的自由软件
- **多用户多任务**：支持多个用户同时使用
- **稳定可靠**：广泛应用于服务器和科学计算
- **命令行界面**：强大的文本界面操作

</div>
<div class="column">

## 为什么生物信息学需要Linux？

- 大多数生物信息学软件运行在Linux上
- 强大的文本处理能力
- 高效的批处理和自动化
- 丰富的开源工具生态

</div>
</div>

---

<!-- _class: multi-column -->
# Linux发行版比较

<div class="columns">
<div class="column">

## 常见发行版
- **Ubuntu**：用户友好，适合初学者
- **CentOS/RHEL**：企业级，稳定性好
- **Debian**：纯开源，包管理优秀
- **Fedora**：新技术，更新频繁

</div>
<div class="column">

## 生物信息学推荐
- **Ubuntu LTS**：长期支持版本
- **CentOS**：服务器环境常用
- **Bio-Linux**：专门的生物信息学发行版
- **Docker容器**：轻量级虚拟化

</div>
</div>

---

<!-- _class: content -->
# Linux系统架构

<div class="columns">
<div class="column">

## 系统层次结构

- **内核（Kernel）**：系统核心，管理硬件资源
- **Shell**：命令解释器，用户与系统交互界面
- **文件系统**：数据存储和组织方式
- **应用程序**：用户使用的各种软件

</div>
<div class="column">

## 文件系统层次标准（FHS）

- `/`：根目录，所有文件的起点
- `/home`：用户主目录
- `/usr`：用户程序和数据
- `/var`：变量数据文件
- `/tmp`：临时文件

</div>
</div>

---

<!-- _class: content -->
# 命令行基础

## Shell简介

- **Bash**：最常用的Shell（Bourne Again Shell）
- **命令提示符**：显示当前用户和路径信息
- **命令格式**：`命令 [选项] [参数]`

## 基本命令结构

```bash
# 命令的基本格式
command -option argument

# 示例
ls -l /home/user
```

---

<!-- _class: code -->
# 基本导航命令

```bash
# 显示当前目录
pwd

# 列出文件和目录
ls
ls -l          # 详细信息
ls -la         # 包含隐藏文件
ls -lh         # 人类可读的文件大小

# 切换目录
cd /path/to/directory
cd ~           # 回到主目录
cd ..          # 上级目录
cd -           # 上次访问的目录
```

---

<!-- _class: code -->
# 文件和目录操作

<div class="columns">
<div class="column">

```bash
# 显示目录树结构
tree
tree -L 2      # 只显示2级深度
```

```bash
# 创建目录
mkdir directory_name
mkdir -p path/to/directory    # 创建多级目录

# 创建文件
touch filename.txt
echo "content" > filename.txt
```

</div>
<div class="column">

## 复制和移动

```bash
# 复制文件
cp source.txt destination.txt
cp -r source_dir dest_dir     # 复制目录

# 移动/重命名
mv old_name.txt new_name.txt
mv file.txt /path/to/directory/
```

---

<!-- _class: code -->
# 文件查看和编辑

```bash
# 查看文件内容
cat filename.txt              # 显示全部内容
less filename.txt             # 分页查看
head filename.txt             # 前10行
head -n 20 filename.txt       # 前20行
tail filename.txt             # 后10行
tail -f logfile.txt           # 实时查看日志

# 文件编辑器
nano filename.txt             # 简单编辑器
vim filename.txt              # 高级编辑器
```

---

<!-- _class: content -->
# 文件权限管理

<div class="columns">
<div class="column">

## 权限概念

- **用户类型**：所有者(u)、组(g)、其他(o)
- **权限类型**：读(r)、写(w)、执行(x)
- **权限表示**：数字(755)或字符(rwxr-xr-x)

</div>
<div class="column">

## 权限操作

```bash
# 查看权限
ls -l filename.txt

# 修改权限
chmod 755 script.sh           # 数字方式
chmod u+x script.sh           # 字符方式
chmod -R 644 directory/       # 递归修改

# 修改所有者
chown user:group filename.txt
```

</div>
</div>

---

<!-- _class: multi-column -->
# 文本处理工具

<div class="columns">
<div class="column">

## 基础工具
```bash
# 搜索文本
grep "pattern" file.txt
grep -i "pattern" file.txt    # 忽略大小写
grep -r "pattern" directory/  # 递归搜索

# 排序和去重
sort file.txt
sort -n numbers.txt           # 数字排序
uniq file.txt                 # 去重
sort file.txt | uniq          # 排序后去重
```

</div>
<div class="column">

## 高级处理
```bash
# 文本统计
wc file.txt                   # 行数、词数、字符数
wc -l file.txt                # 只统计行数

# 文本替换
sed 's/old/new/g' file.txt    # 替换文本
sed -i 's/old/new/g' file.txt # 直接修改文件

# 字段处理
cut -f1,3 -d',' file.csv      # 提取CSV字段
awk '{print $1, $3}' file.txt # 打印特定列
```

</div>
</div>

---

<!-- _class: content -->
# 管道和重定向

<div class="columns">
<div class="column">

## 重定向操作

```bash
# 输出重定向
command > output.txt          # 覆盖写入
command >> output.txt         # 追加写入
command 2> error.log          # 错误重定向
command > output.txt 2>&1     # 合并输出和错误

# 输入重定向
command < input.txt
```

</div>
<div class="column">

## 管道操作

```bash
# 管道连接命令
ls -l | grep ".txt"           # 列出txt文件
cat file.txt | sort | uniq    # 排序去重
ps aux | grep python          # 查找Python进程
```

</div>
</div>

---

<!-- _class: code -->
# 进程管理

```bash
# 查看进程
ps                            # 当前终端进程
ps aux                        # 所有进程详细信息
ps aux | grep process_name    # 查找特定进程

# 进程控制
command &                     # 后台运行
nohup command &               # 后台运行，忽略挂断信号
jobs                          # 查看后台任务
fg %1                         # 将后台任务调到前台
bg %1                         # 将暂停任务放到后台

# 终止进程
kill PID                      # 终止进程
kill -9 PID                   # 强制终止
killall process_name          # 按名称终止
```

---

<!-- _class: content -->
# 系统监控

## 系统资源监控

```bash
# CPU和内存使用
top                           # 实时系统监控
htop                          # 增强版top
free -h                       # 内存使用情况
df -h                         # 磁盘使用情况
du -sh directory/             # 目录大小

# 网络监控
netstat -tuln                 # 网络连接状态
ss -tuln                      # 现代版netstat
```

---

<!-- _class: content -->
# 环境变量和配置

## 环境变量

```bash
# 查看环境变量
env                           # 所有环境变量
echo $PATH                    # 查看PATH变量
echo $HOME                    # 用户主目录

# 设置环境变量
export VAR_NAME=value         # 临时设置
echo 'export VAR_NAME=value' >> ~/.bashrc  # 永久设置
source ~/.bashrc              # 重新加载配置
```

---

<!-- _class: code -->
## 配置文件

- `~/.bashrc`：Bash配置文件
- `~/.profile`：登录时执行的配置
- `/etc/environment`：系统级环境变量

---

<!-- _class: content -->
# 软件包管理

## Ubuntu/Debian系统

```bash
# 更新包列表
sudo apt update

# 安装软件
sudo apt install package_name
sudo apt install -y package_name      # 自动确认

# 搜索软件包
apt search keyword

# 卸载软件
sudo apt remove package_name
sudo apt purge package_name           # 完全删除包括配置
```

---

<!-- _class: content -->
# 压缩和解压

## 常用压缩格式

```bash
# tar格式
tar -czf archive.tar.gz directory/    # 压缩
tar -xzf archive.tar.gz               # 解压

# zip格式
zip -r archive.zip directory/         # 压缩
unzip archive.zip                     # 解压

# gzip格式
gzip file.txt                         # 压缩单个文件
gunzip file.txt.gz                    # 解压
```

---

<!-- _class: content -->
# 网络操作

## 文件传输

```bash
# 下载文件
wget http://example.com/file.txt
curl -O http://example.com/file.txt

# 远程连接
ssh user@hostname                     # SSH连接
scp file.txt user@hostname:/path/     # 复制文件到远程
rsync -av local/ user@hostname:remote/ # 同步目录
```

---

<!-- _class: code -->
# 网络监控

## 网络诊断

```bash
# 网络连通性
ping hostname
traceroute hostname
nslookup hostname
```

---

<!-- _class: content -->
# 脚本编写基础

## Bash脚本结构

```bash
#!/bin/bash
# 脚本说明

# 变量定义
VAR_NAME="value"

# 条件判断
if [ condition ]; then
    echo "条件为真"
else
    echo "条件为假"
fi

# 循环
for file in *.txt; do
    echo "处理文件: $file"
done
```

---

<!-- _class: code -->
# 实用脚本示例

```bash
#!/bin/bash
# 批量重命名文件

# 将所有.txt文件重命名为.bak
for file in *.txt; do
    if [ -f "$file" ]; then
        mv "$file" "${file%.txt}.bak"
        echo "重命名: $file -> ${file%.txt}.bak"
    fi
done

# 批量创建目录
mkdir -p data/{raw,processed,results}
mkdir -p scripts/{analysis,visualization}
mkdir -p logs

echo "目录结构创建完成"
```

---

<!-- _class: content -->
# 生物信息学中的Linux应用

## 常见使用场景

- **数据预处理**：质量控制、格式转换
- **序列分析**：比对、组装、注释
- **批量处理**：处理大量样本数据
- **流程自动化**：编写分析流水线

---

## 重要技能

- 熟练使用命令行操作大文件
- 编写脚本自动化重复任务
- 理解文件格式和数据结构
- 掌握远程服务器操作

---

<!-- _class: content -->
# 常用快捷键

## 命令行快捷键

- `Ctrl + C`：终止当前命令
- `Ctrl + Z`：暂停当前命令
- `Ctrl + L`：清屏
- `Ctrl + A`：光标移到行首
- `Ctrl + E`：光标移到行尾
- `Ctrl + R`：搜索历史命令
- `Tab`：自动补全
- `↑/↓`：浏览历史命令

---

<!-- _class: summary -->
# 本次课程总结

## 主要内容回顾
- Linux系统架构和特点
- 基本命令行操作和文件管理
- 文本处理工具的使用
- 进程管理和系统监控
- 脚本编写基础

## **作业/练习：**
- 完成Linux基础操作练习
- 编写简单的文件处理脚本

---

<!-- _class: end -->
# 谢谢大家！

**有问题请联系：**
- 邮箱：wangys@hunau.edu.cn
- 办公室：16教420室
