# Vector-Borne Disease Modelling: SEIR-SEI & Rt Estimation (R)

This project provides a modular R implementation for simulating and analyzing the transmission dynamics of vector-borne diseases（如登革热、基孔肯雅热）using a coupled SEIR model for humans and SEI model for vectors, and for estimating the real-time reproduction number (Rt) using EpiEstim.

## 文件结构

- `seir_dengue_model.R`         — 主脚本，包含模型调用与示例
- `model_seir_SEI.R`            — SEIR-SEI动力学模型定义
- `fit_incubation_distribution.R` — 潜伏期分布拟合与参数提取
- `estimate_rt_epiestim.R`      — 用EpiEstim估算Rt
- `requirements.txt`            — 依赖R包列表
- `README.md`                   — 项目说明文档

## 快速开始

1. 安装依赖包：

   ```r
   install.packages(c('deSolve', 'ggplot2', 'reshape2', 'MASS', 'EpiEstim'))
   ```

2. 拟合潜伏期分布（如有原始数据）：

   ```r
   source('fit_incubation_distribution.R')
   # 输出均值、标准差、shape、rate等参数
   ```

3. 运行动力学模型：

   ```r
   source('model_seir_SEI.R')
   # 生成人群和蚊群的SEIR-SEI仿真结果
   ```

4. 估算实时再生数Rt：

   ```r
   source('estimate_rt_epiestim.R')
   # 生成Rt曲线
   ```

## 主要功能说明

### 1. 潜伏期分布拟合

- 支持伽马分布、对数正态分布等
- 自动输出均值、标准差、shape、rate等参数
- 结果可直接用于动力学模型参数设定

### 2. SEIR-SEI动力学模型

- 人群：SEIR结构
- 蚊群：SEI结构
- 参数可自定义，支持与分布拟合结果联动
- 输出各舱室人数随时间变化曲线

### 3. 实时再生数Rt估算

- 支持用模型输出或实际病例数据
- 采用EpiEstim包，支持自定义代间分布
- 输出Rt随时间变化曲线

## 参考资料

- [EpiEstim R包](https://cran.r-project.org/web/packages/EpiEstim/index.html)
- [deSolve R包](https://cran.r-project.org/web/packages/deSolve/index.html)
- [登革热/基孔肯雅热流行病学文献]

## 联系方式

如有问题或建议，请联系项目维护者。
