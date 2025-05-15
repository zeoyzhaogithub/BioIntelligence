# 基于深度学习的生物大分子结构与功能预测系统demo


## 项目结构
protein_structure_prediction/
├── data/                     # 数据集（需预处理）
│   ├── PDB/                  # 蛋白质结构数据（来自RCSB PDB）
│   ├── DNA_protein_pairs/    # DNA/RNA-蛋白质结合数据（来自BioLiP）
│   ├── raw/                  # 原始数据
│   ├── processed/            # 预处理数据
│   └── splits/               # 数据集划分
├── models/
│   ├── protein_structure/    # 蛋白质结构预测模型
│   ├── dna_binding/          # DNA-蛋白质结合预测模型
│   ├── transformer_models.py # Transformer模块
│   └── ensemble.py           # 模型集成
├── train.py                  # 训练主程序
├── configs/
│   └── params.yaml           # 超参数配置
├── utils/
│   ├── data_loader.py        # 数据加载
│   ├── metrics.py            # 评估指标
│   ├── pdb_parser.py         # PDB文件解析工具
│   ├── data_loader.py        # 数据加载与特征工程
│   └── metrics.py            # 生物特异性评估指标
│   └── visualization.py      # 结果可视化
└── results/
│   ├── checkpoints/          # 模型保存
│   └── predictions/          # 预测结果
├── train_structure.py      # 蛋白质结构训练脚本
├── train_binding.py        # 结合亲和力训练脚本
├── requirements.txt        # 依赖库（包含Biopython、PyTorch Geometric）
└── README.md               # 项目说明（含可视化示例）

### 安装 PyTorch Geometric 依赖库
根据 PyTorch 和 CUDA 版本替换 ${CUDA} 和 ${TORCH}：

${CUDA}：cu117（CUDA 11.7）、cu118（CUDA 11.8）或 cpu。

${TORCH}：如 2.0.0。

```bash
pip install torch_scatter torch_sparse torch_cluster torch_spline_conv -f https://data.pyg.org/whl/torch-${TORCH}+${CUDA}.html

pip install torch_scatter torch_sparse torch_cluster torch_spline_conv -f https://data.pyg.org/whl/torch-1.13.1+cpu.html
```

4. 安装 torch_geometric
```bash
pip install torch_geometric
```