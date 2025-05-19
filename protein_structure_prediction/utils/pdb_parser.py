#! /user/bin/python3.7
# coding: utf-8
# @Time: 2025-04-19
# @Author: 

"""
蛋白质CIF/PDB文件解析工具
功能：将蛋白质结构文件转换为图数据（节点=氨基酸，边=空间临近关系）
用法：python 
"""

import os
import numpy as np 

import torch

# PyTorch Geometric是一个用于图神经网络的库
# torch_geometric.data.Data 是 PyTorch Geometric（PyG） 库中定义图数据结构（单张图）的核心类，
# 用于表示和存储图数据（节点、边、特征等）。它是构建图神经网络（GNN）模型的基础，
# 几乎所有 PyG 的模型和数据处理工具都依赖此模块。
from torch_geometric.data import Data

# Data 类的主要属性
# |属性名|作用|示例值类型|
# |-----|----|--------|
# |x|节点特征矩阵（形状为 [num_nodes, num_node_features]）|torch.FloatTensor|
# |edge_index|边的连接关系（形状为 [2, num_edges]，存储边的起点和终点索引）|torch.LongTensor|
# |edge_attr|边特征矩阵（形状为 [num_edges, num_edge_features]）|torch.FloatTensor|
# |y|标签（可以是节点级、边级或图级标签）|torch.Tensor|
# |pos|节点的坐标（常用于图结构数据，如分子结构、点云）|torch.FloatTensor|
# |batch|批处理索引（用于将多个图合并成一个批量）|torch.LongTensor|

from Bio.PDB import MMCIFParser,PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning

import warnings   # 管理和控制程序运行过程中产生的警告信息

# Union类型提示系统中实现联合类型的核心工具
from typing import Union # 类型注解（Type Hints）的工具

# 当Biopython解析PDB文件遇到非关键问题（如格式轻微不符或缺失数据）时，
# 会抛出PDBConstructionWarning。这类警告可能频繁出现，尤其是在处理大量文件时，禁用它们可使输出更简洁。

# 禁用Biopython库在处理PDB文件时产生的特定警告（PDBConstructionWarning）
warnings.filterwarnings("ignore", category=PDBConstructionWarning)



# ==================常量定义===========================

# one-hot编码是一种将离散型分类变量转换为二进制向量表示的常用方法
# 独热性：每个氨基酸对应一个唯一的20维向量（因为有20种标准氨基酸）
# 二进制表示：向量中只有1个元素为1（表示存在），其余19个元素均为0

# 20种标准氨基酸的one-hot编码字典
AMINO_ACID_IDS = {
    'ALA': 0, 'ARG': 1, 'ASN': 2, 'ASP': 3, 'CYS': 4,
    'GLN': 5, 'GLU': 6, 'GLY': 7, 'HIS': 8, 'ILE': 9,
    'LEU': 10, 'LYS': 11, 'MET': 12, 'PHE': 13, 'PRO': 14,
    'SER': 15, 'THR': 17, 'TRP': 17, 'TYR': 18, 'VAL': 19
}


# ==================主解析函数===========================
def parser_cif_graph(
    cif_path: str,
    k_neighbors: int = 10,
    max_distance: float = 10.0,
    use_calpha: bool = True
) -> Union[Data, None]:
    """
    将CIF文件解析为PyTorch Geometric的图数据
    参数：
       cif_path（str）: CIF文件路径
       k_neighbors(int): 每个节点的最近邻数量
       max_distance(float): 边连接的最大距离阈值(Å)
       use_calpha(bool):是否仅使用Cα原子（True）或所有原子（False）
    返回:
        Data: PyG图数据对象，包含节点特征、边索引、边属性、坐标标签
    """
    try:
        # ---------------1、解析CIF文件----------------------------------
        # 创建解析mmCIF格式文件的解析器对象
        # QUIET=True：抑制解析过程中的警告和错误输出，保持静默
        parser = MMCIFParser(QUIET = True) 
        
        # 从指定路径的CIF文件中解析蛋白质结构数据
        # "protein"：为结构赋予的标识符（可自定义）
        # 返回值：structure是一个层级结构对象，
        # 包含模型（Model）、链（Chain）、残基（Residue）、原子（Atom）的嵌套信息
        structure = parser.get_structure("protein", cif_path)
        
        # 提取所有Cα原子
        atoms = []
        # 遍历层级结构提取原子
        # 结构层级：
        # Model：代表不同的构象（如NMR结构中的多个模型，X射线通常只有一个）。
        # Chain：蛋白质的一条多肽链或核酸链（如链A、B）。
        # Residue：单个氨基酸或核苷酸残基（如第10号丙氨酸）。
        for model in structure:
            for chain in model:
                for residue in chain:
                    if use_calpha:
                        # "CA"是Cα原子的标准名称
                        if "CA" in residue:
                            atoms.append(residue['CA'])
                        else:
                            # residue.get_atoms() 返回一个生成器（Generator），遍历该残基（Residue）中包含的所有原子（Atom）
                            # atoms.extend(...) 将可迭代对象（如列表、生成器）中的元素逐个追加到列表末尾
                            atoms.extend(residue.get_atoms())

        if len(atoms) == 0:
            raise ValueError("没有找到Cα原子！检查文件格式，或者使用use_calpha=False")

        # -------------------- 2. 提取坐标和特征 --------------------
        # 提取每个原子的三维坐标 
        # 最终得到一个形状为 (N, 3) 的 NumPy 数组，其中 N 是原子数量
        corrds = np.array([atom.get_coord() for atom in atoms]) 
        # 获取残基名
        # 假设所有原子属于同一个链atoms[0].get_parent().parent 获取链对象
        # atoms[0].get_parent()应该得到该原子所属的残基对象
        residue_names = [residue.resname for residue in atoms[0].get_parent().parent]

        # one-hot 氨基酸类型
        
        # 创建一个全零张量，形状为 (残基数量, 氨基酸类别数)
        node_feat = torch.zeros(len(residue_names), len(AMINO_ACID_IDS))
        print(node_feat)
        # 对每个残基名，若存在于预设的 AMINO_ACID_IDS 字典（如 {"ALA": 0, "ARG": 1, ...}），
        # 则在对应位置赋值为 1，实现 One-Hot 编码
        for i, resname in enumerate(residue_names):
            if resname in AMINO_ACID_IDS:
                node_feat[i, AMINO_ACID_IDS[resname]] = 1.0
            else:
                # 处理非标准氨基酸（用全零向量表示）
                # 或者在AMINO_ACID_IDS字典中添加未知类别
                pass
        

        # -------------------- 3. 构建图连接关系 --------------------
    except Exception as e:
        print(f"解析文件{os.path.basename(cif_path)}失败：{str(e)}")
        return None


# ==================批量处理与保存===========================


# ==================示例用法===========================
if __name__ == '__main__':
    # 单文件测试
    cif_path = './protein_structure_prediction/data/PDB/1TUP.cif'
    data = parser_cif_graph(cif_path, k_neighbors=15)
    print( data)
    # print("节点数：", data.num_nodes)
    # print("边数：", data.edge_index.shape[1])

    