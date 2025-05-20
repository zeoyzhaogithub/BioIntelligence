#! /user/bin/python3.7
# coding: utf-8
# @Time: 2025-04-19
# @Author: zeoyzhao

"""
蛋白质结构预测训练入口
用法：python train.py --config configs/train_config.yaml
"""

import argparse  # 为脚本定义、解析和处理命令行参数，自动生成帮助信息
import yaml      # 解析和生成 YAML 格式数据
import torch     # 开源的深度学习框架，专注于提供灵活的张量计算和神经网络构建功能

# 高效加载数据，将数据集包装成一个可迭代的对象，支持批量加载、多进程加速、数据打乱等操作
from torch.utils.data import DataLoader 


from utils.data_loader import ProteinDataset

def main(config):
    # -------------一、初始化系统配置--------------
    # 1.设备设置

    # 当前mac不支持cuda
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    print(f"使用设备:{device}")
    
    # -------------二、数据集加载------------------
    # 训练数据集
    train_dataset = ProteinDataset(root=config['data']['train_dir'])

    # 验证数据集
    val_dataset = ProteinDataset(root=config['data']['val_dir'])
    



if __name__ == '__main__':
    print('hello protein structure prediction')

    # 命令行参数解析
    parser = argparse.ArgumentParser()  # 为脚本定义、解析和处理命令行参数，自动生成帮助信息
    # parser.add_argument("--config", type=str,required=True,help="path to config file")
    parser.add_argument("--config", type=str,help="path to config file",default='./protein_structure_prediction/configs/params.yaml')
    
    args = parser.parse_args()
    # config_path = './protein_structure_prediction/configs/params.yaml' if args.config == '' else args.config 
    # print(args)
    # 加载配置文件
    with open(args.config, "r") as f:
        config = yaml.safe_load(f)
    
    main(config)

