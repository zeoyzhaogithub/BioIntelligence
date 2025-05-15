#! /user/bin/python3.7
# coding: utf-8
# @Time: 2025-04-21
# @Author: 

"""
数据下载
用法：python 
"""

from Bio.PDB import PDBList
import os
import requests
import socket     # 实现不同设备间进程通信

import Bio
print("Biopython 版本:", Bio.__version__)  # Biopython 版本: 1.79
  


def download_data(pdb_id,p_dir):
    os.makedirs(p_dir,exist_ok=True)   # 当目标目录 p_dir 已经存在时，不引发错误，而是静默忽略
    pdbl = PDBList()


if __name__ == '__main__':
    # 蛋白质id
    pdb_id = '1TUP'
