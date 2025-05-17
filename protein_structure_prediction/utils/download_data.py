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

# import Bio
# print("Biopython 版本:", Bio.__version__)  # Biopython 版本: 1.79


# 检查服务器连接
# def check_server():
#     try:
#         socket.create_connection(("files.rcsb.org", 80), 5)
#         return True
#     except OSError:
#         return False

#     if not check_server():
#         raise ConnectionError("无法连接至 PDB 服务器")
    


def download_data(pdb_id,p_dir):
    os.makedirs(p_dir,exist_ok=True)   # 当目标目录 p_dir 已经存在时，不引发错误，而是静默忽略
    pdbl = PDBList()
    # 检查 PDBList 支持的格式列表
    # print(pdbl.__format__)
    # print(pdbl.download_pdb_files)
    # 打印 PDBList 类的属性和方法
#     print(dir(PDBList))  # 查找类似 'format_map' 或 'file_format' 的字段
# ['PDB_REF', '__class__', '__delattr__', '__dict__', '__dir__', '__doc__', 
# '__eq__', '__format__', '__ge__', '__getattribute__', '__gt__', '__hash__',
# '__init__', '__init_subclass__', '__le__', '__lt__', '__module__', '__ne__', 
# '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', 
# '__str__', '__subclasshook__', '__weakref__', '_print_default_format_warning', 
# 'download_entire_pdb', 'download_obsolete_entries', 'download_pdb_files', 
# 'get_all_entries', 'get_all_obsolete', 'get_recent_changes', 'get_seqres_file', 
# 'get_status_list', 'retrieve_pdb_file', 'update_pdb']
    # 方法一：使用Biopython下载
    retrieve = pdbl.retrieve_pdb_file(pbd_id, file_format='cif', pdir=p_dir)  # 超小蛋白（胰蛋白酶抑制剂）
    try:
        retrieve = pdbl.retrieve_pdb_file(
            pdb_id, 
            file_format='pdb', 
            pdir=p_dir)  # 超小蛋白（胰蛋白酶抑制剂）
        
        # 这里下载失败，但是并没有抛出异常，需要手动检查
        if not os.path.exists(retrieve):
            raise FileNotFoundError(f'f"文件未生成: {cif_path}"')
        return retrieve
    except (Exception, FileNotFoundError) as e:
        print(f'下载失败,{str(e)}')
        # 备用方案：直接通过url下载
        url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        resp = requests.get(url, timeout=10)
        resp.raise_for_status()  # 自动触发 HTTPError

        cif_path = os.path.join(p_dir, f'{pdb_id}.cif')
        if resp.status_code == 200:
            with open(cif_path, 'wb') as f:
                f.write(resp.content)
                print('手动下载成功')
                return cif_path
        else:
            print("服务器返回错误：文件不存在！")

if __name__ == '__main__':
    # 蛋白质id
    pdb_id = '1TUP'
    # 存储路径
    p_dir = './protein_structure_prediction/data/PDB/'
    try:
        path = download_data(pdb_id, p_dir)
        print(f'最终下载路径：{path}')
    except Exception as e:
        print(f'全局错误：{str(e)}')
        