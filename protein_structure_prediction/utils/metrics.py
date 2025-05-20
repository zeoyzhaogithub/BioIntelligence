#! /user/bin/python3.7
# coding: utf-8
# @Time: 2025-04-19
# @Author: 

"""
生物特异性评价指标
用法：python 
"""

import torch
from torch.utils.data import Dataset,Dataloader

from Bio.PDB import PDBParser, DSSP
import numpy as np
import os

from torch_geometric.data import Batch,Data