import os
import numpy as np
from scipy.io import loadmat, savemat
from concurrent.futures import ThreadPoolExecutor

# 加载数据
frequency_k = loadmat('../data/processing/k.mat')
k = frequency_k['k']
omega_dict = loadmat('../data/source/omega.mat')
omega_1024 = omega_dict['omega_1024']
omega_512 = omega_1024.squeeze()

numRows, numCols = k.shape
blockSize = 256
numBlocks = numRows // blockSize

# 确保保存文件的文件夹存在
os.makedirs('../data/processing/per', exist_ok=True)

def process_block(b):
    rowStart = (b - 1) * blockSize
    rowEnd = b * blockSize
    k_block = k[rowStart:rowEnd, :]
    omega_e_block = np.einsum('ij,k->ijk', k_block, omega_512)
    per = np.round(2 * np.pi / omega_e_block, decimals=1)
    per = per.reshape((blockSize, numCols, -1))
    filename = f'../data/processing/per/per_{rowStart}-{rowEnd-1}.mat'
    per_to_save = np.ascontiguousarray(per)
    savemat(filename, {'per': per_to_save}, format='5')
    print(f'Saved {filename} with shape {per_to_save.shape}')

with ThreadPoolExecutor(max_workers=4) as executor:
    executor.map(process_block, range(1, numBlocks + 1))
